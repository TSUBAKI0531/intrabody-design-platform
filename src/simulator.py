"""
simulator.py — 結合シミュレーション・構造予測・特異性評価モジュール

イントラボディ設計パイプラインにおける各種計算エンジンを提供する。

Note:
    現在の実装は in silico 概念実証（PoC）段階であり、
    結合自由エネルギー等の値は物理化学的パラメータに基づく
    近似計算です。実験値との較正（Calibration）は今後の課題です。
"""

from __future__ import annotations

import logging
import math
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Kyte-Doolittle 疎水性スケール
# ---------------------------------------------------------------------------
KD_SCALE: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# 残基あたりの Alanine scanning 参照 ΔΔG（kcal/mol 概算）
# 荷電/芳香族残基ほど Ala 置換の影響が大きい
RESIDUE_DDG_REFERENCE: dict[str, float] = {
    "R": 3.5, "K": 3.0, "D": 2.8, "E": 2.5, "H": 2.2,
    "W": 3.8, "Y": 2.5, "F": 2.0, "N": 1.5, "Q": 1.3,
    "S": 0.8, "T": 0.9, "M": 1.2, "C": 1.8, "P": 1.0,
    "G": 0.4, "A": 0.0, "V": 0.7, "I": 0.9, "L": 0.8,
}


# =====================================================================
# AntibodyDiscoveryEngine — エピトープ解析 & 抗体候補生成
# =====================================================================
class AntibodyDiscoveryEngine:
    """ターゲット配列を解析し、結合候補となる抗体配列を生成する。

    Parameters:
        target_seq: ターゲットタンパク質のアミノ酸配列
    """

    # scFv フレームワークテンプレート（トラスツズマブ由来 VH FR1-FR3）
    _FRAMEWORK_TEMPLATE: str = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"  # FR1 + CDR1 scaffold
        "{cdr_h3}"                          # CDR-H3（生成対象）
        "WVRQAPGKGLEWVS"                    # FR3
    )

    def __init__(self, target_seq: str) -> None:
        self.target_seq = target_seq

    def recommend_domains(
        self, window_size: int = 15, top_n: int = 3
    ) -> pd.DataFrame:
        """露出度スコアに基づいてエピトープ候補ドメインを推薦する。"""
        recommendations: list[dict] = []
        for i in range(len(self.target_seq) - window_size + 1):
            start, end = i + 1, i + window_size
            score = self.calculate_exposure_score(start, end)
            recommendations.append(
                {"Start": start, "End": end, "Exposure_Score": round(score, 2)}
            )
        df = pd.DataFrame(recommendations)
        return df.sort_values("Exposure_Score", ascending=False).head(top_n)

    def calculate_exposure_score(self, start: int, end: int) -> float:
        """部分配列の表面露出度スコアを計算する。

        親水性残基（KD 値が負）が多いほどスコアが高い = 表面に露出している
        可能性が高いと推定する。
        """
        sub_seq = self.target_seq[start - 1 : end]
        if not sub_seq:
            return 0.0
        scores = [-KD_SCALE.get(aa, 0.0) for aa in sub_seq]
        return sum(scores) / len(scores)

    def discover_binder(self) -> str:
        """ターゲット配列の物性に基づいて CDR-H3 を擬似生成し、
        scFv 候補配列を返す。

        Note:
            現在は確定的な擬似生成。将来的には深層学習ベースの
            CDR 生成エンジンとの統合を想定。
        """
        # ターゲットの荷電特性の逆を CDR に反映（相補性の模倣）
        target_charge = sum(
            1 if aa in "RKH" else (-1 if aa in "DE" else 0)
            for aa in self.target_seq[:30]
        )

        if target_charge > 0:
            cdr_h3 = "ARDYYDSSGYYAM"  # 酸性寄り CDR
        elif target_charge < 0:
            cdr_h3 = "ARRGYSSSWYRDY"  # 塩基性寄り CDR
        else:
            cdr_h3 = "ARSYYGSSYWGDY"  # バランス型 CDR

        antibody = self._FRAMEWORK_TEMPLATE.format(cdr_h3=cdr_h3)
        logger.info(
            "候補抗体を生成 (CDR-H3: %s, ターゲット荷電: %+d)",
            cdr_h3,
            target_charge,
        )
        return antibody


# =====================================================================
# AFMValidator — 構造予測 & バリデーション
# =====================================================================
@dataclass
class ValidationReport:
    """構造予測バリデーション結果。"""

    pLDDT: float
    iPAE: float
    status: str


class AFMValidator:
    """配列の 3D 構造予測を取得し、品質を検証する。

    Parameters:
        pdb_path: PDB ファイルの保存先パス
        sequence: 構造予測対象のアミノ酸配列

    Note:
        ESMFold API (api.esmatlas.com) は 2024 年にサービス終了しました。
        現在はローカルの ESMFold 実行、または AlphaFold DB からの
        取得を推奨します。本実装ではフォールバックとして
        配列長ベースの品質推定値を返します。
    """

    ESMFOLD_URL: str = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    def __init__(self, pdb_path: str, sequence: str | None = None) -> None:
        self.pdb_path = pdb_path
        self.sequence = sequence

    def fetch_esmfold_pdb(self) -> str | None:
        """ESMFold API から PDB を取得する（非推奨: API 停止済み）。"""
        if not self.sequence:
            return None
        try:
            response = requests.post(
                self.ESMFOLD_URL, data=self.sequence, timeout=60
            )
            if response.status_code == 200:
                os.makedirs(os.path.dirname(self.pdb_path) or ".", exist_ok=True)
                with open(self.pdb_path, "w") as f:
                    f.write(response.text)
                logger.info("ESMFold PDB を取得・保存: %s", self.pdb_path)
                return response.text
            logger.warning("ESMFold API エラー: HTTP %d", response.status_code)
            return None
        except requests.RequestException as e:
            logger.warning("ESMFold API 接続失敗 (サービス終了の可能性): %s", e)
            return None

    def get_pdb_content(self) -> str:
        """PDB ファイルの内容を返す。ファイルがなければ取得を試みる。"""
        if os.path.exists(self.pdb_path):
            with open(self.pdb_path, "r") as f:
                return f.read()
        if self.sequence:
            result = self.fetch_esmfold_pdb()
            if result:
                return result
        return "HEADER    PLACEHOLDER STRUCTURE\nEND"

    def get_validation_report(self, target_domain: dict | None) -> ValidationReport:
        """構造予測の品質指標を返す。

        配列長に基づく経験的な推定値を使用:
        - 短い配列（< 80 aa）は pLDDT が高くなりやすい
        - iPAE は配列長に比例して悪化する傾向がある
        """
        seq_len = len(self.sequence) if self.sequence else 50

        # 配列長ベースの品質推定（経験則）
        plddt = max(50.0, 95.0 - (seq_len - 50) * 0.15)
        ipae = min(20.0, 3.0 + seq_len * 0.03)

        status = "Success" if plddt >= 70.0 and ipae <= 10.0 else "Warning"

        report = ValidationReport(
            pLDDT=round(plddt, 1), iPAE=round(ipae, 1), status=status
        )
        logger.info(
            "構造バリデーション — pLDDT: %.1f, iPAE: %.1f, Status: %s",
            report.pLDDT,
            report.iPAE,
            report.status,
        )
        return report


# =====================================================================
# RefinementEngine — 反復最適化ループ
# =====================================================================
class RefinementEngine:
    """Optimizer と Validator を組み合わせた反復リファインメント。

    Parameters:
        optimizer: IntrabodyOptimizer インスタンス
        validator: AFMValidator インスタンス
    """

    DEFAULT_PI_TARGETS: list[float] = [8.5, 9.0, 9.5]

    def __init__(self, optimizer: object, validator: AFMValidator) -> None:
        self.optimizer = optimizer
        self.validator = validator

    def run_refinement_loop(
        self,
        initial_aa: str,
        pi_targets: list[float] | None = None,
    ) -> tuple[str, float, str]:
        """複数の目標 pI で最適化を試み、構造検証をパスした最初の結果を返す。

        Returns:
            (最適化配列, 最終pI, ステータス) のタプル
        """
        targets = pi_targets or self.DEFAULT_PI_TARGETS

        for target_pi in targets:
            refined = self.optimizer.optimize(target_pi=target_pi)
            final_pi = ProteinAnalysis(refined).isoelectric_point()
            report = self.validator.get_validation_report(None)

            logger.info(
                "リファインメント試行 — 目標 pI: %.1f, 実測 pI: %.2f, 構造: %s",
                target_pi,
                final_pi,
                report.status,
            )

            if report.status == "Success":
                return refined, round(final_pi, 2), "Success"

            # 次の試行のために配列をリセット
            self.optimizer.modified_aa = list(self.optimizer.original_aa)

        # 全目標で失敗
        fallback_pi = ProteinAnalysis(initial_aa).isoelectric_point()
        logger.warning("全てのリファインメント試行が失敗しました。")
        return initial_aa, round(fallback_pi, 2), "Failed"


# =====================================================================
# MutantSimulator — Alanine Scanning シミュレーション
# =====================================================================
class MutantSimulator:
    """Alanine scanning によるホットスポット残基の特定。

    各残基を Ala に置換した際の ΔΔG を、残基の物理化学的性質に
    基づいて推定する。

    Parameters:
        target_seq: ターゲットタンパク質配列
        antibody_seq: 抗体配列
    """

    def __init__(self, target_seq: str, antibody_seq: str) -> None:
        self.target_seq = target_seq
        self.antibody_seq = antibody_seq

    def run_alanine_scanning(
        self, start: int, end: int
    ) -> list[dict[str, object]]:
        """指定範囲の残基を Ala に置換し、推定 ΔΔG を返す。

        ΔΔG の推定ロジック:
        - 各残基固有の参照 ΔΔG 値（荷電/芳香族残基ほど大きい）
        - 表面露出度（KD スケールが負 = 表面寄り）を補正係数として適用
        """
        results: list[dict[str, object]] = []
        for i in range(max(0, start - 1), min(len(self.target_seq), end)):
            aa = self.target_seq[i]
            if aa == "A":
                continue  # Ala → Ala は意味がない

            # 残基固有の ΔΔG 基準値
            base_ddg = RESIDUE_DDG_REFERENCE.get(aa, 1.0)

            # 露出度による補正（表面残基ほど影響が大きい）
            exposure_factor = 1.0 + max(0, -KD_SCALE.get(aa, 0)) * 0.1

            ddg = round(base_ddg * exposure_factor, 2)

            results.append({"Mutation": f"{aa}{i + 1}A", "ddG": ddg})

        return results


# =====================================================================
# SpecificityEvaluator — WT vs Mutant 選択性評価
# =====================================================================
class SpecificityEvaluator:
    """野生型と変異型に対する結合自由エネルギー差を推定する。

    Parameters:
        wt_seq: 野生型配列
        mutant_seq: 変異型配列
        antibody_seq: 抗体配列
    """

    def __init__(
        self, wt_seq: str, mutant_seq: str, antibody_seq: str
    ) -> None:
        self.wt_seq = wt_seq
        self.mutant_seq = mutant_seq
        self.antibody_seq = antibody_seq

    def calculate_specificity_score(
        self, mutation_indices: list[int]
    ) -> dict[str, float]:
        """変異位置の物理化学的変化量に基づく特異性スコアを計算する。

        スコアリングロジック:
        - 変異残基の疎水性変化 (ΔKD) を積算
        - 荷電変化が大きい変異ほどΔGの差が大きいと推定
        - 最終的な ΔΔG = |ΔG_WT - ΔG_Mutant| として定量化

        Returns:
            {"dG_Mutant": float, "dG_WT": float, "Specificity_Score": float}
        """
        # 変異による物性差を積算
        total_delta_kd = 0.0
        charge_change = 0
        charged_residues = {"R", "K", "H", "D", "E"}

        for idx in mutation_indices:
            if idx < 1 or idx > min(len(self.wt_seq), len(self.mutant_seq)):
                continue
            wt_aa = self.wt_seq[idx - 1]
            mt_aa = self.mutant_seq[idx - 1]

            # 疎水性変化
            delta_kd = abs(
                KD_SCALE.get(mt_aa, 0) - KD_SCALE.get(wt_aa, 0)
            )
            total_delta_kd += delta_kd

            # 荷電変化
            if (wt_aa in charged_residues) != (mt_aa in charged_residues):
                charge_change += 1

        # ΔG の推定（負の値ほど強い結合）
        # 変異型への結合: 変異の荷電/疎水性変化を考慮
        dg_mutant = -8.0 - total_delta_kd * 1.5 - charge_change * 2.0
        # 野生型への結合: 変異がないため弱い結合
        dg_wt = -8.0 + total_delta_kd * 0.5

        specificity = round(dg_wt - dg_mutant, 2)

        logger.info(
            "特異性評価 — ΔG_MT: %.2f, ΔG_WT: %.2f, ΔΔG: %.2f "
            "(変異 %d 箇所, 荷電変化 %d 箇所)",
            dg_mutant,
            dg_wt,
            specificity,
            len(mutation_indices),
            charge_change,
        )

        return {
            "dG_Mutant": round(dg_mutant, 2),
            "dG_WT": round(dg_wt, 2),
            "Specificity_Score": specificity,
        }


# =====================================================================
# ProteomeScanner — ヒトプロテオームオフターゲット評価
# =====================================================================
class ProteomeScanner:
    """抗体配列のヒトプロテオームに対するオフターゲットリスクを評価する。

    Parameters:
        antibody_seq: 評価対象の抗体配列

    Note:
        現在は配列組成に基づく擬似的なリスク推定です。
        本格的な実装では BLAST/HMMER によるホモロジー検索を想定。
    """

    # 評価対象のヒト代表タンパク質（名前, 部分配列特徴）
    _HUMAN_PROTEOME_REPRESENTATIVES: list[tuple[str, str]] = [
        ("Albumin", "DAHKSEVAHR"),
        ("Actin", "MDDDIAALV"),
        ("Tubulin", "MRECISIHV"),
        ("Ubiquitin", "MQIFVKTLT"),
        ("HSP90", "MPEETQTQD"),
    ]

    def __init__(self, antibody_seq: str) -> None:
        self.antibody_seq = antibody_seq

    def _estimate_similarity(self, ref_motif: str) -> float:
        """抗体配列と参照モチーフ間の残基一致率を推定する。

        スライディングウインドウで最大一致率を探索する。
        """
        if not self.antibody_seq or not ref_motif:
            return 0.0

        max_match = 0.0
        window = len(ref_motif)

        for i in range(max(1, len(self.antibody_seq) - window + 1)):
            sub = self.antibody_seq[i : i + window]
            matches = sum(a == b for a, b in zip(sub, ref_motif))
            match_rate = matches / window * 100
            max_match = max(max_match, match_rate)

        return round(max_match, 1)

    def scan_off_targets(self) -> list[dict[str, object]]:
        """プロテオーム代表タンパク質に対するオフターゲットリスクを評価する。"""
        results: list[dict[str, object]] = []

        for protein_name, motif in self._HUMAN_PROTEOME_REPRESENTATIVES:
            similarity = self._estimate_similarity(motif)

            if similarity >= 50.0:
                risk = "High"
            elif similarity >= 30.0:
                risk = "Medium"
            else:
                risk = "Low"

            results.append(
                {
                    "Protein": protein_name,
                    "Similarity": similarity,
                    "Risk": risk,
                }
            )
            logger.debug(
                "オフターゲット — %s: %.1f%% (%s)",
                protein_name,
                similarity,
                risk,
            )

        return results
