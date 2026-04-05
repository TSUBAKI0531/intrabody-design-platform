"""
optimizer.py — イントラボディの細胞内安定性を確保するための pI 最適化エンジン

細胞内（還元環境）での凝集を防ぐため、フレームワーク領域（FR）の
アミノ酸残基を置換し、等電点（pI）を目標範囲に調整する。
CDR（相補性決定領域）は抗原結合に直結するため、置換対象から保護する。
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 定数
# ---------------------------------------------------------------------------
# pI を上げたい場合に FR 領域で置換する候補（塩基性残基）
BASIC_SUBSTITUTES: list[str] = ["K", "R"]
# pI を下げたい場合に FR 領域で置換する候補（酸性残基）
ACIDIC_SUBSTITUTES: list[str] = ["D", "E"]
# 置換対象にできる疎水性 FR 残基
HYDROPHOBIC_TARGETS: set[str] = set("LIVFMA")
# pI 許容誤差
PI_TOLERANCE: float = 0.3


@dataclass
class MutationRecord:
    """1 回の置換操作を記録するデータクラス。"""

    position: int
    original: str
    substituted: str
    pi_before: float
    pi_after: float


class IntrabodyOptimizer:
    """FR 領域のアミノ酸を置換して pI を目標値に近づける最適化エンジン。

    Parameters:
        sequence: 最適化対象のアミノ酸配列
    """

    def __init__(self, sequence: str) -> None:
        self.original_aa: str = sequence
        self.modified_aa: list[str] = list(sequence)
        self.cdr_indices: set[int] = self._identify_cdrs()
        self.mutation_log: list[MutationRecord] = []

    def _identify_cdrs(self) -> set[int]:
        """Cys 残基の位置に基づいて CDR 領域を推定し、保護対象インデックスを返す。

        Kabat 番号付けの経験則:
        - CDR1: 1番目の Cys から +23 ～ +35
        - CDR3: 最後の Cys から +2 ～ +15

        Cys が 2 個未満の場合は、配列の 25–35% と 85–95% を概算 CDR とする
        フォールバック処理を行う。
        """
        indices: set[int] = set()
        cys_positions = [m.start() for m in re.finditer("C", self.original_aa)]

        if len(cys_positions) >= 2:
            # 標準的な推定
            indices.update(range(cys_positions[0] + 23, cys_positions[0] + 35))
            indices.update(range(cys_positions[-1] + 2, cys_positions[-1] + 15))
        else:
            # フォールバック: 配列長ベースの概算
            seq_len = len(self.original_aa)
            indices.update(range(int(seq_len * 0.25), int(seq_len * 0.35)))
            indices.update(range(int(seq_len * 0.85), int(seq_len * 0.95)))
            logger.warning(
                "Cys が %d 個しか検出されませんでした。フォールバック CDR 推定を使用します。",
                len(cys_positions),
            )

        # 範囲外のインデックスを除去
        valid_range = set(range(len(self.original_aa)))
        return indices & valid_range

    def _current_pi(self) -> float:
        """現在の配列の等電点を計算する。"""
        return ProteinAnalysis("".join(self.modified_aa)).isoelectric_point()

    def optimize(self, target_pi: float = 9.0, max_mutations: int = 15) -> str:
        """等電点を目標値に近づけるための置換を実行する。

        Parameters:
            target_pi: 目標 pI（デフォルト 9.0）
            max_mutations: 最大置換回数（過剰改変の防止）

        Returns:
            最適化後のアミノ酸配列
        """
        self.mutation_log.clear()
        current_pi = self._current_pi()
        logger.info("最適化開始 — 初期 pI: %.2f, 目標 pI: %.2f", current_pi, target_pi)

        mutation_count = 0

        for i in range(len(self.modified_aa)):
            if mutation_count >= max_mutations:
                logger.info("最大置換回数 (%d) に到達。", max_mutations)
                break

            current_pi = self._current_pi()
            if abs(current_pi - target_pi) <= PI_TOLERANCE:
                logger.info("目標 pI 範囲内に到達: %.2f", current_pi)
                break

            # CDR 領域はスキップ
            if i in self.cdr_indices:
                continue

            # 置換対象の残基でなければスキップ
            if self.modified_aa[i] not in HYDROPHOBIC_TARGETS:
                continue

            # pI を上げるか下げるかで置換残基を選択
            pi_before = current_pi
            original_aa = self.modified_aa[i]

            if current_pi < target_pi:
                # 塩基性残基に置換して pI を上昇させる
                substitute = BASIC_SUBSTITUTES[mutation_count % len(BASIC_SUBSTITUTES)]
            else:
                # 酸性残基に置換して pI を低下させる
                substitute = ACIDIC_SUBSTITUTES[mutation_count % len(ACIDIC_SUBSTITUTES)]

            self.modified_aa[i] = substitute
            pi_after = self._current_pi()
            mutation_count += 1

            self.mutation_log.append(
                MutationRecord(
                    position=i + 1,
                    original=original_aa,
                    substituted=substitute,
                    pi_before=round(pi_before, 2),
                    pi_after=round(pi_after, 2),
                )
            )
            logger.debug(
                "置換 #%d: %s%d%s (pI: %.2f → %.2f)",
                mutation_count,
                original_aa,
                i + 1,
                substitute,
                pi_before,
                pi_after,
            )

        final_pi = self._current_pi()
        logger.info(
            "最適化完了 — 最終 pI: %.2f, 置換数: %d", final_pi, mutation_count
        )
        return "".join(self.modified_aa)

    def get_mutation_summary(self) -> list[dict]:
        """置換履歴をリスト[dict]で返す（レポート用）。"""
        return [
            {
                "Position": m.position,
                "Original": m.original,
                "Substituted": m.substituted,
                "pI_before": m.pi_before,
                "pI_after": m.pi_after,
            }
            for m in self.mutation_log
        ]
