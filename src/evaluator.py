"""
evaluator.py — 候補配列のバッチ評価エンジン

複数のイントラボディ候補配列を物理化学的指標で一括評価し、
ランキングする。溶解度（GRAVY）と等電点（pI）を組み合わせた
総合スコアで、細胞内環境に適した候補を選別する。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

logger = logging.getLogger(__name__)


@dataclass
class ScoringWeights:
    """バッチ評価のスコアリング重み。

    Attributes:
        gravy_weight: GRAVY（疎水性）への重み。
            負の GRAVY は親水性 = 細胞質内での溶解性向上を示す。
        pi_deviation_weight: 目標 pI からの偏差への重み。
        reference_pi: 目標とする等電点。
            細胞質の pH（~7.2）付近を基準にする場合と、
            イントラボディの推奨範囲（8.5–9.5）を基準にする場合がある。
    """

    gravy_weight: float = 2.0
    pi_deviation_weight: float = 1.0
    reference_pi: float = 9.0  # イントラボディ推奨範囲の中央値


class BatchEvaluator:
    """候補配列群の一括評価とランキング。

    Parameters:
        candidate_dict: {候補ID: アミノ酸配列} の辞書
        weights: スコアリング重み設定（省略時はデフォルト）

    Example:
        >>> candidates = {"CAN_001": "EVQLVESGG...", "CAN_002": "DIQMTQSPS..."}
        >>> evaluator = BatchEvaluator(candidates)
        >>> df = evaluator.evaluate_all()
        >>> print(df.head())
    """

    def __init__(
        self,
        candidate_dict: dict[str, str],
        weights: ScoringWeights | None = None,
    ) -> None:
        self.candidates = candidate_dict
        self.weights = weights or ScoringWeights()
        self._results: list[dict] = []

    def evaluate_all(self) -> pd.DataFrame:
        """全候補を評価し、総合スコアの降順でソートした DataFrame を返す。"""
        self._results.clear()

        for cid, seq in self.candidates.items():
            try:
                analysis = ProteinAnalysis(seq)
                pi = analysis.isoelectric_point()
                gravy = analysis.gravy()
                mw = analysis.molecular_weight()

                # 総合スコア:
                #   親水性が高い（GRAVY が低い）ほど高評価
                #   pI が目標値に近いほど高評価
                solubility_score = -gravy * self.weights.gravy_weight
                pi_score = -abs(pi - self.weights.reference_pi) * self.weights.pi_deviation_weight
                total_score = solubility_score + pi_score

                self._results.append(
                    {
                        "ID": cid,
                        "Sequence": seq,
                        "Length": len(seq),
                        "MW_kDa": round(mw / 1000, 1),
                        "pI": round(pi, 2),
                        "GRAVY": round(gravy, 3),
                        "Solubility_Score": round(solubility_score, 2),
                        "pI_Score": round(pi_score, 2),
                        "Total_Score": round(total_score, 2),
                    }
                )
            except Exception as e:
                logger.warning("候補 %s の評価に失敗: %s", cid, e)
                self._results.append(
                    {
                        "ID": cid,
                        "Sequence": seq,
                        "Length": len(seq),
                        "MW_kDa": None,
                        "pI": None,
                        "GRAVY": None,
                        "Solubility_Score": None,
                        "pI_Score": None,
                        "Total_Score": float("-inf"),
                    }
                )

        df = pd.DataFrame(self._results).sort_values(
            "Total_Score", ascending=False
        )
        logger.info(
            "バッチ評価完了: %d 候補中 %d 件を正常評価",
            len(self.candidates),
            sum(1 for r in self._results if r["pI"] is not None),
        )
        return df
