"""
設計パイプラインの定数・閾値を一元管理するモジュール。

マジックナンバーをコード中に散在させず、ここで管理することで
閾値の変更や実験条件の切り替えを容易にする。
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# アミノ酸関連
# ---------------------------------------------------------------------------
VALID_AMINO_ACIDS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle 疎水性スケール
KD_HYDROPHOBICITY: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# ---------------------------------------------------------------------------
# pI 最適化パラメータ
# ---------------------------------------------------------------------------
PI_TARGET_LOWER: float = 8.5
PI_TARGET_UPPER: float = 9.5
PI_TARGET_DEFAULT: float = 9.5
PI_REFINEMENT_TARGETS: list[float] = [8.5, 9.0, 9.5]
PI_FAILED_FALLBACK: float = 7.0

# pI を上げるための置換候補（優先順: 置換のインパクトが小さい順）
PI_INCREASE_DONORS: str = "LIVFMA"       # 疎水性 / 小さい残基
PI_INCREASE_ACCEPTORS: list[str] = ["K", "R", "H"]  # 塩基性残基

# pI を下げるための置換候補
PI_DECREASE_DONORS: str = "RK"
PI_DECREASE_ACCEPTORS: list[str] = ["Q", "N"]

# ---------------------------------------------------------------------------
# 特異性評価
# ---------------------------------------------------------------------------
SPECIFICITY_DDG_THRESHOLD: float = 8.0

# ---------------------------------------------------------------------------
# 構造予測バリデーション
# ---------------------------------------------------------------------------
PLDDT_THRESHOLD: float = 70.0
IPAE_THRESHOLD: float = 10.0

# ESMFold API（Meta が 2024 年にサービス縮小 — フォールバック考慮必須）
ESMFOLD_API_URL: str = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_TIMEOUT_SEC: int = 60

# ---------------------------------------------------------------------------
# バッチ評価
# ---------------------------------------------------------------------------
EVALUATOR_PI_REFERENCE: float = 7.2  # 中性 pH 基準
EVALUATOR_GRAVY_WEIGHT: float = 2.0

# ---------------------------------------------------------------------------
# ベクター構築
# ---------------------------------------------------------------------------
DEFAULT_VECTOR_NAME: str = "pLenti-ITB"

# ヒト最適化コドンテーブル（高頻度コドンを選択）
# 参考: Kazusa Codon Usage Database — Homo sapiens
HUMAN_CODON_TABLE: dict[str, str] = {
    "A": "GCC", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "TCC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}

# Kozak コンセンサス配列（翻訳開始効率向上）
KOZAK_SEQUENCE: str = "GCCACC"

# プロモーター配列辞書
PROMOTER_SEQUENCES: dict[str, str] = {
    "EF1a": (
        "GGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCC"
        "GAGAAGTTGGGGGGAGGGGTCG"
    ),
    "CMV": (
        "TAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATAT"
        "ATGGAGTTCCGCGTTACATAACTTACGG"
    ),
}

# ---------------------------------------------------------------------------
# プロテオームスキャン
# ---------------------------------------------------------------------------
OFF_TARGET_IDENTITY_CUTOFF: float = 0.3

# ---------------------------------------------------------------------------
# レポート生成
# ---------------------------------------------------------------------------
REPORT_TITLE: str = "Intrabody Studio Analysis Report"
REPORT_FONT_TITLE: int = 20
REPORT_FONT_SECTION: int = 14
REPORT_FONT_BODY: int = 11
REPORT_FONT_SMALL: int = 10
REPORT_FONT_CODE: int = 9

# ---------------------------------------------------------------------------
# デフォルトサンプル配列 (KRAS)
# ---------------------------------------------------------------------------
KRAS_WT: str = (
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)

KRAS_G12D: str = (
    "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)
