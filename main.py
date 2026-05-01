"""
Intrabody Design Pipeline — CLI エントリポイント

Usage:
    python main.py
    python main.py --wt data/wt.txt --mt data/mt.txt --output results/
    python main.py --help
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import dataclass, field

import pandas as pd

from src.optimizer import IntrabodyOptimizer
from src.simulator import (
    AntibodyDiscoveryEngine,
    AFMValidator,
    MutantSimulator,
    ProteomeScanner,
    RefinementEngine,
    SpecificityEvaluator,
)
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

# ---------------------------------------------------------------------------
# ロガー設定
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 定数
# ---------------------------------------------------------------------------
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

DEFAULT_WT = (
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)

DEFAULT_MT = (
    "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)


# ---------------------------------------------------------------------------
# データクラス
# ---------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    """パイプラインの実行設定。"""

    wt_seq: str
    mt_seq: str
    output_dir: str = "results"
    project_name: str = "KRAS-G12D"
    top_domains: int = 1


@dataclass
class PipelineResult:
    """パイプラインの実行結果。"""

    optimized_aa: str = ""
    final_pi: float = 0.0
    status: str = ""
    specificity_score: float = 0.0
    off_targets: list = field(default_factory=list)
    target_start: int = 0
    target_end: int = 0


# ---------------------------------------------------------------------------
# ユーティリティ
# ---------------------------------------------------------------------------
def load_sequence(path: str) -> str:
    """テキスト / FASTA ファイルから配列を読み込む。"""
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # FASTA 形式: ヘッダ行 ('>') を除外
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines).upper()


def validate_sequence(seq: str, label: str = "配列") -> str:
    """アミノ酸配列を検証し、正規化して返す。"""
    cleaned = seq.strip().upper().replace(" ", "").replace("\n", "")
    if not cleaned:
        raise ValueError(f"{label}が空です。")
    invalid = set(cleaned) - VALID_AA
    if invalid:
        raise ValueError(
            f"{label}に無効な文字が含まれています: {', '.join(sorted(invalid))}"
        )
    return cleaned


# ---------------------------------------------------------------------------
# パイプライン
# ---------------------------------------------------------------------------
def run_pipeline(config: PipelineConfig) -> PipelineResult:
    """イントラボディ設計パイプラインのメインロジック。"""
    result = PipelineResult()

    # --- 1. 変異箇所の検出 ---
    diff_indices = [
        i + 1
        for i, (a, b) in enumerate(zip(config.wt_seq, config.mt_seq))
        if a != b
    ]
    logger.info(
        "変異 %d 箇所を検出: 位置 %s", len(diff_indices), diff_indices
    )

    # --- 2. ドメイン推薦 & 露出度評価 ---
    discovery = AntibodyDiscoveryEngine(config.mt_seq)
    recommendations = discovery.recommend_domains(top_n=config.top_domains)
    target = recommendations.iloc[0]
    result.target_start = int(target["Start"])
    result.target_end = int(target["End"])
    logger.info(
        "推薦ドメイン: %d–%d (Exposure Score: %s)",
        result.target_start,
        result.target_end,
        target["Exposure_Score"],
    )

    # --- 3. De Novo 抗体生成 ---
    raw_antibody = discovery.discover_binder()
    logger.info("初期候補生成完了: %s...", raw_antibody[:20])

    # --- 4. 自動最適化 & リファインメント ---
    optimizer = IntrabodyOptimizer(raw_antibody)
    validator = AFMValidator("results/temp.pdb")
    refiner = RefinementEngine(optimizer, validator)
    result.optimized_aa, result.final_pi, result.status = (
        refiner.run_refinement_loop(raw_antibody)
    )
    logger.info(
        "最適化完了 — Status: %s, pI: %.2f",
        result.status,
        result.final_pi,
    )

    # --- 5. 特異性評価 (WT vs Mutant) ---
    spec_eval = SpecificityEvaluator(
        config.wt_seq, config.mt_seq, result.optimized_aa
    )
    spec_res = spec_eval.calculate_specificity_score(diff_indices)
    result.specificity_score = spec_res["Specificity_Score"]
    logger.info("特異性スコア (ΔΔG): %s", result.specificity_score)

    # --- 6. ヒトプロテオーム安全性スキャン ---
    scanner = ProteomeScanner(result.optimized_aa)
    result.off_targets = scanner.scan_off_targets()
    logger.info("オフターゲット候補: %d 件", len(result.off_targets))

    # --- 7. ベクター構築 & レポート生成 ---
    os.makedirs(config.output_dir, exist_ok=True)

    builder = VectorDesigner()
    gb_path = os.path.join(config.output_dir, "Final_Vector.gb")
    builder.build_genbank("Final_ITB", result.optimized_aa, gb_path)
    logger.info("GenBank ファイル生成: %s", gb_path)

    report_data = {
        "sequence": result.optimized_aa,
        "pI": result.final_pi,
        "energy": -11.5,
        "specificity": result.specificity_score,
        "off_targets": result.off_targets,
    }
    report_gen = AnalysisReport(
        {
            "name": config.project_name,
            "start": result.target_start,
            "end": result.target_end,
        },
        report_data,
    )
    pdf_path = os.path.join(config.output_dir, "Full_Analysis_Report.pdf")
    report_gen.generate(pdf_path)
    logger.info("PDF レポート生成: %s", pdf_path)

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """コマンドライン引数をパースする。"""
    parser = argparse.ArgumentParser(
        description="Intrabody Design Pipeline — CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "例:\n"
            "  python main.py\n"
            "  python main.py --wt data/wt.fasta --mt data/mt.fasta\n"
            "  python main.py --output results/experiment_01\n"
        ),
    )
    parser.add_argument(
        "--wt",
        type=str,
        default=None,
        help="WT 配列ファイル (FASTA / テキスト)。省略時はデフォルト KRAS-WT を使用。",
    )
    parser.add_argument(
        "--mt",
        type=str,
        default=None,
        help="MT 配列ファイル (FASTA / テキスト)。省略時はデフォルト KRAS-G12D を使用。",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results",
        help="出力ディレクトリ (デフォルト: results/)",
    )
    parser.add_argument(
        "--project-name",
        type=str,
        default="KRAS-G12D",
        help="プロジェクト名 (レポートに記載)",
    )
    return parser.parse_args()


def main() -> int:
    """エントリポイント。成功時 0、失敗時 1 を返す。"""
    args = parse_args()

    logger.info("🚀 Intrabody Design Pipeline を開始します")

    # --- 配列の読み込み ---
    try:
        wt_seq = (
            validate_sequence(load_sequence(args.wt), "WT 配列")
            if args.wt
            else validate_sequence(DEFAULT_WT, "WT 配列")
        )
        mt_seq = (
            validate_sequence(load_sequence(args.mt), "MT 配列")
            if args.mt
            else validate_sequence(DEFAULT_MT, "MT 配列")
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error("配列の読み込みに失敗しました: %s", e)
        return 1

    if len(wt_seq) != len(mt_seq):
        logger.error(
            "WT (%d aa) と MT (%d aa) の配列長が一致しません。",
            len(wt_seq),
            len(mt_seq),
        )
        return 1

    # --- パイプライン実行 ---
    config = PipelineConfig(
        wt_seq=wt_seq,
        mt_seq=mt_seq,
        output_dir=args.output,
        project_name=args.project_name,
    )

    try:
        result = run_pipeline(config)
    except Exception as e:
        logger.error("パイプラインの実行中にエラーが発生しました: %s", e)
        logger.debug("Traceback:", exc_info=True)
        return 1

    logger.info(
        "✅ パイプライン完了! 結果は '%s/' に保存されました。", config.output_dir
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
