"""
report_generator.py — PDF 技術レポートの自動生成

パイプラインの解析結果を、配列情報・物性評価・安全性評価・
ベクターマップを含む構造化された PDF レポートとして出力する。
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from datetime import datetime

import matplotlib
matplotlib.use("Agg")  # GUI バックエンドを使わない
import matplotlib.pyplot as plt

from fpdf import FPDF

logger = logging.getLogger(__name__)


@dataclass
class TargetInfo:
    """レポートに記載するターゲット情報。"""

    name: str
    start: int
    end: int


@dataclass
class ReportData:
    """レポートに記載する解析結果データ。

    dict ではなく dataclass にすることで、
    必須項目の欠落をコード上で検出しやすくする。
    """

    sequence: str
    dna_sequence: str = ""
    pI: str = "N/A"
    specificity: float = 0.0
    off_targets: list[dict] = field(default_factory=list)
    energy: float = 0.0


class AnalysisReport:
    """イントラボディ設計パイプラインの技術レポートを PDF で生成する。

    Parameters:
        target_info: ターゲット情報（dict または TargetInfo）
        results: 解析結果データ（dict または ReportData）

    Example:
        >>> target = {"name": "KRAS-G12D", "start": 10, "end": 15}
        >>> data = {"sequence": "EVQL...", "dna_sequence": "GAGG...",
        ...         "pI": "9.12", "specificity": 8.3, "off_targets": []}
        >>> report = AnalysisReport(target, data)
        >>> report.generate("output/report.pdf")
    """

    def __init__(
        self,
        target_info: dict | TargetInfo,
        results: dict | ReportData,
    ) -> None:
        # dict でもデータクラスでも受け取れるようにする
        if isinstance(target_info, dict):
            self.target = TargetInfo(**target_info)
        else:
            self.target = target_info

        if isinstance(results, dict):
            self.data = ReportData(**{
                k: v for k, v in results.items()
                if k in ReportData.__dataclass_fields__
            })
        else:
            self.data = results

    def generate(self, output_path: str, gb_path: str | None = None) -> None:
        """PDF レポートを生成してファイルに出力する。

        Parameters:
            output_path: PDF の出力先パス
            gb_path: GenBank ファイルのパス（ベクターマップ描画用、任意）
        """
        pdf = FPDF()
        pdf.set_auto_page_break(auto=True, margin=15)
        pdf.add_page()

        # --- タイトル ---
        pdf.set_font("Helvetica", "B", 20)
        pdf.cell(0, 15, "Intrabody Studio Analysis Report", ln=True, align="C")
        pdf.set_font("Helvetica", size=10)
        pdf.cell(
            0, 10,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            ln=True, align="R",
        )
        pdf.ln(5)

        # --- 1. ターゲット設定 ---
        self._section_title(pdf, "1. Target Configuration")
        pdf.set_font("Helvetica", size=11)
        pdf.cell(0, 7, f"Target Name: {self.target.name}", ln=True)
        pdf.cell(
            0, 7,
            f"Mutation Domain: {self.target.start} - {self.target.end}",
            ln=True,
        )
        pdf.ln(5)

        # --- 2. 配列情報 ---
        self._section_title(pdf, "2. Sequence Information")
        self._sequence_block(pdf, "Optimized Amino Acid Sequence:", self.data.sequence)
        if self.data.dna_sequence:
            self._sequence_block(pdf, "Constructed DNA Sequence:", self.data.dna_sequence)
        pdf.ln(5)

        # --- 3. ベクターマップ（任意） ---
        if gb_path:
            self._add_vector_map(pdf, gb_path)

        # --- 4. 物性・安全性評価 ---
        self._section_title(pdf, "4. Biophysical & Safety Analysis")
        pdf.set_font("Helvetica", size=11)
        pdf.cell(0, 7, f"Final Isoelectric Point (pI): {self.data.pI}", ln=True)
        pdf.cell(
            0, 7,
            f"Specificity Score (ddG): {self.data.specificity}",
            ln=True,
        )
        pdf.ln(3)

        # オフターゲット
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(0, 7, "Off-target Assessment:", ln=True)
        if self.data.off_targets:
            for risk in self.data.off_targets:
                pdf.set_font("Helvetica", size=10)
                protein = risk.get("Protein", "Unknown")
                similarity = risk.get("Similarity", "N/A")
                risk_level = risk.get("Risk", "Unknown")
                pdf.cell(
                    0, 6,
                    f"  - {protein}: {risk_level} (Similarity: {similarity}%)",
                    ln=True,
                )
        else:
            pdf.set_font("Helvetica", size=10)
            pdf.cell(0, 6, "  No significant off-targets detected.", ln=True)
        pdf.ln(10)

        # --- 5. 判定基準 ---
        self._section_title(pdf, "5. Interpretation Criteria")
        pdf.set_font("Helvetica", size=10)
        criteria = (
            "Success Criteria for Intracellular Functional Scaffolds:\n"
            "- Isoelectric Point (pI): Target range 8.5 - 9.5 "
            "to ensure cytoplasmic solubility.\n"
            "- Specificity Score (ddG): Values > 8.0 indicate "
            "high selectivity for Mutant vs WT.\n"
            "- Confidence (pLDDT): Scores > 70 confirm "
            "high-quality structural prediction.\n"
            "- Docking Accuracy (iPAE): Distances < 10.0 A "
            "indicate precise interface alignment."
        )
        pdf.multi_cell(0, 6, criteria)

        # --- 出力 ---
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        pdf.output(output_path)
        logger.info("PDF レポート生成: %s", output_path)

    # ------------------------------------------------------------------
    # ヘルパーメソッド
    # ------------------------------------------------------------------
    @staticmethod
    def _section_title(pdf: FPDF, title: str) -> None:
        """セクション見出しを描画する。"""
        pdf.set_font("Helvetica", "B", 14)
        pdf.set_fill_color(240, 240, 240)
        pdf.cell(0, 10, title, ln=True, fill=True)
        pdf.ln(2)

    @staticmethod
    def _sequence_block(pdf: FPDF, label: str, sequence: str) -> None:
        """配列情報のブロックを描画する。"""
        pdf.set_font("Helvetica", "B", 10)
        pdf.cell(0, 7, label, ln=True)
        pdf.set_font("Courier", size=9)
        # 80 文字ごとに改行して見やすくする
        formatted = "\n".join(
            sequence[i : i + 80] for i in range(0, len(sequence), 80)
        )
        pdf.multi_cell(0, 5, formatted)
        pdf.ln(3)

    def _add_vector_map(self, pdf: FPDF, gb_path: str) -> None:
        """GenBank ファイルからベクターマップ画像を生成して PDF に挿入する。"""
        if not os.path.exists(gb_path):
            logger.warning("GenBank ファイルが見つかりません: %s", gb_path)
            return

        try:
            from dna_features_viewer import BiopythonTranslator

            graphic_record = BiopythonTranslator().translate_record(gb_path)
            ax, _ = graphic_record.plot(figure_width=8)
            map_img_path = gb_path.replace(".gb", "_map.png")
            ax.figure.savefig(map_img_path, bbox_inches="tight", dpi=150)
            plt.close(ax.figure)

            self._section_title(pdf, "3. Vector Map Visualization")
            pdf.image(map_img_path, x=15, w=180)
            pdf.ln(5)
            logger.info("ベクターマップ画像を生成: %s", map_img_path)

        except ImportError:
            logger.warning("dna-features-viewer が未インストールのため、ベクターマップをスキップします。")
        except Exception as e:
            logger.warning("ベクターマップ生成中にエラー: %s", e)
