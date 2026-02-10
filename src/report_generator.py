from fpdf import FPDF
from datetime import datetime
from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
import os

class AnalysisReport:
    def __init__(self, target_info, results):
        self.target_info = target_info
        self.results = results
        self.pdf = FPDF()
        self.pdf.set_auto_page_break(auto=True, margin=15)
        self.pdf.add_page()

    def generate(self, output_path, gb_path=None):
        self.pdf.set_font("Helvetica", 'B', 20)
        self.pdf.cell(0, 15, "Intrabody Studio Analysis Report", ln=True, align='C')
        self.pdf.set_font("Helvetica", size=10)
        self.pdf.cell(0, 10, f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", ln=True, align='R')
        self.pdf.ln(5)

        self._add_section_title("1. Target Configuration")
        self.pdf.set_font("Helvetica", size=11)
        self.pdf.cell(0, 7, f"Target Name: {self.target_info['name']}", ln=True)
        self.pdf.cell(0, 7, f"Mutation Domain: {self.target_info['start']}-{self.target_info['end']}", ln=True)
        self.pdf.ln(5)

        self._add_section_title("2. Sequence Information")
        self.pdf.set_font("Helvetica", 'B', 10)
        self.pdf.cell(0, 7, "Optimized Amino Acid Sequence:", ln=True)
        self.pdf.set_font("Courier", size=9)
        self.pdf.multi_cell(0, 5, self.results['sequence'])
        self.pdf.ln(3)
        self.pdf.set_font("Helvetica", 'B', 10)
        self.pdf.cell(0, 7, "Constructed DNA Sequence:", ln=True)
        self.pdf.set_font("Courier", size=9)
        self.pdf.multi_cell(0, 5, self.results['dna_sequence'])
        self.pdf.ln(5)

        if gb_path and os.path.exists(gb_path):
            self._add_section_title("3. Vector Map Visualization")
            graphic_record = BiopythonTranslator().translate_record(gb_path)
            ax, _ = graphic_record.plot(figure_width=8)
            map_img_path = "results/vector_map_temp.png"
            ax.figure.savefig(map_img_path, bbox_inches='tight')
            plt.close(ax.figure)
            self.pdf.image(map_img_path, x=15, w=180)
            self.pdf.ln(5)

        self._add_section_title("4. Biophysical & Safety Analysis")
        self.pdf.set_font("Helvetica", size=11)
        self.pdf.cell(0, 7, f"Final Isoelectric Point (pI): {self.results['pI']}", ln=True)
        self.pdf.cell(0, 7, f"Specificity Score (ddG): {self.results['specificity']}", ln=True)
        self.pdf.ln(3)
        self.pdf.set_font("Helvetica", 'B', 10)
        self.pdf.cell(0, 7, "Off-target Assessment:", ln=True)
        for risk in self.results['off_targets']:
            self.pdf.set_font("Helvetica", size=10)
            self.pdf.cell(0, 6, f"- {risk['Protein']}: {risk['Risk']} (Similarity: {risk['Similarity']}%)", ln=True)
        self.pdf.ln(10)

        # --- 5. 判定基準の説明を追記 ---
        self._add_section_title("5. Interpretation Criteria")
        self.pdf.set_font("Helvetica", size=10)
        criteria_text = (
            "Success Criteria for Intracellular Functional Scaffolds:\n"
            "- Isoelectric Point (pI): Target range 8.5 - 9.5 to ensure cytoplasmic solubility.\n"
            "- Specificity Score (ddG): Values > 8.0 indicate high selectivity for Mutant vs WT.\n"
            "- Confidence (pLDDT): Scores > 70 confirm high-quality structural prediction.\n"
            "- Docking Accuracy (iPAE): Distances < 10.0 A indicate precise interface alignment."
        )
        self.pdf.multi_cell(0, 6, criteria_text)

        self.pdf.output(output_path)

    def _add_section_title(self, title):
        self.pdf.set_font("Helvetica", 'B', 14)
        self.pdf.set_fill_color(240, 240, 240)
        self.pdf.cell(0, 10, title, ln=True, fill=True)
        self.pdf.ln(2)