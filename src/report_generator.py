from fpdf import FPDF
from datetime import datetime

class AnalysisReport:
    def __init__(self, target_info, results):
        self.target_info = target_info
        self.results = results
        self.pdf = FPDF()
        self.pdf.add_page()

    def generate(self, output_path):
        self.pdf.set_font("Helvetica", 'B', 16)
        self.pdf.cell(0, 10, "Intrabody Studio Analysis Report", ln=True, align='C')
        self.pdf.ln(10)
        
        self.pdf.set_font("Helvetica", size=11)
        self.pdf.cell(0, 7, f"Target: {self.target_info['name']} (Domain: {self.target_info['start']}-{self.target_info['end']})", ln=True)
        self.pdf.cell(0, 7, f"Final pI: {self.results['pI']}", ln=True)
        self.pdf.cell(0, 7, f"Specificity Score (ddG): {self.results['specificity']}", ln=True)
        self.pdf.ln(5)
        
        self.pdf.set_font("Helvetica", 'B', 12)
        self.pdf.cell(0, 7, "Off-target Risks:", ln=True)
        self.pdf.set_font("Helvetica", size=10)
        for ot in self.results['off_targets']:
            self.pdf.cell(0, 6, f"- {ot['Protein']}: {ot['Risk']} (Similarity: {ot['Similarity']}%)", ln=True)
        
        self.pdf.output(output_path)