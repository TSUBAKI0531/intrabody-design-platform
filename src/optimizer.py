import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class IntrabodyOptimizer:
    def __init__(self, sequence):
        self.original_aa = sequence
        self.modified_aa = list(sequence)
        self.cdr_indices = self._identify_cdrs()

    def _identify_cdrs(self):
        """Kabat番号の経験則に基づきCDRを保護"""
        indices = set()
        cys_positions = [m.start() for m in re.finditer('C', self.original_aa)]
        if len(cys_positions) >= 2:
            indices.update(range(cys_positions[0]+23, cys_positions[0]+35)) # CDR1
            indices.update(range(cys_positions[-1]+2, cys_positions[-1]+15)) # CDR3
        return indices

    def optimize(self, target_pi=9.5):
        """等電点を目標値に近づけるための置換ロジック"""
        analysis = ProteinAnalysis("".join(self.modified_aa))
        current_pi = analysis.isoelectric_point()
        
        # pIを上げるために、FR領域の疎水性残基をR(Arg)に置換
        if current_pi < target_pi:
            for i in range(len(self.modified_aa)):
                if i not in self.cdr_indices and self.modified_aa[i] in 'LIVF':
                    self.modified_aa[i] = 'R'
                    if ProteinAnalysis("".join(self.modified_aa)).isoelectric_point() >= target_pi:
                        break
        return "".join(self.modified_aa)