import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class IntrabodyOptimizer:
    def __init__(self, sequence):
        self.original_aa = sequence
        self.modified_aa = list(sequence)
        self.kd_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        self.cdr_indices = self._identify_cdrs()

    def _identify_cdrs(self):
        """Kabatの経験則に基づきCDR領域を特定して保護する"""
        indices = set()
        cys_positions = [m.start() for m in re.finditer('C', self.original_aa)]
        if len(cys_positions) >= 2:
            indices.update(range(cys_positions[0]+23, cys_positions[0]+35)) # CDR1推定
            indices.update(range(cys_positions[-1]+2, cys_positions[-1]+15)) # CDR3推定
        return indices

    def optimize(self, target_pi=9.0):
        """等電点を目標値に近づけつつ、疎水性を上げないように置換を行う"""
        current_pi = ProteinAnalysis("".join(self.modified_aa)).isoelectric_point()
        
        # pIを上げるために、FR領域の疎水性残基(L, I, V, F)をR(Arg)に置換
        if current_pi < target_pi:
            for i in range(len(self.modified_aa)):
                if i not in self.cdr_indices and self.modified_aa[i] in 'LIVF':
                    self.modified_aa[i] = 'R'
                    if ProteinAnalysis("".join(self.modified_aa)).isoelectric_point() >= target_pi:
                        break
        return "".join(self.modified_aa)