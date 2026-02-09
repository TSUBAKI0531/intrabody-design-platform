import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class BatchEvaluator:
    def __init__(self, candidate_dict):
        self.candidates = candidate_dict
        self.results = []

    def evaluate_all(self):
        for cid, seq in self.candidates.items():
            analysis = ProteinAnalysis(seq)
            pi = analysis.isoelectric_point()
            gravy = analysis.gravy() # 疎水性指標
            
            # 総合スコア: 溶解度が高い(GRAVYが低い)ほど、かつpIが7.2から遠いほど高評価
            score = (-gravy * 2) + abs(pi - 7.2)
            
            self.results.append({
                'ID': cid,
                'Sequence': seq,
                'pI': round(pi, 2),
                'GRAVY': round(gravy, 2),
                'Total_Score': round(score, 2)
            })
            
        return pd.DataFrame(self.results).sort_values(by='Total_Score', ascending=False)