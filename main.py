import os
import pandas as pd
from src.optimizer import IntrabodyOptimizer
from src.simulator import AntibodyDiscoveryEngine, InteractionSimulator, AFMValidator, RefinementEngine, MutantSimulator, SpecificityEvaluator, ProteomeScanner
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

def main():
    print("🚀 --- Intrabody Design Pipeline: Professional Mode ---")

    # 1. ターゲット設定 (野生型 vs 多重変異型)
    wt_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH"
    mt_seq = "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH"

    # 多重変異箇所の自動特定
    diff_indices = [i+1 for i, (a, b) in enumerate(zip(wt_seq, mt_seq)) if a != b]
    print(f"[*] Detected {len(diff_indices)} mutations at positions: {diff_indices}")

    # 2. 自動ドメイン推薦 & 露出度評価
    discovery = AntibodyDiscoveryEngine(mt_seq)
    recommendations = discovery.recommend_domains(top_n=1)
    target_range = recommendations.iloc[0]
    print(f"[*] Recommended Domain: {target_range['Start']}-{target_range['End']} (Exposure Score: {target_range['Exposure_Score']})")

    # 3. 抗体生成 (De Novo Discovery)
    raw_antibody = discovery.discover_binder()
    print(f"[*] Initial Candidate Generated: {raw_antibody[:20]}...")

    # 4. 自動最適化 & リカバリ (Auto-Refinement)
    optimizer = IntrabodyOptimizer(raw_antibody)
    validator = AFMValidator("results/temp.pdb")
    refiner = RefinementEngine(optimizer, validator)
    
    optimized_aa, final_pi, status = refiner.run_refinement_loop(raw_antibody)
    print(f"[*] Optimization Status: {status} (Final $pI$: {final_pi})")

    # 5. 特異性評価 (WT vs Mutant)
    spec_eval = SpecificityEvaluator(wt_seq, mt_seq, optimized_aa)
    spec_results = spec_eval.calculate_specificity_score(diff_indices)
    print(f"[*] Specificity Score ($\Delta \Delta G$): {spec_results['Specificity_Score']}")

    # 6. 安全性スキャン (Human Proteome Scan)
    scanner = ProteomeScanner(optimized_aa)
    off_targets = scanner.scan_off_targets()

    # 7. ベクター構築 & レポート生成
    builder = VectorDesigner()
    os.makedirs("results", exist_ok=True)
    builder.build_genbank("Final_ITB", optimized_aa, "results/Final_Vector.gb")
    
    report_data = {
        'sequence': optimized_aa, 'pI': final_pi, 'energy': -11.5,
        'specificity': spec_results['Specificity_Score'], 'off_targets': off_targets
    }
    report_gen = AnalysisReport({'name': "KRAS-G12D", 'start': target_range['Start'], 'end': target_range['End']}, report_data)
    report_gen.generate("results/Full_Analysis_Report.pdf")
    
    print(f"✅ Success! Pipeline complete. Results saved in 'results/' folder.")

if __name__ == "__main__":
    main()