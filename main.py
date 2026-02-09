from src.optimizer import IntrabodyOptimizer
from src.evaluator import BatchEvaluator
from src.simulator import InteractionSimulator
from src.vector_builder import VectorDesigner
import os

def main():
    # 1. ターゲットと抗体候補の設定
    target_protein_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH" # KRAS example
    
    candidates = {
        "ITB-Alpha": "MAEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS",
        "ITB-Beta": "MAEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSLLLL"
    }

    print("--- Step 1: Batch Evaluation ---")
    evaluator = BatchEvaluator(candidates)
    ranking = evaluator.evaluate_all()
    print(ranking)

    # 最上位候補を選択
    top_id = ranking.iloc[0]['ID']
    top_seq = ranking.iloc[0]['Sequence']

    print(f"\n--- Step 2: Optimizing {top_id} ---")
    optimizer = IntrabodyOptimizer(top_seq)
    optimized_aa = optimizer.optimize(target_pi=9.2)
    print(f"Optimized AA: {optimized_aa[:30]}...")

    print("\n--- Step 3: Binding Simulation ---")
    simulator = InteractionSimulator(target_protein_seq, optimized_aa)
    energy = simulator.predict_binding_energy()
    print(f"Predicted Binding Energy (ΔG): {energy} kcal/mol")

    print("\n--- Step 4: Generating GenBank File ---")
    builder = VectorDesigner()
    out_file = f"results/{top_id}_vector.gb"
    builder.build_genbank(top_id, optimized_aa, out_file)
    print(f"Done! File saved to: {out_file}")

if __name__ == "__main__":
    main()