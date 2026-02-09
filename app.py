import streamlit as st
import pandas as pd
import os
from src.optimizer import IntrabodyOptimizer
from src.evaluator import BatchEvaluator
from src.simulator import InteractionSimulator
from src.vector_builder import VectorDesigner

st.set_page_config(page_title="Intrabody Designer", page_icon="🧬")

st.title("🧬 Intrabody Development Platform")
st.markdown("抗体配列の最適化、物性評価、およびベクター構築をブラウザ上で行います。")

# サイドバー：ターゲット設定
st.sidebar.header("Target Settings")
target_seq = st.sidebar.text_area("Target Protein Sequence", 
    value="MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH")

# メイン画面：抗体候補の入力
st.header("1. Input Antibody Candidates")
input_data = st.text_area("Candidate ID and Sequence (CSV format: ID, Sequence)", 
    value="ITB-Alpha, MAEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS\nITB-Beta, MAEVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSLLLL")

if st.button("Run Pipeline"):
    # データ整形
    candidate_dict = {}
    for line in input_data.strip().split("\n"):
        cid, seq = line.split(",")
        candidate_dict[cid.strip()] = seq.strip()

    # --- Step 1: Evaluation ---
    st.header("2. Batch Evaluation & Ranking")
    evaluator = BatchEvaluator(candidate_dict)
    ranking = evaluator.evaluate_all()
    st.dataframe(ranking)

    # --- Step 2: Optimization ---
    top_candidate = ranking.iloc[0]
    st.header(f"3. Optimization Result ({top_candidate['ID']})")
    
    optimizer = IntrabodyOptimizer(top_candidate['Sequence'])
    optimized_aa = optimizer.optimize(target_pi=9.2)
    
    col1, col2 = st.columns(2)
    col1.metric("Original pI", top_candidate['pI'])
    col2.metric("Optimized AA Length", len(optimized_aa))
    st.code(optimized_aa, language="text")

    # --- Step 3: Simulation ---
    st.header("4. Binding Simulation")
    simulator = InteractionSimulator(target_seq, optimized_aa)
    energy = simulator.predict_binding_energy()
    st.metric(r"Predicted Binding Energy ($\Delta G$)", f"{energy} kcal/mol")

    # --- Step 4: GenBank Generation ---
    st.header("5. Download Vector File")
    builder = VectorDesigner()
    os.makedirs("results", exist_ok=True)
    out_path = f"results/{top_candidate['ID']}_vector.gb"
    builder.build_genbank(top_candidate['ID'], optimized_aa, out_path)
    
    with open(out_path, "rb") as f:
        st.download_button(
            label="Download GenBank File",
            data=f,
            file_name=f"{top_candidate['ID']}_vector.gb",
            mime="text/plain"
        )