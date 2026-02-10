import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.express as px
from stmol import showmol
import py3Dmol

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData

from src.optimizer import IntrabodyOptimizer
from src.simulator import AntibodyDiscoveryEngine, InteractionSimulator, AFMValidator, RefinementEngine, MutantSimulator, SpecificityEvaluator, ProteomeScanner
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

st.set_page_config(page_title="Intrabody Studio Pro", layout="wide")

st.title("🧬 Intrabody Studio Pro")
st.markdown("De Novo設計から多重変異特異性評価、オフターゲットスキャンまでを網羅")
st.markdown("---")

# 1. デュアルターゲット入力
st.header("1. Target Strategy (Wild-type vs Mutant)")
col_wt, col_mt = st.columns(2)
with col_wt:
    wt_seq = st.text_area("Wild-type (WT) Sequence", height=150, value="MTEYKLVVVGAGGVGKSALTI...")
with col_mt:
    mt_seq = st.text_area("Mutant (MT) Sequence", height=150, value="MTEYKLVVVGAGDVGKSALTI...")

# 2. 変異検知 & ドメイン推薦
diff_indices = [i+1 for i, (a, b) in enumerate(zip(wt_seq, mt_seq)) if a != b]
discovery = AntibodyDiscoveryEngine(mt_seq)

st.header("2. Domain Selection & Exposure Score")
if st.button("🌟 Recommend Best Binding Sites"):
    recs = discovery.recommend_domains()
    st.table(recs)

col1, col2, col3 = st.columns(3)
start_pos = col1.number_input("Start Pos", value=min(diff_indices)-2 if diff_indices else 1)
end_pos = col2.number_input("End Pos", value=max(diff_indices)+2 if diff_indices else 20)
exposure = discovery.calculate_exposure_score(start_pos, end_pos)
col3.metric("Exposure Score", f"{exposure:.2f}")

# 3. 設計パイプライン実行
if st.button("🚀 Run Full Design Pipeline"):
    with st.spinner("Executing advanced simulation..."):
        # Discovery & Optimization
        raw_aa = discovery.discover_binder()
        optimizer = IntrabodyOptimizer(raw_aa)
        validator = AFMValidator("results/temp.pdb")
        refined_aa, final_pi, status = RefinementEngine(optimizer, validator).run_refinement_loop(raw_aa)
        
        # Specficity & Safety
        spec_eval = SpecificityEvaluator(wt_seq, mt_seq, refined_aa)
        spec_res = spec_eval.calculate_specificity_score(diff_indices)
        scanner = ProteomeScanner(refined_aa)
        off_targets = scanner.scan_off_targets()

        # 結果表示
        st.header("3. Results & Validation")
        c1, c2, c3 = st.columns(3)
        c1.metric("Specificity ($\Delta \Delta G$)", spec_res['Specificity_Score'])
        c2.metric("Final $pI$", final_pi)
        c3.metric("Validation Status", status)

        # 4. 可視化
        tab1, tab2, tab3 = st.tabs(["📊 Alanine Scan", "📈 Hydropathy", "🧊 3D Model"])
        with tab1:
            mut_sim = MutantSimulator(mt_seq, refined_aa)
            scan_df = pd.DataFrame(mut_sim.run_alanine_scanning(start_pos, end_pos))
            st.plotly_chart(px.bar(scan_df, x='Mutation', y='ddG', color='ddG', color_continuous_scale='Reds'))
        with tab2:
            st.subheader("Hydropathy Profile (Kyte-Doolittle)")
            analysis = ProteinAnalysis(refined_aa)
            
            # 引数をキーワード指定(window=...)ではなく、順番通りに渡します
            # 第1引数: window (9)
            # 第2引数: edge (0.4)
            # 第3引数: scale (ProtParamData.kd)
            chart_data = analysis.protein_scale(9, 0.4, ProtParamData.kd)
            
            # グラフ表示
            st.line_chart(chart_data)
            st.info("💡 スコアが 0 より大きいと疎水性（ベタつき）、小さいと親水性（溶けやすさ）を示します。")
        with tab3:
            view = py3Dmol.view(query='pdb:1B27', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}}); view.addSurface(py3Dmol.VDW, {'opacity': 0.3})
            showmol(view, height=500)

        # 5. レポート出力
        st.header("4. Export & Safety Report")
        st.table(pd.DataFrame(off_targets))
        
        # PDF生成
        report_gen = AnalysisReport({'name': "Project Mutant", 'start': start_pos, 'end': end_pos}, 
                                   {'sequence': refined_aa, 'pI': final_pi, 'energy': -11.5, 
                                    'specificity': spec_res['Specificity_Score'], 'off_targets': off_targets})
        report_path = "results/report.pdf"
        report_gen.generate(report_path)
        
        with open(report_path, "rb") as f:
            st.download_button("📄 Download PDF Analysis Report", f, file_name="analysis_report.pdf")