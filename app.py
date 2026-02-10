import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.express as px
from stmol import showmol
import py3Dmol

# Biopython
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData

# 自作モジュール
from src.optimizer import IntrabodyOptimizer
from src.simulator import AntibodyDiscoveryEngine, AFMValidator, RefinementEngine, MutantSimulator, SpecificityEvaluator, ProteomeScanner
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

st.set_page_config(page_title="Intrabody Studio Pro", layout="wide")

# セッション状態の初期化
if 'analysis_results' not in st.session_state:
    st.session_state.analysis_results = None

st.title("🧬 Intrabody Studio Pro")
st.markdown("---")

# 1. ターゲット入力
st.header("1. Target Strategy")
col_wt, col_mt = st.columns(2)
with col_wt:
    wt_seq = st.text_area("Wild-type (WT) Sequence", height=150, value="MTEYKLVVVGAGGVGKSALTI...")
with col_mt:
    mt_seq = st.text_area("Mutant (MT) Sequence", height=150, value="MTEYKLVVVGAGDVGKSALTI...")

# 変異検知
diff_indices = [i+1 for i, (a, b) in enumerate(zip(wt_seq, mt_seq)) if a != b]
discovery = AntibodyDiscoveryEngine(mt_seq)

# 2. ドメイン選択
st.header("2. Domain Selection")
col1, col2, col3 = st.columns(3)
start_pos = col1.number_input("Start Position", value=min(diff_indices)-2 if diff_indices else 1)
end_pos = col2.number_input("End Position", value=max(diff_indices)+2 if diff_indices else 20)
exposure = discovery.calculate_exposure_score(start_pos, end_pos)
col3.metric("Exposure Score", f"{exposure:.2f}")

# 設計実行ボタン
if st.button("🚀 Run Full Design Pipeline"):
    with st.spinner("Analyzing and Optimizing..."):
        raw_aa = discovery.discover_binder()
        optimizer = IntrabodyOptimizer(raw_aa)
        validator = AFMValidator("results/temp.pdb")
        refined_aa, final_pi, status = RefinementEngine(optimizer, validator).run_refinement_loop(raw_aa)
        
        spec_eval = SpecificityEvaluator(wt_seq, mt_seq, refined_aa)
        spec_res = spec_eval.calculate_specificity_score(diff_indices)
        scanner = ProteomeScanner(refined_aa)
        off_targets = scanner.scan_off_targets()

        # 結果をメモリに保存 
        st.session_state.analysis_results = {
            'refined_aa': refined_aa,
            'final_pi': final_pi,
            'status': status,
            'spec_res': spec_res,
            'off_targets': off_targets,
            'start_pos': start_pos,
            'end_pos': end_pos
        }

# 結果の表示
if st.session_state.analysis_results:
    res = st.session_state.analysis_results
    
    st.header("3. Design Results & Validation")
    c1, c2, c3 = st.columns(3)
    c1.metric("Specificity Score ($\Delta \Delta G$)", res['spec_res']['Specificity_Score'])
    c2.metric("Final $pI$", f"{res['final_pi']:.2f}")
    c3.metric("Status", res['status'])

    # 解説パネル 
    with st.expander("📊 判定基準とスコアの見方について"):
        st.markdown("""
        ### **Success判定のガイドライン**
        | 指標 | Success基準 | 役割 |
        | :--- | :--- | :--- |
        | **等電点 (pI)** | **8.5 - 9.5** | 細胞内での凝集を防ぎ、溶解性を維持するために必要です。 |
        | **特異性スコア (ΔΔG)** | **8.0以上** | 野生型を避け、変異型のみを狙い撃ちできているかを示します。 |
        | **構造信頼度 (pLDDT)** | **70以上** | 構造予測が信頼できるレベルかを示します。 |
        | **界面精度 (iPAE)** | **10.0Å以下** | ターゲットとの配置が正確に予測できているかを示します。 |
        """)

    tab1, tab2, tab3 = st.tabs(["📊 Alanine Scan", "📈 Hydropathy Profile", "🧊 3D Structure"])
    
    with tab1:
        mut_sim = MutantSimulator(mt_seq, res['refined_aa'])
        scan_df = pd.DataFrame(mut_sim.run_alanine_scanning(res['start_pos'], res['end_pos']))
        st.plotly_chart(px.bar(scan_df, x='Mutation', y='ddG', color='ddG', color_continuous_scale='Reds'))

    with tab2:
        analysis = ProteinAnalysis(res['refined_aa'])
        # 正しい引数順序: (scale, window, edge)
        chart_data = analysis.protein_scale(ProtParamData.kd, 9, 0.4)
        st.line_chart(chart_data)

    with tab3:
        view = py3Dmol.view(query='pdb:1B27', width=800, height=500)
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        showmol(view, height=500)

    # 4. エクスポートセクション 
    st.header("4. Export Report & 3D Model")
    btn_col1, btn_col2 = st.columns(2)
    
    with btn_col1:
        if st.button("📄 Generate Full Technical Report (PDF)"):
            builder = VectorDesigner()
            dna_seq = "".join([builder.codon_table.get(aa, "NNN") for aa in res['refined_aa']])
            os.makedirs("results", exist_ok=True)
            gb_path = "results/ITB_Vector_Final.gb"
            builder.build_genbank("ITB_Custom", res['refined_aa'], gb_path)
            
            report_data = {
                'sequence': res['refined_aa'], 'dna_sequence': dna_seq,
                'pI': f"{res['final_pi']:.2f}", 'specificity': res['spec_res']['Specificity_Score'],
                'off_targets': res['off_targets']
            }
            report_gen = AnalysisReport({'name': "Intrabody Project", 'start': res['start_pos'], 'end': res['end_pos']}, report_data)
            pdf_path = "results/Technical_Report.pdf"
            report_gen.generate(pdf_path, gb_path=gb_path)
            
            with open(pdf_path, "rb") as f:
                st.download_button("Download PDF", f, file_name="Intrabody_Report.pdf")

    with btn_col2:
        # PDBファイルのダウンロード 
        validator = AFMValidator("results/temp.pdb")
        pdb_data = validator.get_pdb_content()
        st.download_button(
            label="🧊 Download 3D Model (.pdb)",
            data=pdb_data,
            file_name="designed_intrabody.pdb",
            mime="text/plain",
            help="PyMOLやCueMol2で解析可能なPDB形式です。"
        )