import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import plotly.express as px
from stmol import showmol
import py3Dmol

# Biopython関連のインポート
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData

# 自作モジュールのインポート
from src.optimizer import IntrabodyOptimizer
from src.simulator import AntibodyDiscoveryEngine, InteractionSimulator, AFMValidator, RefinementEngine, MutantSimulator, SpecificityEvaluator, ProteomeScanner
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

st.set_page_config(page_title="Intrabody Studio Pro", layout="wide")

st.title("🧬 Intrabody Studio Pro")
st.markdown("De Novo設計から多重変異特異性評価、安全性スキャン、技術レポート作成まで")
st.markdown("---")

# サイドバー設定
st.sidebar.header("Pipeline Settings")
use_refinement = st.sidebar.checkbox("Enable Auto-Refinement", value=True)
target_pi = st.sidebar.slider("Target pI", 4.0, 11.0, 9.5)

# 1. ターゲット戦略入力
st.header("1. Dual-Target Strategy (Wild-type vs Mutant)")
col_wt, col_mt = st.columns(2)
with col_wt:
    wt_seq = st.text_area("Wild-type (WT) Sequence", height=150, value="MTEYKLVVVGAGGVGKSALTI...")
with col_mt:
    mt_seq = st.text_area("Mutant (MT) Sequence", height=150, value="MTEYKLVVVGAGDVGKSALTI...")

# 変異箇所の自動検知
diff_indices = [i+1 for i, (a, b) in enumerate(zip(wt_seq, mt_seq)) if a != b]
discovery = AntibodyDiscoveryEngine(mt_seq)

# 2. ドメイン選択
st.header("2. Domain Selection & Exposure Score")
if st.button("🌟 Recommend Best Binding Sites"):
    recs = discovery.recommend_domains()
    st.table(recs)

col1, col2, col3 = st.columns(3)
start_pos = col1.number_input("Start Pos", value=min(diff_indices)-2 if diff_indices else 1)
end_pos = col2.number_input("End Pos", value=max(diff_indices)+2 if diff_indices else 20)
exposure = discovery.calculate_exposure_score(start_pos, end_pos)
col3.metric("Exposure Score", f"{exposure:.2f}")

# 3. 設計実行
st.markdown("---")
if st.button("🚀 Run Full Design Pipeline"):
    with st.spinner("Executing advanced simulation and optimization..."):
        # Discovery & Optimization
        raw_aa = discovery.discover_binder()
        optimizer = IntrabodyOptimizer(raw_aa)
        validator = AFMValidator("results/temp.pdb")
        
        if use_refinement:
            refined_aa, final_pi, status = RefinementEngine(optimizer, validator).run_refinement_loop(raw_aa)
        else:
            refined_aa = optimizer.optimize(target_pi=target_pi)
            final_pi = ProteinAnalysis(refined_aa).isoelectric_point()
            status = "Manual"
        
        # 特異性 & 安全性解析
        spec_eval = SpecificityEvaluator(wt_seq, mt_seq, refined_aa)
        spec_res = spec_eval.calculate_specificity_score(diff_indices)
        scanner = ProteomeScanner(refined_aa)
        off_targets = scanner.scan_off_targets()

        # 結果ダッシュボード
        st.header("3. Design Results & Validation")
        c1, c2, c3 = st.columns(3)
        c1.metric("Specificity Score ($\Delta \Delta G$)", spec_res['Specificity_Score'])
        c2.metric("Final $pI$", f"{final_pi:.2f}")
        c3.metric("Validation Status", status)

        # 4. 視覚化セクション
        tab1, tab2, tab3 = st.tabs(["📊 Alanine Scan", "📈 Hydropathy Profile", "🧊 3D Structure"])
        
        with tab1:
            mut_sim = MutantSimulator(mt_seq, refined_aa)
            scan_df = pd.DataFrame(mut_sim.run_alanine_scanning(start_pos, end_pos))
            fig_scan = px.bar(scan_df, x='Mutation', y='ddG', color='ddG', color_continuous_scale='Reds',
                             labels={'ddG': '$\Delta \Delta G$ (kcal/mol)'})
            st.plotly_chart(fig_scan, use_container_width=True)

        with tab2:
            st.subheader("Hydropathy Profile (Kyte-Doolittle)")
            analysis = ProteinAnalysis(refined_aa)
            # 正しい引数順序: (scale, window, edge)
            chart_data = analysis.protein_scale(ProtParamData.kd, 9, 0.4)
            st.line_chart(chart_data)
            st.info("💡 正の値は疎水性（凝集リスク）、負の値は親水性を示します。")

        with tab3:
            view = py3Dmol.view(query='pdb:1B27', width=800, height=500)
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.addSurface(py3Dmol.VDW, {'opacity': 0.3, 'color': 'white'})
            showmol(view, height=500)

        # 5. 安全性解析 & レポート出力
        st.header("4. Safety Assessment & Technical Report")
        st.table(pd.DataFrame(off_targets))
        
        # --- レポート生成ロジックの統合 ---
        if st.button("📄 Generate Full Technical Report (PDF)"):
            builder = VectorDesigner()
            # アミノ酸からDNA配列を生成
            dna_seq = "".join([builder.codon_table.get(aa, "NNN") for aa in refined_aa])
            
            # ベクターファイルの構築
            os.makedirs("results", exist_ok=True)
            gb_path = "results/ITB_Vector_Final.gb"
            builder.build_genbank("ITB_Custom", refined_aa, gb_path)
            
            # レポート用データの集計
            report_data = {
                'sequence': refined_aa,
                'dna_sequence': dna_seq,
                'pI': f"{final_pi:.2f}",
                'specificity': spec_res['Specificity_Score'],
                'off_targets': off_targets
            }
            
            target_info = {
                'name': "Intrabody Project Mutant",
                'start': start_pos,
                'end': end_pos
            }
            
            # PDFの生成
            pdf_path = "results/Technical_Report.pdf"
            report_gen = AnalysisReport(target_info, report_data)
            report_gen.generate(pdf_path, gb_path=gb_path)
            
            with open(pdf_path, "rb") as f:
                st.download_button(
                    label="Download Complete Technical Report",
                    data=f,
                    file_name="Intrabody_Technical_Report.pdf",
                    mime="application/pdf"
                )
            st.success("全ての情報（配列、図、解析結果）を含むレポートが生成されました。")