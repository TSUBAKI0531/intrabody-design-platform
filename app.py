"""
Intrabody Studio Pro — Streamlit Web UI

細胞内抗体（イントラボディ）の設計パイプラインを
インタラクティブに実行するための Web アプリケーション。
"""

import os
import traceback

import streamlit as st
import pandas as pd
import plotly.express as px
import streamlit.components.v1 as components
import py3Dmol
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData

from src.optimizer import IntrabodyOptimizer
from src.simulator import (
    AntibodyDiscoveryEngine,
    AFMValidator,
    RefinementEngine,
    MutantSimulator,
    SpecificityEvaluator,
    ProteomeScanner,
)
from src.vector_builder import VectorDesigner
from src.report_generator import AnalysisReport

# ---------------------------------------------------------------------------
# 定数
# ---------------------------------------------------------------------------
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

EXAMPLE_WT = (
    "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)
EXAMPLE_MT = (
    "MTEYKLVVVGAGDVGKSALTIQLIQNHFVDEYDPTIEDSY"
    "RKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCV"
    "FAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAA"
    "RTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREI"
    "RQH"
)

# pI の許容範囲（Success 判定用）
PI_LOWER = 8.5
PI_UPPER = 9.5
SPECIFICITY_THRESHOLD = 8.0


# ---------------------------------------------------------------------------
# ユーティリティ
# ---------------------------------------------------------------------------
def validate_sequence(seq: str, label: str = "配列") -> str:
    """アミノ酸配列のバリデーション。

    Returns:
        正規化（空白除去・大文字化）された配列文字列。
    Raises:
        ValueError: 不正な文字が含まれる場合。
    """
    cleaned = seq.strip().upper().replace(" ", "").replace("\n", "")
    if not cleaned:
        raise ValueError(f"{label}が空です。アミノ酸配列を入力してください。")
    invalid_chars = set(cleaned) - VALID_AA
    if invalid_chars:
        raise ValueError(
            f"{label}に無効な文字が含まれています: {', '.join(sorted(invalid_chars))}  \n"
            f"使用可能: {', '.join(sorted(VALID_AA))}"
        )
    return cleaned


def detect_mutations(wt: str, mt: str) -> list[int]:
    """WT と MT 配列を比較し、変異位置（1-indexed）のリストを返す。"""
    return [i + 1 for i, (a, b) in enumerate(zip(wt, mt)) if a != b]


# ---------------------------------------------------------------------------
# パイプライン実行
# ---------------------------------------------------------------------------
def run_pipeline(
    wt_seq: str,
    mt_seq: str,
    start_pos: int,
    end_pos: int,
    diff_indices: list[int],
) -> dict:
    """設計パイプラインを実行し、結果辞書を返す。"""
    progress = st.progress(0, text="抗体候補を生成中...")

    discovery = AntibodyDiscoveryEngine(mt_seq)
    raw_aa = discovery.discover_binder()
    progress.progress(20, text="pI 最適化中...")

    optimizer = IntrabodyOptimizer(raw_aa)
    validator = AFMValidator("results/temp.pdb", sequence=raw_aa)
    refined_aa, final_pi, status = RefinementEngine(
        optimizer, validator
    ).run_refinement_loop(raw_aa)
    progress.progress(50, text="特異性評価中...")

    spec_eval = SpecificityEvaluator(wt_seq, mt_seq, refined_aa)
    spec_res = spec_eval.calculate_specificity_score(diff_indices)
    progress.progress(75, text="プロテオームスキャン中...")

    scanner = ProteomeScanner(refined_aa)
    off_targets = scanner.scan_off_targets()
    progress.progress(100, text="完了!")

    return {
        "refined_aa": refined_aa,
        "final_pi": final_pi,
        "status": status,
        "spec_res": spec_res,
        "off_targets": off_targets,
        "start_pos": start_pos,
        "end_pos": end_pos,
    }


# ---------------------------------------------------------------------------
# 結果表示
# ---------------------------------------------------------------------------
def display_results(res: dict, mt_seq: str) -> None:
    """パイプライン結果をタブ・メトリクスで表示する。"""
    st.header("3. 設計結果 & バリデーション")

    c1, c2, c3 = st.columns(3)
    c1.metric("Specificity Score", res["spec_res"]["Specificity_Score"])
    c2.metric("最終 pI", f"{res['final_pi']:.2f}")
    c3.metric("ステータス", res["status"])

    # --- 判定基準 ---
    with st.expander("📊 判定基準とスコアの見方"):
        st.markdown(
            """
| 指標 | Success 基準 | 役割 |
|:---|:---|:---|
| 等電点 (pI) | 8.5 – 9.5 | 細胞内での溶解性維持 |
| 特異性スコア (ΔΔG) | 8.0 以上 | 変異型への選択性 |
| 構造信頼度 (pLDDT) | 70 以上 | 構造予測の信頼性 |
| 界面精度 (iPAE) | 10.0 Å 以下 | 配置予測の正確性 |
"""
        )

    # --- 解析タブ ---
    tab1, tab2, tab3 = st.tabs(
        ["📊 Alanine Scan", "📈 Hydropathy Profile", "🧊 3D Structure"]
    )

    with tab1:
        mut_sim = MutantSimulator(mt_seq, res["refined_aa"])
        scan_df = pd.DataFrame(
            mut_sim.run_alanine_scanning(res["start_pos"], res["end_pos"])
        )
        if not scan_df.empty:
            fig = px.bar(
                scan_df,
                x="Mutation",
                y="ddG",
                color="ddG",
                color_continuous_scale="Reds",
                title="Alanine Scanning — ΔΔG per Residue",
            )
            fig.update_layout(xaxis_tickangle=-45)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Alanine Scan の結果がありません。")

    with tab2:
        try:
            analysis = ProteinAnalysis(res["refined_aa"])
            kd_values = analysis.protein_scale(ProtParamData.kd, 9, 0.4)
            st.line_chart(kd_values)
        except Exception as e:
            st.error(f"Hydropathy 解析中にエラーが発生しました: {e}")

    with tab3:
        view = py3Dmol.view(query="pdb:1B27", width=800, height=500)
        view.setStyle({"cartoon": {"color": "spectrum"}})
        components.html(view._make_html(), height=520, scrolling=False)
        st.caption("※ デモ用構造 (PDB: 1B27) を表示しています。")


def export_section(res: dict) -> None:
    """レポート・PDB ダウンロードセクション。"""
    st.header("4. エクスポート")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("📄 技術レポート (PDF) を生成"):
            try:
                builder = VectorDesigner()
                dna_seq = "".join(
                    builder.codon_table.get(aa, "NNN") for aa in res["refined_aa"]
                )
                os.makedirs("results", exist_ok=True)
                gb_path = "results/ITB_Vector_Final.gb"
                builder.build_genbank("ITB_Custom", res["refined_aa"], gb_path)

                report_data = {
                    "sequence": res["refined_aa"],
                    "dna_sequence": dna_seq,
                    "pI": f"{res['final_pi']:.2f}",
                    "specificity": res["spec_res"]["Specificity_Score"],
                    "off_targets": res["off_targets"],
                }
                report_gen = AnalysisReport(
                    {
                        "name": "Intrabody Project",
                        "start": res["start_pos"],
                        "end": res["end_pos"],
                    },
                    report_data,
                )
                pdf_path = "results/Technical_Report.pdf"
                report_gen.generate(pdf_path, gb_path=gb_path)

                with open(pdf_path, "rb") as f:
                    st.download_button(
                        "⬇ PDF をダウンロード",
                        f,
                        file_name="Intrabody_Report.pdf",
                        mime="application/pdf",
                    )
            except Exception as e:
                st.error(f"レポート生成に失敗しました: {e}")

    with col2:
        if st.button("🧊 3D モデル (.pdb) をダウンロード"):
            try:
                validator = AFMValidator(
                    "results/designed.pdb", sequence=res["refined_aa"]
                )
                with st.spinner("構造データを取得中..."):
                    pdb_data = validator.get_pdb_content()
                st.download_button(
                    label="⬇ PDB をダウンロード",
                    data=pdb_data,
                    file_name="designed_intrabody.pdb",
                    mime="text/plain",
                )
            except Exception as e:
                st.error(f"PDB 取得に失敗しました: {e}")


# ---------------------------------------------------------------------------
# メインアプリ
# ---------------------------------------------------------------------------
def main() -> None:
    st.set_page_config(
        page_title="Intrabody Studio Pro",
        page_icon="🧬",
        layout="wide",
    )

    # --- セッション初期化 ---
    if "analysis_results" not in st.session_state:
        st.session_state.analysis_results = None

    # --- サイドバー ---
    with st.sidebar:
        st.title("🧬 Intrabody Studio Pro")
        st.markdown("---")
        st.markdown(
            "細胞内抗体のデザインから\nベクター構築までを自動化する\n統合パイプラインです。"
        )
        st.markdown("---")
        st.markdown(
            "**技術スタック**  \n"
            "BioPython · AlphaFold2 · Plotly · Streamlit"
        )

    # --- メインコンテンツ ---
    st.title("🧬 Intrabody Studio Pro")
    st.markdown("---")

    # === Step 1: 配列入力 ===
    st.header("1. ターゲット配列の入力")

    col_wt, col_mt = st.columns(2)
    with col_wt:
        wt_raw = st.text_area(
            "Wild-type (WT) 配列",
            value=EXAMPLE_WT,
            height=150,
            help="標準的な20アミノ酸 (A-Y) で入力してください。",
        )
    with col_mt:
        mt_raw = st.text_area(
            "Mutant (MT) 配列",
            value=EXAMPLE_MT,
            height=150,
            help="変異を含むターゲット配列を入力してください。",
        )

    # バリデーション
    try:
        wt_seq = validate_sequence(wt_raw, "WT 配列")
        mt_seq = validate_sequence(mt_raw, "MT 配列")
    except ValueError as e:
        st.error(str(e))
        st.stop()

    if len(wt_seq) != len(mt_seq):
        st.error(
            f"WT ({len(wt_seq)} aa) と MT ({len(mt_seq)} aa) の配列長が一致しません。"
        )
        st.stop()

    # 変異検出
    diff_indices = detect_mutations(wt_seq, mt_seq)
    if diff_indices:
        st.success(f"✅ {len(diff_indices)} 箇所の変異を検出: 位置 {diff_indices}")
    else:
        st.warning("⚠️ WT と MT の間に変異が見つかりません。配列を確認してください。")

    # === Step 2: ドメイン選択 ===
    st.header("2. ドメイン選択")

    discovery = AntibodyDiscoveryEngine(mt_seq)
    col1, col2, col3 = st.columns(3)

    default_start = min(diff_indices) - 2 if diff_indices else 1
    default_end = max(diff_indices) + 2 if diff_indices else 20

    start_pos = col1.number_input(
        "開始位置", value=max(1, default_start), min_value=1
    )
    end_pos = col2.number_input(
        "終了位置",
        value=min(len(mt_seq), default_end),
        min_value=int(start_pos) + 1,
    )
    exposure = discovery.calculate_exposure_score(start_pos, end_pos)
    col3.metric("Exposure Score", f"{exposure:.2f}")

    # === パイプライン実行 ===
    st.markdown("---")
    if st.button("🚀 パイプラインを実行", type="primary", use_container_width=True):
        try:
            results = run_pipeline(wt_seq, mt_seq, start_pos, end_pos, diff_indices)
            st.session_state.analysis_results = results
        except Exception as e:
            st.error(f"パイプラインの実行中にエラーが発生しました:\n\n```\n{e}\n```")
            with st.expander("詳細なエラー情報"):
                st.code(traceback.format_exc())

    # === 結果表示 ===
    if st.session_state.analysis_results:
        res = st.session_state.analysis_results
        st.markdown("---")
        display_results(res, mt_seq)
        st.markdown("---")
        export_section(res)


if __name__ == "__main__":
    main()
