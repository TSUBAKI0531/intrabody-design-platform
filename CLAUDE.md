# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install dependencies
pip install -r requirements.txt

# Run Streamlit Web UI
streamlit run app.py
# → http://localhost:8501

# Run CLI pipeline (default: KRAS-G12D demo)
python main.py

# CLI with custom FASTA inputs
python main.py --wt data/example_sequences/kras_wt.fasta \
               --mt data/example_sequences/kras_g12d.fasta \
               --output results/experiment_01
```

There are no tests or linters configured.

## Architecture

### Module structure

All modules import from `src/config.py` as the single source of truth for constants and thresholds. The flow is:

```
app.py / main.py
  ├─ simulator.py   → AntibodyDiscoveryEngine, RefinementEngine, SpecificityEvaluator, ProteomeScanner
  ├─ optimizer.py   → IntrabodyOptimizer
  ├─ evaluator.py   → BatchEvaluator
  ├─ vector_builder.py → VectorDesigner
  └─ report_generator.py → AnalysisReport
         ↑ all depend on ↑
              config.py
```

### Pipeline logic

`RefinementEngine.run_refinement_loop()` is the core bidirectional optimization loop:
1. **Forward** — `IntrabodyOptimizer.optimize()` substitutes FR residues to reach a target pI (CDR regions protected via `_identify_cdrs()`)
2. **Backward** — `AFMValidator.get_validation_report()` checks structural quality (pLDDT ≥ 70, iPAE ≤ 10)
3. If structure fails, the sequence is reset and the next pI target in `[8.5, 9.0, 9.5]` is tried

`IntrabodyOptimizer._identify_cdrs()` uses Cys-position-based CDR estimation (Kabat heuristic). Sequences with fewer than 2 Cys residues fall back to length-based estimation.

### Key design decisions

- **pI target range 8.5–9.5**: Intrabodies function in the reducing cytoplasm; higher pI promotes electrostatic repulsion to prevent aggregation
- **CDR protection**: Only FR (framework) residues are substitution targets; CDR residues at antigen-binding positions are locked
- **Streamlit 3D viewer**: Uses `streamlit.components.v1.html(view._make_html())` — **not** `stmol`. The stmol library was replaced to fix Streamlit Cloud compatibility
- **ESMFold API**: The `AFMValidator.fetch_esmfold_pdb()` method calls `api.esmatlas.com` which shut down in 2024. The validator falls back to sequence-length-based quality estimation

### Known issue

`main.py` imports `InteractionSimulator` from `src.simulator`, but this class is not defined there. Running `python main.py` raises `ImportError`. The class was removed during refactoring but the import was not cleaned up.

### Output files (written to `results/`)

| File | Format | Contents |
|---|---|---|
| `Final_Vector.gb` | GenBank | Lentiviral vector with Kozak + human-codon-optimized CDS |
| `Full_Analysis_Report.pdf` | PDF | Sequence, vector map, biophysical metrics, off-target assessment |
