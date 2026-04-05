# src/ コアモジュール — コードレビュー＆改善レポート

## 総合所見

6ファイル全体を通して、**生物学的な設計意図は的確**です。pI最適化・ΔΔG特異性評価・プロテオームスキャン・ベクター構築という一連のパイプラインは、イントラボディ研究の実際のワークフローをよく反映しています。一方で、**ソフトウェアとしての品質面**に改善余地が大きいため、以下に整理します。

---

## 1. simulator.py（最重大）

### 問題点

**ファイルが巨大 & 責務過剰（God Object パターン）**
- 1ファイルに6クラスが詰め込まれている
- AntibodyDiscoveryEngine, AFMValidator, RefinementEngine, MutantSimulator, SpecificityEvaluator, ProteomeScanner はそれぞれ独立した責務
- ポートフォリオ評価者は「設計力」をまず見るため、ここが最も印象を左右する

**プレースホルダーの明示不足**
- `discover_binder()` がハードコード文字列を返すが、docstring やコメントで「ここは将来的にde novo生成エンジンに置き換える設計」と明記されていない
- `SpecificityEvaluator.calculate_specificity_score()` がハードコード値を返す → 実装中なのか意図的なモックなのか判別不能
- `ProteomeScanner.scan_off_targets()` も同様
- **ポートフォリオにおいて「ハードコードされた値」は"手抜き"ではなく"段階的開発の設計判断"として説明する必要がある**

**ESMFold APIのURL**
- `https://api.esmatlas.com/foldSequence/v1/pdb/` は Meta が2024年にサービス終了
- 現在動作しない可能性が高い → フォールバック処理が不十分

### 改善方針
- ファイルを責務別に分割（discovery.py / validation.py / refinement.py / specificity.py）
  → ただしポートフォリオ的にはファイル数が増えすぎるのも良くないため、
    **最低限 simulator.py 内にセクションコメントを入れ、将来の分割を示唆する docstring を追加**
- 全プレースホルダーに `# TODO:` コメントと docstring で設計意図を明記
- ESMFold URL を定数化し、フォールバックを改善

---

## 2. optimizer.py

### 問題点

**CDR同定ロジックが脆弱**
- `_identify_cdrs()` が Cys 残基位置のみに依存 → Cys がない/位置が異なる配列で破綻
- 実際の Kabat/IMGT 番号付けとは異なる独自ヒューリスティック

**最適化が一方向のみ**
- pI を上げる処理（LIVF → R）のみで、下げる処理がない
- `target_pi` を下回る場合の処理がない

**置換戦略が単純すぎる**
- 全ての疎水性残基を R (Arg) のみに置換 → 配列の多様性を損なう
- K (Lys), H (His) なども塩基性残基として選択肢に入るべき

**ガード条件なし**
- 入力配列が空/不正な場合のチェックなし

### 改善方針
- docstring で「簡易版CDR同定であり、プロダクション利用にはANARCI等を推奨」と明記
- pI を下げる方向の置換も追加（R/K → Q/N）
- 置換候補をリスト化して戦略パターンを明確化
- 入力バリデーション追加

---

## 3. evaluator.py

### 問題点
- **app.py / main.py のどちらからもインポートされていない** → 未使用モジュール
- スコアリング関数 `(-gravy * 2) + abs(pi - 7.2)` の根拠が不明
- pI 7.2 からの距離を評価しているが、イントラボディのpI目標は 8.5–9.5 であり矛盾

### 改善方針
- config.py の定数（PI_TARGET）を参照するよう修正
- docstring でスコアリング根拠を説明
- main.py / app.py から呼び出せるように統合ポイントを明示

---

## 4. vector_builder.py

### 問題点

**ベクター要素が不足**
- Kozak 配列（GCCACC）なし → 翻訳開始効率に影響
- 終止コドンの明示的追加なし
- PolyA シグナル、WPRE 要素なし
- プロモーター選択肢が EF1a のみ

**コドンテーブルが最小限**
- ヒト最適化コドンとしては妥当だが、コドン使用頻度の根拠が不明

### 改善方針
- Kozak + 終止コドンの自動付加
- プロモーター選択オプション（CMV / EF1a / CAG）
- docstring でコドン選択の根拠を明記
- ベクター要素の段階的拡張が可能な設計に

---

## 5. report_generator.py

### 問題点
- `dna_sequence` キーが `app.py` から渡されるが `main.py` からは渡されない → KeyError リスク
- フォントサイズやマージンがハードコード散在
- ページ番号なし
- エラーハンドリングなし（ファイル書き込み失敗時など）

### 改善方針
- `results.get('dna_sequence', 'N/A')` で安全にアクセス
- レイアウト定数をクラス変数に集約
- ページ番号追加
- try-except による保護

---

## 6. 横断的な改善事項

| 項目 | 現状 | 改善後 |
|:---|:---|:---|
| docstring | なし | Google Style で全クラス・メソッドに追加 |
| 型ヒント | なし | 全関数の引数・返り値に付与 |
| logging | なし（print文） | logging モジュール使用 |
| 定数管理 | 各ファイルにハードコード | config.py に集約 |
| エラー処理 | なし | try-except + 意味あるエラーメッセージ |
| プレースホルダー表記 | 不明確 | `# TODO:` + docstring で設計意図を明記 |
