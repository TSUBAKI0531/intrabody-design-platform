# Intrabody Development Platform 🧬

**Intrabody Development Platform** は、細胞内抗体（イントラボディ）のデザインからベクター構築、結合シミュレーションまでを自動化する統合ツールキットです。

## 🌟 主な機能
- **Intracellular Optimization**: 細胞内（還元環境）での凝集を防ぐための $pI$ 調整と親水化置換。
- **Batch Evaluation**: 複数の候補配列から、溶解度・電荷・発現効率に基づいたエリート配列の自動選別。
- **Binding Simulation**: AlphaFold2 等の予測に基づいた結合自由エネルギー ($\Delta G$) の評価。
- **Vector Construction**: ヒト細胞用コドン最適化を施したレンチウイルスベクター配列（GenBank形式）の自動生成。

## 🚀 インストール
```bash
git clone [https://github.com/YourUsername/intrabody-design-platform.git](https://github.com/YourUsername/intrabody-design-platform.git)
cd intrabody-design-platform
pip install -r requirements.txt
🛠 使い方
main.py に抗体候補配列を入力し、実行します。

Python
python main.py
📊 出力結果
results/ フォルダに最適化された GenBankファイル が生成されます。

同時に、ベクターマップ図（PNG）も出力されます。

📜 ライセンス
MIT License


---

## 3. VS Code から GitHub への公開手順

ターミナルを使わなくても、VS Code の機能だけで完結できます。

1.  **Gitのリポジトリ初期化:**
    * VS Codeの左メニューにある「ソース管理」アイコン（枝分かれしたアイコン）をクリック。
    * 「リポジトリを初期化する」を選択。
2.  **ファイルのステージングとコミット:**
    * 変更されたファイル一覧の横にある「+」ボタンを押し、すべてのファイルを「ステージング」します。
    * 上のメッセージ欄に `Initial commit` と入力し、「コミット」ボタンを押します。
3.  **GitHubへ公開:**
    * 「ブランチを公開」または「GitHubに公開」ボタンをクリック。
    * `Public`（公開）か `Private`（非公開）かを選択。
    * 公開が完了すると、ブラウザで自分のGitHubリポジトリを確認できます。



---

## 4. .gitignore の設定（重要）

Python特有のキャッシュや、一時的な結果ファイルをアップロードしないために、`.gitignore` ファイルを作成し、以下を記述してください。

```text
__pycache__/
*.pyc
results/*
!results/.gitkeep
.env
venv/
.DS_Store