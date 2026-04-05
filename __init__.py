"""
src — Intrabody Design Platform コアモジュール

細胞内抗体（イントラボディ）の設計パイプラインを構成するモジュール群。

Modules:
    optimizer:        pI 最適化・親水化エンジン（FR 領域の残基置換）
    simulator:        結合シミュレーション・構造予測・特異性評価
    evaluator:        候補配列のバッチ評価・ランキング
    vector_builder:   GenBank 形式ベクター構築（コドン最適化対応）
    report_generator: PDF 技術レポートの自動生成
    config:           定数・閾値の一元管理
"""
