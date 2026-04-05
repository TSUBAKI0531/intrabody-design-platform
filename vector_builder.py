"""
vector_builder.py — レンチウイルス発現ベクターの自動構築

イントラボディ配列をヒト細胞で発現させるためのベクター DNA 配列を
設計し、GenBank 形式で出力する。コドン最適化（ヒト最頻コドン）、
プロモーター選択、Kozak 配列の自動挿入に対応。
"""

from __future__ import annotations

import logging
from enum import Enum

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class Promoter(Enum):
    """利用可能なプロモーター。"""

    EF1A = "EF1a"
    CMV = "CMV"
    CAG = "CAG"


# ---------------------------------------------------------------------------
# プロモーター配列（代表的な 5' 末端領域）
# ---------------------------------------------------------------------------
PROMOTER_SEQUENCES: dict[Promoter, str] = {
    Promoter.EF1A: (
        "GGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCG"
    ),
    Promoter.CMV: (
        "TAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAAC"
    ),
    Promoter.CAG: (
        "CTCGACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCCAT"
    ),
}

# Kozak 配列（翻訳開始効率を向上）
KOZAK_SEQUENCE: str = "GCCACC"

# ヒト最頻コドンテーブル
HUMAN_CODON_TABLE: dict[str, str] = {
    "A": "GCC", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC",
    "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "TCC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


class VectorDesigner:
    """レンチウイルス発現ベクターの設計と GenBank ファイル出力。

    Parameters:
        promoter: 使用するプロモーター（デフォルト: EF1a）
        add_kozak: Kozak 配列を CDS 上流に挿入するか
        add_stop_codon: 終止コドンを自動付与するか

    Example:
        >>> designer = VectorDesigner(promoter=Promoter.CMV)
        >>> designer.build_genbank("MyITB", "EVQLVESGG...", "output.gb")
    """

    def __init__(
        self,
        promoter: Promoter = Promoter.EF1A,
        add_kozak: bool = True,
        add_stop_codon: bool = True,
    ) -> None:
        self.promoter = promoter
        self.add_kozak = add_kozak
        self.add_stop_codon = add_stop_codon
        self.codon_table = HUMAN_CODON_TABLE

    def _translate_to_dna(self, aa_sequence: str) -> str:
        """アミノ酸配列をヒト最頻コドンで DNA に変換する。"""
        codons = [self.codon_table.get(aa, "NNN") for aa in aa_sequence]
        unknown = [aa for aa in aa_sequence if aa not in self.codon_table]
        if unknown:
            logger.warning("未知のアミノ酸が含まれています: %s", unknown)
        return "".join(codons)

    def build_genbank(
        self,
        name: str,
        aa_sequence: str,
        output_path: str,
        description: str = "Intrabody Expression Vector",
    ) -> str:
        """GenBank 形式のベクター配列ファイルを生成する。

        Parameters:
            name: コンストラクト名
            aa_sequence: イントラボディのアミノ酸配列
            output_path: 出力ファイルパス
            description: SeqRecord の説明文

        Returns:
            生成された DNA 配列全長
        """
        # --- DNA 配列の構築 ---
        promoter_seq = PROMOTER_SEQUENCES[self.promoter]
        kozak_seq = KOZAK_SEQUENCE if self.add_kozak else ""
        cds_dna = self._translate_to_dna(aa_sequence)

        # 開始コドン（ATG）が先頭にない場合は追加
        if not cds_dna.startswith("ATG"):
            cds_dna = "ATG" + cds_dna

        # 終止コドンの付与
        if self.add_stop_codon and not cds_dna.endswith(("TAA", "TAG", "TGA")):
            cds_dna += self.codon_table["*"]

        full_dna = promoter_seq + kozak_seq + cds_dna

        # --- SeqRecord の構築 ---
        record = SeqRecord(
            Seq(full_dna),
            id=name,
            name=name[:16],
            description=description,
        )
        record.annotations["molecule_type"] = "DNA"

        # Feature: プロモーター
        promoter_end = len(promoter_seq)
        record.features.append(
            SeqFeature(
                FeatureLocation(0, promoter_end),
                type="promoter",
                qualifiers={"label": self.promoter.value},
            )
        )

        # Feature: Kozak
        if self.add_kozak:
            kozak_start = promoter_end
            kozak_end = kozak_start + len(kozak_seq)
            record.features.append(
                SeqFeature(
                    FeatureLocation(kozak_start, kozak_end),
                    type="misc_feature",
                    qualifiers={"label": "Kozak", "note": "Translation initiation"},
                )
            )

        # Feature: CDS
        cds_start = promoter_end + len(kozak_seq)
        record.features.append(
            SeqFeature(
                FeatureLocation(cds_start, len(full_dna)),
                type="CDS",
                qualifiers={
                    "label": "Intrabody",
                    "codon_start": 1,
                    "translation": aa_sequence,
                },
            )
        )

        # --- 書き出し ---
        SeqIO.write(record, output_path, "genbank")
        logger.info(
            "GenBank ファイル生成: %s (プロモーター: %s, %d bp)",
            output_path,
            self.promoter.value,
            len(full_dna),
        )
        return full_dna
