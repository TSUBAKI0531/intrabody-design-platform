from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

class VectorDesigner:
    def __init__(self):
        self.codon_table = {
            'A': 'GCC', 'R': 'CGC', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
            'Q': 'CAG', 'E': 'GAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
            'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTC', 'P': 'CCC',
            'S': 'TCC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTG', '*': 'TGA'
        }

    def build_genbank(self, name, aa_sequence, output_path):
        dna_body = "".join([self.codon_table.get(aa, "NNN") for aa in aa_sequence])
        
        promoter = "GGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCG"
        v5_tag = "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG"
        
        full_dna = promoter + dna_body + v5_tag
        record = SeqRecord(Seq(full_dna), id=name, name=name[:16], description="Intrabody Expression Vector")
        
        # --- ここを追加：GenBank出力に必須の情報 ---
        record.annotations["molecule_type"] = "DNA"
        # ---------------------------------------
        
        record.features.append(SeqFeature(FeatureLocation(0, len(promoter)), type="promoter", qualifiers={"label": "EF1a"}))
        record.features.append(SeqFeature(FeatureLocation(len(promoter), len(promoter)+len(dna_body)), type="CDS", qualifiers={"label": "Intrabody"}))
        
        # 保存
        SeqIO.write(record, output_path, "genbank")
        return output_path