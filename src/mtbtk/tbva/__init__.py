from .vcf_parser import VCFParser  # 假设你的 VCFParser 类在 vcf_parser.py 文件中
from .lineage_identifier import LineageIdentifier
from .annotator import VariantAnnotator
from .SNP_effect_annotator import SnpEffectAnnotator
from .config import genes

__all__ = ['genes','VCFParser','VariantAnnotator','LineageIdentifier','SnpEffectAnnotator'] 