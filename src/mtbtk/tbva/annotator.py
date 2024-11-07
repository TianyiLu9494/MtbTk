import pandas as pd
from intervaltree import IntervalTree
from .SNP_effect_annotator import SnpEffectAnnotator
from .lineage_identifier import LineageIdentifier
from .config import genes, domains, epitopes

class VariantAnnotator:
    def __init__(self, SNPdf):
        self.SNPdf = SNPdf
        self.genes = genes
        self.annotator = SnpEffectAnnotator()
        self.gene_intervals = self._create_gene_intervals()
        self.domains = domains
        self.domain_intervals = self._create_domain_intervals()
        self.epitopes = epitopes
        self.epitope_intervals= self._create_epitope_intervals()
        
    def _create_gene_intervals(self):
        gene_intervals = IntervalTree()
        for gene_name, gene_info in self.genes.items():
            if gene_info['biotype'] in ['protein_coding', 'pseudogene', 'other']:
                gene_intervals[gene_info['start']:gene_info['end']] = gene_name
        return gene_intervals
    
    def _reverse_complement(self,codon):
        # 定义碱基互补对
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # 计算反向互补序列
        reverse_complement_codon = ''.join(complement[base] for base in reversed(codon))
        
        return reverse_complement_codon
    
    def _is_in_same_codon(self, snp1, snp2):
        overlapping_genes = self.gene_intervals[snp1.POS] & self.gene_intervals[snp2.POS]
        for interval in overlapping_genes:
            gene_name = interval.data
            gene_info = self.genes[gene_name]
            start, end = gene_info['start'], gene_info['end']
            strand = gene_info['strand']
            
            if strand == '+':
                codon1_pos = (snp1.POS - start) // 3
                codon2_pos = (snp2.POS - start) // 3
                ref_codon = gene_info['sequence'][codon1_pos*3:codon1_pos*3 + 3]
            else:  # strand == '-'
                codon1_pos = (end - snp1.POS) // 3
                codon2_pos = (end - snp2.POS) // 3
                ref_codon = self._reverse_complement(gene_info['sequence'][codon1_pos*3:codon1_pos*3 + 3])
            if codon1_pos == codon2_pos:
                return [True, gene_info, ref_codon]
        return [False, None, None]
    
    def _create_domain_intervals(self):
        domain_intervals = IntervalTree()
        for domain in self.domains:
            domain_intervals[domain['genome_pos_start']:domain['genome_pos_end']] = domain
        return domain_intervals

    def _is_in_domain(self, snp_pos):
        overlapping_domains = self.domain_intervals[snp_pos]
        if overlapping_domains:
            return list(overlapping_domains)
        return None
    
    def _create_epitope_intervals(self):
        epitope_intervals = IntervalTree()
        for epitope in self.epitopes:
            epitope_intervals[epitope['start']:epitope['end']] = epitope
        return epitope_intervals
    
    def _is_in_epitope(self, snp_pos):
        overlapping_epitopes = self.epitope_intervals[snp_pos]
        if overlapping_epitopes:
            return list(overlapping_epitopes)
        return None
    
    def LineageAlleleIdentification(self):
        LeangeAlleles = LineageIdentifier(self.SNPdf)
        return pd.DataFrame(LeangeAlleles).sort_values(by='lineage')

    def SnpEffectAnnotation(self):
        AnnotatedSNP = []
        rows = self.SNPdf.itertuples(index=False)
        current_snp = next(rows, None)

        while current_snp is not None:
            next_snp = next(rows, None)
            if next_snp:
                in_same_codon, gene, ref_codon = self._is_in_same_codon(current_snp, next_snp)
            else:
                in_same_codon = False
            if in_same_codon:
                # 计算两个 SNP 之间的距离
                distance = next_snp.POS - current_snp.POS
                if distance > 1:
                    # 如果两个 SNP 之间有间隔，说明它们中间还有一个碱基未发生变异
                    pos = current_snp.POS
                    ref = ref_codon
                    alt = current_snp.ALT + ref_codon[1] + next_snp.ALT
                    current_snp = next_snp
                else:
                    pos = current_snp.POS
                    ref = current_snp.REF + next_snp.REF
                    alt = current_snp.ALT + next_snp.ALT
                    
                    # 检查是否还有第三个SNP在同一个codon
                    third_snp = next(rows, None)
                    if third_snp:
                        in_same_third, gene_2, ref_codon_third = self._is_in_same_codon(current_snp, third_snp)
                    else:
                        in_same_third = False
                    if in_same_third and (gene_2 == gene):
                        ref += third_snp.REF
                        alt += third_snp.ALT
                        current_snp = next(rows, None)
                    else:
                        current_snp = third_snp
            else:
                pos = current_snp.POS
                ref = current_snp.REF
                alt = current_snp.ALT
                current_snp = next_snp

            # 注释
            try:
                gene_name, variant, annotation, variant_impact = self.annotator.annotate(int(pos), ref, alt)
            except Exception as e:
                print(f"Error occurred at position: {pos}")
                print(f"Current SNP: POS={pos}, REF={ref}, ALT={alt}")
                print(f"Error details: {str(e)}")
                raise
            
            # 存储注释结果
            AnnotatedSNP.append({
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'GENE': gene_name,
                'VARIANT': variant,
                'VARIANT_TYPE': annotation,
                'IMPACT': variant_impact
            })
        self.AnnotatedSNP = pd.DataFrame(AnnotatedSNP)
        return self.AnnotatedSNP
    
    def EpitopeAnnotation(self):
        epitope_annotations = []
        missense_variants = self.AnnotatedSNP[self.AnnotatedSNP['VARIANT_TYPE'] == 'missense_variant']
        
        for _, variant in missense_variants.iterrows():
            epitopes = self._is_in_epitope(variant['POS'])
            if epitopes:
                for interval in epitopes:
                    annotation = variant.to_dict()  # 首先创建包含 variant 信息的字典
                    annotation.update(interval.data)  # 然后添加 epitope 信息
                    epitope_annotations.append(annotation)
        
        return pd.DataFrame(epitope_annotations)

    def DomainAnnotation(self):
        domain_annotations = []
        missense_variants = self.AnnotatedSNP[self.AnnotatedSNP['VARIANT_TYPE'] == 'missense_variant']
        
        for _, variant in missense_variants.iterrows():
            domains = self._is_in_domain(variant['POS'])
            if domains:
                for interval in domains:
                    annotation = variant.to_dict()  # 首先创建包含 variant 信息的字典
                    annotation.update(interval.data)  # 然后添加 domain 信息
                    domain_annotations.append(annotation)
        
        return pd.DataFrame(domain_annotations)
    
    def SnpEffectAnnotation_test(self):
        AnnotatedSNP = []
        rows = self.SNPdf.itertuples(index=False)
        current_snp = next(rows, None)
        previous_snps = []  # 用于存储之前的SNP信息

        while current_snp is not None:
            next_snp = next(rows, None)
            
            if next_snp and self._is_in_same_codon(current_snp, next_snp):
                # 如果在同一个codon，合并SNP
                pos = current_snp.POS
                ref = current_snp.REF + next_snp.REF
                alt = current_snp.ALT + next_snp.ALT
                
                # 检查是否还有第三个SNP在同一个codon
                third_snp = next(rows, None)
                if third_snp and self._is_in_same_codon(current_snp, third_snp):
                    ref += third_snp.REF
                    alt += third_snp.ALT
                    current_snp = next(rows, None)
                else:
                    current_snp = third_snp
            else:
                pos = current_snp.POS
                ref = current_snp.REF
                alt = current_snp.ALT
                current_snp = next_snp

            # 注释
            try:
                gene_name, variant, annotation, variant_impact = self.annotator.annotate(int(pos), ref, alt)
            except Exception as e:
                print(f"Error occurred at position: {pos}")
                print(f"Current SNP: POS={pos}, REF={ref}, ALT={alt}")
                
                # 打印之前的三个SNP（如果有的话）
                for i, prev_snp in enumerate(previous_snps[-3:], 1):
                    print(f"Previous SNP {i}: POS={prev_snp['POS']}, REF={prev_snp['REF']}, ALT={prev_snp['ALT']}")
                
                print(f"Error details: {str(e)}")
                raise  # 重新抛出异常
            
            # 存储注释结果
            annotated_snp = {
                'POS': pos,
                'REF': ref,
                'ALT': alt,
                'GENE': gene_name,
                'VARIANT': variant,
                'ANNOTATION': annotation,
                'IMPACT': variant_impact
            }
            AnnotatedSNP.append(annotated_snp)
            
            # 保存当前SNP信息到previous_snps列表
            previous_snps.append(annotated_snp)
            if len(previous_snps) > 3:
                previous_snps.pop(0)  # 只保留最近的3个SNP
        self.AnnotatedSNP = pd.DataFrame(AnnotatedSNP)
        return self.AnnotatedSNP