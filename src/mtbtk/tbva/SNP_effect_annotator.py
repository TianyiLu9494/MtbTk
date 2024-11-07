import json
import os
from .config import genes
class SnpEffectAnnotator:
    def __init__(self):
        self.genes = genes
        self.positions = {}
        for gene_name,gene_property in self.genes.items():
            self.positions[gene_property['start']] = f'{gene_name}-start'
            self.positions[gene_property['end']] = f'{gene_name}-end'
        self.sort_positions = sorted(self.positions.keys())
        self.variant_impact = {
            "missense_variant": "MODERATE",
            "stop_gained": "HIGH",
            "loss_of_function;stop_gained": "HIGH",
            "synonymous_variant": "LOW",
            "stop_lost": "HIGH",
            "stop_retained_variant": "LOW",
            "upstream_gene_variant": "MODIFIER",
            "start_lost": "HIGH",
            "loss_of_function;start_lost": "HIGH",
            "initiator_codon_variant": "LOW",
            "start_retained_variant": "LOW",
            "downstream_gene_variant": "MODIFIER",
            "non_coding_transcript_exon_variant": "LOW",
        }
        self.aa_to_one_letter = {
            "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
            "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G",
            "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
            "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
            "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
            "Ter": "*"
        }
        self.codon_table = {
            "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
            "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
            "TAT": "Tyr", "TAC": "Tyr", "TAA": "*", "TAG": "*",
            "TGT": "Cys", "TGC": "Cys", "TGA": "*", "TGG": "Trp",
            "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
            "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
            "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
            "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
            "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
            "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
            "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
            "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
            "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
            "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
            "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
            "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
        }
        self.initiation_codons = ["TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"]
        self.LoF_gene = ['fbiA','fbiB','fbiC','fgd1','mmpL5','tlyA','pncA','eis','pepQ','Rv2983','ddn','ethA','gid','katG']

    def find_closest_position(self,target):
        '''
        输入基因组中所有基因的start和end position（以list形式排好序）；
        和一个目标位置；
        返回离目标位置最近的两个基因的start和end位置。
        '''
        if not self.sort_positions:
            return None, None  # 返回None如果数组为空

        # 边界情况处理
        if target <= self.sort_positions[0] or target >= self.sort_positions[-1]:
            closest_left = self.sort_positions[-1]
            closest_right = self.sort_positions[0]  # 目标数大于等于数组最后一个元素
            second_left = self.sort_positions[-2]
            second_right = self.sort_positions[1]
            return closest_left, closest_right, second_left, second_right

        # 二分查找
        left, right = 0, len(self.sort_positions) - 1
        while left <= right:
            mid = left + (right - left) // 2
            if target == self.sort_positions[mid]:
                break  # 找到目标，跳出循环
            elif target < self.sort_positions[mid]:
                right = mid - 1
            else:
                left = mid + 1

        # 找到最接近的两个数
        closest_left = self.sort_positions[right] if right >= 0 else None
        closest_right = self.sort_positions[left] if left < len(self.sort_positions) else 1
        # 找到其次近的两个数
        second_left = self.sort_positions[right - 1]
        second_right = self.sort_positions[left + 1] if left < len(self.sort_positions) else (left + 1 if left + 1 < len(self.sort_positions) else 2)
        return closest_left, closest_right, second_left, second_right

    @staticmethod
    def get_reverse_complementary(bases):
        base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        complementary_bases = []
        for base in bases:
            if base in base_dict:
                complementary_bases.append(base_dict[base])
            else:
                raise ValueError('Invalid base')
        return ''.join(complementary_bases[::-1])

    @staticmethod
    def get_alt_codon(ref_codon, pos_in_codon, alt):
        # 确保位置在正确的范围内
        if pos_in_codon not in [0, 1, 2]:
            raise ValueError("位置必须是0,1或2")
    
        # 确保替换的长度不会超过密码子的长度
        if len(alt) + pos_in_codon > 3:
            raise ValueError(f"替换的碱基超过了一个codon的长度。input:{ref_codon, pos_in_codon, alt}")
    
        # 替换碱基
        new_codon = ref_codon[:pos_in_codon] + alt + ref_codon[pos_in_codon + len(alt):]
    
        return new_codon

    def generate_aminoacid_variant(self,gene, pos, ref, alt):
        strand = gene['strand']
        if strand == '-':
            ref = self.get_reverse_complementary(ref)
            alt = self.get_reverse_complementary(alt)
            pos_in_gene = gene['end'] - pos - len(ref) + 1
        else:
            pos_in_gene = pos - gene['start']
        ref_codon_pos = pos_in_gene // 3
        pos_in_codon = pos_in_gene % 3
        ref_codon = gene['sequence'][ref_codon_pos*3:ref_codon_pos*3+3]
        # check if the nucleotide matches the reference
        if gene['sequence'][ref_codon_pos*3+pos_in_codon:ref_codon_pos*3+pos_in_codon+len(ref)] != ref:
            raise ValueError(f'The reference nucleotide {ref}({pos}) does not match the genome sequence at {gene["ID"]} {gene["sequence"][ref_codon_pos*3+pos_in_codon:ref_codon_pos*3+pos_in_codon+len(ref)]}')
        # Change the ref to alt in the codon at pos_in_codon
        alt_codon = self.get_alt_codon(ref_codon, pos_in_codon, alt)
        aa_ref = self.codon_table[ref_codon]
        aa_alt = self.codon_table[alt_codon]
        if ref_codon_pos == 0:
            if alt_codon in self.initiation_codons or aa_ref == aa_alt:
                return ['start_retained_variant',[ref_codon_pos + 1,pos_in_gene + 1]]
            else:
                return ['start_lost',[ref_codon_pos + 1,pos_in_gene + 1]]
        elif aa_ref == aa_alt:
            return ['synonymous_variant',[pos_in_gene + 1, ref_codon_pos + 1, aa_ref]]
        else:
            return ['missense_variant',[ref_codon_pos + 1, aa_ref, aa_alt]]

    def coding_gene_annotator(self,gene_name, pos, ref, alt):
        gene = self.genes[gene_name]
        variant_type, properties = self.generate_aminoacid_variant(gene, pos, ref, alt)
        if variant_type == 'missense_variant':
            aa_pod,aa_ref,aa_alt = properties
            if aa_ref == "*" and aa_alt != "*":
                annotation = 'stop_lost'
                hgvs = f'p.Ter{aa_pod}{aa_alt}ext*?'
            elif aa_ref != "*" and aa_alt == "*":
                if gene_name in self.LoF_gene:
                    annotation = 'loss_of_function;stop_gained'
                    hgvs = f'LoF'
                else:
                    annotation = 'stop_gained'
                    hgvs = f'p.{aa_ref}{aa_pod}{aa_alt}'
            else:
                annotation = 'missense_variant'
                hgvs = f'p.{aa_ref}{aa_pod}{aa_alt}'
        elif variant_type == 'synonymous_variant':
            nuc_pos, codon_pos, aa = properties
            if aa == "*":
                annotation = 'stop_retained_variant'
            else:
                annotation = 'synonymous_variant'
            if gene['strand'] == '-':
                ref = self.get_reverse_complementary(ref)
                alt = self.get_reverse_complementary(alt)
            hgvs = f'c.{nuc_pos}{ref}>{alt}'
        elif variant_type == 'start_retained_variant':
            aa_pos, nuc_pos = properties
            annotation = variant_type
            if gene['strand'] == '-':
                ref = self.get_reverse_complementary(ref)
                alt = self.get_reverse_complementary(alt)
            hgvs = f'c.{nuc_pos}{ref}>{alt}'
        elif variant_type == 'start_lost':
            if gene_name in self.LoF_gene:
                annotation = 'loss_of_function;start_lost'
                hgvs = f'LoF'
            else:
                annotation = variant_type
                hgvs = f'p.Met1?'
        variant = f'{gene_name}_{hgvs}'
        variant_impact = self.variant_impact[annotation]
        return gene_name, variant, annotation, variant_impact

    def rna_annotator(self,gene_name, pos, ref, alt):
        """
        本函数用于处理非蛋白编码基因的SNP
        :param gene_name: 
        :param pos: 
        :param ref: 
        :param alt: 
        :return: 
        """
        gene = self.genes[gene_name]
        if gene['strand'] == '-':
            ref = self.get_reverse_complementary(ref)
            alt = self.get_reverse_complementary(alt)
            pos_in_gene = gene['end'] - pos + 1
        else:
            pos_in_gene = pos - gene['start'] + 1
        annotation = 'non_coding_transcript_exon_variant'
        variant = f'{gene_name}_n.{pos_in_gene}{ref}>{alt}'
        variant_impact = self.variant_impact[annotation]
        return gene_name, variant, annotation, variant_impact

    def upstream_gene_annotator(self,gene_name,pos,ref,alt):
        gene = self.genes[gene_name]
        if gene['strand'] == '+':
            if pos < gene['start']:
                pos_in_gene = pos - gene['start']
            else:
                # 应对特殊情况，考虑到基因组是环状所以，也就是SNP在最后一个基因的右侧在第一个基因的左侧时，
                # 这时候pos到gene['start']的上游多远就要用pos减去（基因组的长度加第一个基因的start）
                pos_in_gene = pos - (4411532 + 1)
        else:
            ref = self.get_reverse_complementary(ref)
            alt = self.get_reverse_complementary(alt)
            pos_in_gene = gene['end'] - pos
        annotation = 'upstream_gene_variant'
        variant = f'{gene_name}_c.{pos_in_gene}{ref}>{alt}'
        if gene['biotype'] in ['rRNA', 'tRNA', 'ncRNA','misc_RNA']:
            variant = variant.replace('c.','n.')
        variant_impact = self.variant_impact[annotation]
        return gene_name, variant, annotation, variant_impact

    def downstream_gene_annotator(self,gene_name,pos,ref,alt):
        gene = self.genes[gene_name]
        if gene['strand'] == '+':
            pos_in_gene = pos - gene['end']
        else:
            ref = self.get_reverse_complementary(ref)
            alt = self.get_reverse_complementary(alt)
            pos_in_gene = gene['start'] - pos
        annotation = 'downstream_gene_variant'
        variant = f'{gene_name}_c.*{pos_in_gene}{ref}>{alt}'
        if gene['biotype'] in ['rRNA', 'tRNA', 'ncRNA','misc_RNA']:
            variant = variant.replace('c.','n.')
        variant_impact = self.variant_impact[annotation]
        return gene_name, variant, annotation, variant_impact


    def irg_annotator(self,left_gene_name,right_gene_name, pos, ref, alt):
        """
        本函数用于处理两个基因之间的non-coding区域,主要根据两个基因的biotype和方向来判断SNP是在哪个基因上游/下游
        :param left_gene_name: 
        :param right_gene_name: 
        :param pos: 
        :param ref: 
        :param alt: 
        :return: 
        """
        left_gene = self.genes[left_gene_name]
        right_gene = self.genes[right_gene_name]
        # 首先判断两个基因的方向如果一致，那么我们选择SNP作为上游调控
        if left_gene['strand'] == '+' and right_gene['strand'] == '+':
            return self.upstream_gene_annotator(right_gene_name,pos,ref,alt)
        elif left_gene['strand'] == '-' and right_gene['strand'] == '-':
            return self.upstream_gene_annotator(left_gene_name,pos,ref,alt)
        # 如果两个基因的方向不一致，那么我们首先根据基因的优先级选择
        biotype_priority = self.compare_gene_biotype(left_gene_name,right_gene_name)
        name_priority = self.compare_gene_name(left_gene_name,right_gene_name)
        if biotype_priority != 'equal':
            priority_gene = biotype_priority
        elif name_priority != 'equal':
            priority_gene = name_priority
        else:
            priority_gene = None
        # 判断priority_gene是否为None，如果存在则根据这个基因的方向来判断SNP是在哪个基因上游/下游
        if priority_gene is not None:
            if priority_gene == left_gene_name:
                if self.genes[left_gene_name]['strand'] == '+':
                    return self.downstream_gene_annotator(left_gene_name,pos,ref,alt)
                else:
                    return self.upstream_gene_annotator(left_gene_name,pos,ref,alt)
            else:
                if self.genes[right_gene_name]['strand'] == '+':
                    return self.upstream_gene_annotator(right_gene_name,pos,ref,alt)
                else:
                    return self.downstream_gene_annotator(right_gene_name,pos,ref,alt)
        else:
            # 方向为一正一负时为下游
            if left_gene['strand'] == "+" and right_gene['strand'] == "-":
                closest_gene = left_gene_name if pos - left_gene['end'] < right_gene['start'] - pos else right_gene_name
                return self.downstream_gene_annotator(closest_gene,pos,ref,alt)
            # 方向为一负一正时为上游
            elif left_gene['strand'] == "-" and right_gene['strand'] == "+":
                closest_gene = left_gene_name if pos - left_gene['end'] < right_gene['start'] - pos else right_gene_name
                return self.upstream_gene_annotator(closest_gene,pos,ref,alt)

    def over_lapping_gene_annotator(self,left_gene_name,left_side,right_gene_name,right_side,pos,ref,alt):
        """
        本函数用于处理两个基因出现overlap的情况，主要根据SNP在两个基因上的情况。
        当恰好在两个基因overlapping的区域则更具优先级选择一个基因进行注释
        :param left_gene_name: 
        :param left_side: 
        :param right_gene_name: 
        :param right_side: 
        :param pos: 
        :param ref: 
        :param alt: 
        :return: 
        """
        if left_side ==  right_side:
            if left_side == 'start':
                return self.coding_gene_annotator(left_gene_name,pos,ref,alt)
            else:
                return self.coding_gene_annotator(right_gene_name,pos,ref,alt)
        elif left_side == 'start' and right_side == 'end':
            priority_gene = self.get_priority_gene(left_gene_name,right_gene_name)
            return self.coding_gene_annotator(priority_gene,pos,ref,alt)

    def compare_gene_biotype(self,left_gene_name,right_gene_name):
        """
        比较gene_biotype的函数 根据这个顺序['protein_coding','rRNA', 'tRNA', 'ncRNA','misc_RNA', 'pseudogene', 'other']
        :param left_gene_name: 
        :param right_gene_name: 
        :return: 
        """
        left_gene_biotype = self.genes[left_gene_name]['biotype']
        right_gene_biotype = self.genes[right_gene_name]['biotype']
        biotype_order = ['protein_coding','rRNA', 'tRNA', 'ncRNA','misc_RNA', 'pseudogene', 'other']
        if biotype_order.index(left_gene_biotype) < biotype_order.index(right_gene_biotype):
            return left_gene_name
        elif biotype_order.index(left_gene_biotype) > biotype_order.index(right_gene_biotype):
            return right_gene_name
        else:
            return 'equal'

    def compare_gene_name(self,left_gene_name,right_gene_name):
        """
        比较gene_name的函数，根据gene_name的字母顺序来判断
        :param left_gene_name: 
        :param right_gene_name: 
        :return: 
        """
        if left_gene_name.startswith('Rv') and not right_gene_name.startswith('Rv'):
            return right_gene_name
        elif not left_gene_name.startswith('Rv') and right_gene_name.startswith('Rv'):
            return left_gene_name
        else:
            return 'equal'

    def get_priority_gene(self,left_gene_name,right_gene_name):
        """
        本函数用于获取两个基因中优先级高的基因
        :param left_gene_name: 
        :param right_gene_name: 
        :return: 
        """
        #首先判断两个基因的biotype
        if self.compare_gene_biotype(left_gene_name,right_gene_name) != 'equal':
            return self.compare_gene_biotype(left_gene_name,right_gene_name)
        #如果两个基因的biotype相同，那么我们优先选择基因名不是Rv的基因
        if self.compare_gene_name(left_gene_name,right_gene_name) != 'equal':
            return self.compare_gene_name(left_gene_name,right_gene_name)
        #最后比较连个基因的长度，选择长度长的基因
        if self.genes[left_gene_name]['end'] - self.genes[left_gene_name]['start'] > self.genes[right_gene_name]['end'] - self.genes[right_gene_name]['start']:
            return left_gene_name
        else:
            return right_gene_name


    def annotate(self,pos,ref,alt):
        """
        本函数的是这个class最终的接口，主要负责第一步的判断，来调用不同的函数进行注释。
        首先根据pos找到最近的两个基因位置，然后根据这两个基因位置的信息，SNP是在一个基因内，两个overlap的基因上，还是在两个基因之间的non-coding区域
        """
        if pos in self.sort_positions:
            gene_name = self.positions[pos].split('-')[0]
            return self.coding_gene_annotator(gene_name,pos,ref,alt)

        # 找到两个位置对应的基因position
        left,right,second_left,second_right = self.find_closest_position(pos)
        try:
            left_gene_name, left_side = self.positions[left].split('-')
            right_gene_name, right_side = self.positions[right].split('-')
        except:
            print(f'left:{left};{self.positions[left].split("-")},right:{right};{self.positions[right].split("-")}')
            raise ValueError(f'Position not found in the gene positions: {pos}')
        #基因名相同说明在一个基因内
        if left_gene_name == right_gene_name:
            # 正常左边是gene的start右侧是gene的end
            if left_side == 'start':
                if self.genes[left_gene_name]['biotype'] == 'protein_coding' or self.genes[left_gene_name]['biotype'] == 'pseudogene':
                    return self.coding_gene_annotator(left_gene_name,pos,ref,alt)
                else:
                    return self.rna_annotator(left_gene_name,pos,ref,alt)
            else:
                raise ValueError(f'The closest gene is not the start of the gene:left:{left} as {left_gene_name}, right:{right} as {right_gene_name}. Position : {pos}, ref :{ref}, alt : {alt}')
        # 说明现在两基因的间区
        elif left_side == 'end' and right_side == 'start':
            # Special case when second_left and second_right are the same gene
            if self.positions[second_left].split('-')[0] == self.positions[second_right].split('-')[0]:
                gene_name = self.positions[second_left].split('-')[0]
                if self.positions[second_left].split('-')[1] == 'start':
                    if self.genes[gene_name]['biotype'] == 'protein_coding' or self.genes[gene_name]['biotype'] == 'pseudogene':
                        return self.coding_gene_annotator(gene_name,pos,ref,alt)
                    else:
                        return self.rna_annotator(gene_name,pos,ref,alt)
                else:
                    raise ValueError(f'The closest gene is not the start of the gene:left:{second_left} right:{second_right} as {gene_name}. Position : {pos}, ref :{ref}, alt : {alt}')
            else:
                return self.irg_annotator(left_gene_name,right_gene_name,pos,ref,alt)
        # 基因名不同且不在两个基因之间，说明出现两个基因overlapping的情况
        else:
            return self.over_lapping_gene_annotator(left_gene_name,left_side,right_gene_name,right_side,pos,ref,alt)