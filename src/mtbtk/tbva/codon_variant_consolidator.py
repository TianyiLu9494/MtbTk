from mtbtk import DATA_DIR
import pandas as pd
def identify_snp_groups(vcf_filename, cds_info):
    """
    Identifies groups of SNPs that are located in the same codon based on the provided VCF file and CDS information.

    Args:
        vcf_filename (str): The path to the VCF file containing SNP information.
        cds_info (list): A list of tuples containing CDS information. Each tuple should contain the following elements:
                         - cds (str): The CDS identifier.
                         - start (int): The start position of the CDS.
                         - end (int): The end position of the CDS.
                         - strand (str): The strand of the CDS.
                         - seq (str): The sequence of the CDS.

    Returns:
        tuple: A tuple containing two elements:
               - snp_groups (dict): A dictionary where the keys are the positions of the first SNP in each group, and the values
                                   are lists of tuples representing the SNPs in the group. Each tuple contains the following elements:
                                   - pos (int): The position of the SNP.
                                   - ref (str): The reference allele.
                                   - alt (str): The alternate allele.
               - merge_positions (set): A set containing the positions of the SNPs that need to be merged.

    """
    snp_groups = {}
    merge_positions = set()  # Record the positions of SNPs that need to be merged
    prev_snps = []
    with open(vcf_filename, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            pos, ref, alt = int(parts[1]), parts[3], parts[4]

            # Check if the current SNP is in the same codon as the previous SNPs
            snp_in_same_codon = False
            for cds, start, end, strand, seq in cds_info:
                if all(start <= p[0] <= end for p in prev_snps + [(pos, ref, alt)]):
                    pos_codons = [(p[0] - start) // 3 for p in prev_snps + [(pos, ref, alt)]]
                    if len(set(pos_codons)) == 1:
                        snp_in_same_codon = True
                        break

            if snp_in_same_codon:
                prev_snps.append((pos, ref, alt))
                if len(prev_snps) >= 2:
                    snp_groups[prev_snps[0][0]] = prev_snps.copy()
                    merge_positions.update([snp[0] for snp in prev_snps[1::]])
            else:
                prev_snps = [(pos, ref, alt)]

    return snp_groups, merge_positions

def merge_and_output_snps(vcf_filename, snp_groups, merge_positions, output_filename):
    """
    Merge SNPs in the same codon and output the modified VCF file.

    Args:
        vcf_filename (str): The path to the input VCF file.
        snp_groups (dict): A dictionary containing SNP groups where the key is the position and the value is a list of SNPs.
        merge_positions (list): A list of positions to merge SNPs.
        output_filename (str): The path to the output file.

    Returns:
        None
    """
    with open(vcf_filename, 'r') as file_in, open(output_filename, 'w') as file_out:
        for line in file_in:
            if line.startswith('#'):
                file_out.write(line)
                continue
            parts = line.strip().split('\t')
            pos = int(parts[1])

            if pos in snp_groups:
                # Merge SNPs
                snps = snp_groups[pos]
                ref = ''.join([snp[1] for snp in snps])
                alt = ''.join([snp[2] for snp in snps])
                parts[3] = ref
                parts[4] = alt
                file_out.write('\t'.join(parts) + '\n')
            elif pos not in merge_positions:
                file_out.write(line)

def main(vcf_file, output_vcf):
    cds_info_csv = './cds_info.csv'  # Define the path to cds_info.csv here
    cds_info = pd.read_csv(cds_info_csv).values.tolist()
    snp_groups, merge_positions = identify_snp_groups(vcf_file, cds_info)
    merge_and_output_snps(vcf_file, snp_groups, merge_positions, output_vcf)
    print("VCF file modification completed.")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description='''
#     This script processes Variant Call Format (VCF) files to merge Single Nucleotide Polymorphisms (SNPs) that occur in the same codon within a Coding DNA Sequence (CDS) region.

#     The script was developed to address a specific issue encountered when annotating VCF files using the SnpEff tool. SnpEff treats SNPs that occur in the same codon as independent events, even when they occur simultaneously and affect the same codon. This script merges such SNPs into a single entry, allowing SnpEff to correctly annotate them.

#     In addition to the original VCF file, the script requires a reference genome in FASTA format and a General Feature Format (GFF) file for the reference genome. The script extracts all CDS regions from the GFF file. It then iterates over the VCF file, checking whether each SNP falls within a CDS region. For SNPs within a CDS region, the script checks whether there are any adjacent SNPs within a range of 3 base pairs. It then checks whether these adjacent SNPs occur in the same codon, by subtracting the start position of the codon and dividing by 3. If the division results in the same number, the SNPs are considered to be in the same codon and are merged.

#     The script outputs a modified VCF file with merged SNPs.

#     Usage:
#         python Merge_SNP_in_same_codon.py -i <input_vcf> -o <output_vcf>

#     Arguments:
#         -i, --input_vcf: Path to the input VCF file.
#         -o, --output_vcf: Path to the output VCF file.
#     ''')
#     parser.add_argument('-i', '--input_vcf', type=str, help='Input VCF file path')
#     parser.add_argument('-o', '--output_vcf', type=str, help='Output VCF file path')
#     args = parser.parse_args()
#     main(args.input_vcf, args.output_vcf)