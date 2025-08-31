#!/usr/bin/env python3

import sys
from collections import defaultdict, Counter
import argparse

class GenomeDecoder:
    def __init__(self, filename):
        self.filename = filename
        self.snps = []
        self.chromosome_counts = defaultdict(int)
        self.genotype_counts = Counter()
        self.total_snps = 0
        self.no_calls = 0
        
    def parse_genome_file(self):
        print(f"Parsing genome file: {self.filename}")
        
        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('#') or not line:
                    continue
                
                if line.startswith('rsid'):
                    continue
                    
                parts = line.split('\t')
                if len(parts) == 4:
                    rsid, chromosome, position, genotype = parts
                    
                    self.snps.append({
                        'rsid': rsid,
                        'chromosome': chromosome,
                        'position': int(position) if position.isdigit() else position,
                        'genotype': genotype
                    })
                    
                    self.chromosome_counts[chromosome] += 1
                    self.genotype_counts[genotype] += 1
                    self.total_snps += 1
                    
                    if genotype == '--':
                        self.no_calls += 1
        
        print(f"Successfully parsed {self.total_snps:,} SNPs")
    
    def generate_summary(self):
        print("\n" + "="*60)
        print("GENOME SUMMARY STATISTICS")
        print("="*60)
        
        print(f"\nTotal SNPs analyzed: {self.total_snps:,}")
        print(f"Successfully genotyped: {self.total_snps - self.no_calls:,} ({100*(self.total_snps - self.no_calls)/self.total_snps:.2f}%)")
        print(f"No-calls (--): {self.no_calls:,} ({100*self.no_calls/self.total_snps:.2f}%)")
        
        print("\n" + "-"*40)
        print("CHROMOSOME DISTRIBUTION")
        print("-"*40)
        
        sorted_chromosomes = []
        for chrom in self.chromosome_counts.keys():
            if chrom.isdigit():
                sorted_chromosomes.append((int(chrom), chrom))
            else:
                sorted_chromosomes.append((100 if chrom == 'X' else 101 if chrom == 'Y' else 102, chrom))
        
        sorted_chromosomes.sort()
        
        for _, chrom in sorted_chromosomes:
            count = self.chromosome_counts[chrom]
            percentage = (count / self.total_snps) * 100
            bar = '█' * int(percentage * 2)
            print(f"Chr {chrom:>2}: {count:>7,} SNPs ({percentage:>5.2f}%) {bar}")
        
        print("\n" + "-"*40)
        print("GENOTYPE DISTRIBUTION")
        print("-"*40)
        
        homozygous = 0
        heterozygous = 0
        
        for genotype, count in self.genotype_counts.items():
            if genotype != '--':
                if len(genotype) == 2:
                    if genotype[0] == genotype[1]:
                        homozygous += count
                    else:
                        heterozygous += count
        
        print(f"Homozygous variants: {homozygous:,} ({100*homozygous/(homozygous+heterozygous):.2f}%)")
        print(f"Heterozygous variants: {heterozygous:,} ({100*heterozygous/(homozygous+heterozygous):.2f}%)")
        
        print("\n" + "-"*40)
        print("TOP 10 MOST COMMON GENOTYPES")
        print("-"*40)
        
        for genotype, count in self.genotype_counts.most_common(10):
            percentage = (count / self.total_snps) * 100
            print(f"{genotype:>3}: {count:>7,} ({percentage:>5.2f}%)")
    
    def analyze_specific_snps(self):
        print("\n" + "-"*40)
        print("GENETIC VARIANT ANALYSIS WITH INTERPRETATIONS")
        print("-"*40)
        
        interesting_snps = {
            'rs1815739': 'ACTN3 gene (muscle performance)',
            'rs1800497': 'DRD2 gene (dopamine receptor)',
            'rs4680': 'COMT gene (dopamine metabolism)',
            'rs1799971': 'OPRM1 gene (opioid receptor)',
            'rs9939609': 'FTO gene (obesity risk)',
            'rs7903146': 'TCF7L2 gene (diabetes risk)',
            'rs429358': 'APOE gene variant (Alzheimer\'s risk)',
            'rs7412': 'APOE gene variant (Alzheimer\'s risk)',
        }
        
        genotype_interpretations = {
            'rs1815739': {
                'CC': 'RR genotype - Enhanced power/sprint performance, fast-twitch muscle fibers',
                'CT': 'RX genotype - Mixed muscle fiber type, intermediate performance',
                'TT': 'XX genotype - Enhanced endurance, higher injury risk, slower sprint performance'
            },
            'rs1800497': {
                'GG': 'Normal dopamine receptor (DRD2) density',
                'AG': 'Slightly reduced dopamine receptor density',
                'AA': 'Reduced dopamine receptor density, may affect reward processing'
            },
            'rs4680': {
                'GG': 'Val/Val - High COMT activity, faster dopamine breakdown ("warrior" variant)',
                'AG': 'Val/Met - Intermediate COMT activity',
                'AA': 'Met/Met - Low COMT activity (~75% reduced), slower dopamine breakdown, higher prefrontal dopamine ("worrier" variant)'
            },
            'rs1799971': {
                'AA': 'Normal opioid receptor (OPRM1) function',
                'AG': 'Altered opioid receptor function, may affect pain sensitivity',
                'GG': 'Significantly altered opioid receptor function, higher pain medication requirements'
            },
            'rs9939609': {
                'TT': 'Lower obesity risk, protective variant',
                'AT': 'Intermediate obesity risk',
                'AA': 'Higher obesity risk, increased tendency for high sugar/fat intake'
            },
            'rs7903146': {
                'CC': 'Lower type 2 diabetes risk',
                'CT': 'Moderate type 2 diabetes risk (~1.4x)',
                'TT': 'Higher type 2 diabetes risk (~2x)'
            }
        }
        
        print("\nSearching for notable SNPs in your genome:")
        found_snps = {}
        apoe_snps = {}
        
        for snp in self.snps:
            if snp['rsid'] in interesting_snps:
                found_snps[snp['rsid']] = snp
                if snp['rsid'] in ['rs429358', 'rs7412']:
                    apoe_snps[snp['rsid']] = snp['genotype']
        
        if not found_snps:
            print("None of the sample notable SNPs were found in this dataset.")
            return
        
        for rsid, snp in found_snps.items():
            if rsid not in ['rs429358', 'rs7412']:
                print(f"\n{rsid}: {interesting_snps[rsid]}")
                print(f"  Location: Chromosome {snp['chromosome']}, Position {snp['position']:,}")
                print(f"  Your genotype: {snp['genotype']}")
                
                if rsid in genotype_interpretations and snp['genotype'] in genotype_interpretations[rsid]:
                    print(f"  Interpretation: {genotype_interpretations[rsid][snp['genotype']]}")
        
        if 'rs429358' in apoe_snps and 'rs7412' in apoe_snps:
            print(f"\nAPOE Status (combined rs429358 and rs7412):")
            rs429358_gt = apoe_snps['rs429358']
            rs7412_gt = apoe_snps['rs7412']
            
            print(f"  rs429358: {rs429358_gt}, rs7412: {rs7412_gt}")
            
            if rs429358_gt == 'TT' and rs7412_gt == 'TT':
                print(f"  APOE type: ε2/ε2 - Lower Alzheimer's risk (protective)")
            elif rs429358_gt == 'TT' and rs7412_gt == 'CT':
                print(f"  APOE type: ε2/ε3 - Lower than average Alzheimer's risk")
            elif rs429358_gt == 'TT' and rs7412_gt == 'CC':
                print(f"  APOE type: ε3/ε3 - Average Alzheimer's risk (most common)")
            elif rs429358_gt == 'CT' and rs7412_gt == 'CC':
                print(f"  APOE type: ε3/ε4 - Moderately increased Alzheimer's risk (~3x)")
            elif rs429358_gt == 'CC' and rs7412_gt == 'CC':
                print(f"  APOE type: ε4/ε4 - Significantly increased Alzheimer's risk (~12-15x)")
            elif rs429358_gt == 'CT' and rs7412_gt == 'CT':
                print(f"  APOE type: ε2/ε4 - Variable Alzheimer's risk")
            else:
                print(f"  APOE type: Unable to determine from these genotypes")
        elif 'rs429358' in apoe_snps or 'rs7412' in apoe_snps:
            print(f"\nAPOE Status: Incomplete (need both rs429358 and rs7412 to determine)")
    
    def export_statistics(self, output_file="genome_stats.txt"):
        with open(output_file, 'w') as f:
            f.write(f"Genome Statistics Report\n")
            f.write(f"{'='*40}\n")
            f.write(f"File: {self.filename}\n")
            f.write(f"Total SNPs: {self.total_snps:,}\n")
            f.write(f"Successfully genotyped: {self.total_snps - self.no_calls:,}\n")
            f.write(f"No-calls: {self.no_calls:,}\n\n")
            
            f.write(f"Chromosome Distribution:\n")
            for chrom, count in sorted(self.chromosome_counts.items()):
                f.write(f"  Chr {chrom}: {count:,} SNPs\n")
        
        print(f"\nStatistics exported to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Decode and analyze 23andMe genome data')
    parser.add_argument('--file', type=str, default='genome.txt',
                       help='Path to the genome file')
    parser.add_argument('--export', action='store_true', help='Export statistics to file')
    
    args = parser.parse_args()
    
    decoder = GenomeDecoder(args.file)
    decoder.parse_genome_file()
    decoder.generate_summary()
    decoder.analyze_specific_snps()
    
    if args.export:
        decoder.export_statistics()

if __name__ == "__main__":
    main()