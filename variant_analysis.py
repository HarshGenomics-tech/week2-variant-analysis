import pandas as pd
from scipy import stats
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv("clinvar.vcf", comment = "#", sep  = "\t", header = None)
print(df.head())
print(len(df))
def extract_clnsig(info_string):
    if "CLNSIG=" in info_string:
        return info_string.split("CLNSIG=")[1].split(";")[0]
    return None
df['clnsig'] = df[7].apply(extract_clnsig)
print(df['clnsig'].value_counts())

chrom_sizes_hg38 = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415
}

total_variants = len(df)
genome_size = sum(chrom_sizes_hg38.values())
chromosome_size = []
p_values = []
variant_counts = df[0].value_counts()

for chrom in chrom_sizes_hg38:
    variant_count = variant_counts.get(chrom, 0)
    chromosome_size.append(chrom_sizes_hg38[chrom])
    expected_probability = chrom_sizes_hg38[chrom]/genome_size
    p_valuer = scipy.stats.binomtest(variant_count, total_variants, expected_probability).pvalue
    p_values.append(p_valuer)

f_d_r = scipy.stats.false_discovery_control(p_values)
print(f_d_r)

df_filtered = df[(df['clnsig'] == "Pathogenic") | (df['clnsig'] == "Benign")]
print(len(df_filtered))

count_per_chrom = df.groupby(df_filtered[0])['clnsig'].value_counts()
print(count_per_chrom)

genome_wide_pathogenic_ratio = len(df_filtered[df_filtered['clnsig'] == "Pathogenic"])/len(df_filtered)

p_values_chrom = []
chromosomes = []
x_axis = []

for chrom in chrom_sizes_hg38:
    chrom_data = count_per_chrom.get(chrom)
    if chrom.isdigit():
      chrom_key = int(chrom)
    else:
      chrom_key = chrom
    
    observed_pathogenic_count = count_per_chrom.get((chrom_key, "Pathogenic"), 0)
    observed_benign_count = count_per_chrom.get((chrom_key, "Benign"), 0)


    expected_pathogenic_count = (observed_pathogenic_count + observed_benign_count) * genome_wide_pathogenic_ratio
    expected_benign_count = (observed_benign_count + observed_pathogenic_count) * (1-genome_wide_pathogenic_ratio)


    if observed_pathogenic_count > 0 and observed_benign_count > 0:
      chi_square_test = scipy.stats.chi2_contingency([[observed_pathogenic_count, observed_benign_count],
                                                     [expected_pathogenic_count, expected_benign_count]])
      x_axis.append(observed_pathogenic_count/expected_pathogenic_count)

      p_values_chrom.append(chi_square_test.pvalue) 
      chromosomes.append(chrom_key)
    



f_d_r_chrom = scipy.stats.false_discovery_control(p_values_chrom)
print(chromosomes, p_values_chrom, f_d_r_chrom)

log2_ratio = [np.log2(val) for val in x_axis]
neg_log10_p = [-np.log10(val) for val in f_d_r_chrom]

plt.scatter(log2_ratio, neg_log10_p)
plt.xlabel("(Observed/Expected Pathogenic Ratio)")
plt.ylabel("FDR-adjusted p-values (-log10)")
plt.title("Pathogenic Variant Enrichment Across Chromosomes")
plt.show()




    

        









    
