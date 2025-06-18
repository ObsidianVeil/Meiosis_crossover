import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pyBigWig as pbw
import os
import fnmatch

barcode = os.environ["barcode"]

filename = f"{barcode}_coverage.bw" #set file name
bigwigfile = pbw.open(f"bigwigs/{filename}") #find the file

Chromlist = bigwigfile.chroms().keys()
Chromlist = list(Chromlist)

pattern_k = "k_chr*"
pattern_c = "c_chr*"
Chromlist_k = fnmatch.filter(Chromlist, pattern_k)
Chromlist_na = [chrom[2:] for chrom in Chromlist_k]

for chrom in Chromlist_na:
    k_chrom = f"k_{chrom}"
    c_chrom = f"c_{chrom}"


    k_chromlen=bigwigfile.chroms(k_chrom)
    print(barcode, chrom, "k chromosome length:", k_chromlen)

    k_maxvalue=bigwigfile.stats(k_chrom, 0, k_chromlen, type="max")
    k_cov = bigwigfile.values(k_chrom, 0, k_chromlen, numpy=False)

    c_chromlen=bigwigfile.chroms(c_chrom)
    print(barcode, chrom, "c chromosome length: ", c_chromlen)

    c_maxvalue=bigwigfile.stats(c_chrom, 0, c_chromlen, type="max")

    c_cov = bigwigfile.values(c_chrom, 0, c_chromlen, numpy=False)

    print(barcode, chrom, "k max value: ", k_maxvalue)
    print(barcode, chrom, "c max value: ", c_maxvalue)

    k_arr = np.array(k_cov)
    c_arr = np.array(c_cov)

    len_k = len(k_arr)
    len_c = len(c_arr)

    if len_k < len_c:
        k_arr = np.pad(k_arr, (0, len_c - len_k), "constant")
    elif len_c < len_k:
        c_arr = np.pad(c_arr, (0, len_k - len_c), 'constant')

    k_greater = k_arr > c_arr
    c_greater = c_arr > k_arr
    equal = k_arr == c_arr

    # Assign 1 where greater, 0.5 where equal, and 0 otherwise
    comparison_result_k = np.zeros_like(k_arr, dtype=float)
    comparison_result_k[k_greater] = 1.0
    comparison_result_k[c_greater] = 0.0  # redundant since array is initialized to 0
    comparison_result_k[equal] = 0.5

    # Suppose comparison is a NumPy array of 1s and 0s (1 = k > c)
    comparison_k = np.sign(k_arr - c_arr)

    # Extract (start, width) for regions where comparison == 1 (i.e. k > c)
    change_points = np.where(np.diff(comparison_k) != 0)[0] + 1
    starts = np.concatenate(([0], change_points))
    ends = np.concatenate((change_points, [len(comparison_k)]))
    widths = ends - starts
    labels = comparison_k[starts]

    comparison_result_k = [(s, w) for s, w, l in zip(starts, widths, labels) if l == 1]

    # Suppose comparison is a NumPy array of 1s and 0s (1 = k > c)
    comparison_c = np.sign(c_arr - k_arr)

    # Extract (start, width) for regions where comparison == 1 (i.e. k > c)
    change_points = np.where(np.diff(comparison_c) != 0)[0] + 1
    starts = np.concatenate(([0], change_points))
    ends = np.concatenate((change_points, [len(comparison_c)]))
    widths = ends - starts
    labels = comparison_c[starts]

    comparison_result_c = [(s, w) for s, w, l in zip(starts, widths, labels) if l == 1]

    plt.figure(figsize=(25, 5), layout="constrained")
    plt.broken_barh(comparison_result_k, (0, 1), facecolors="orange", label="Chromosome K")
    plt.broken_barh(comparison_result_c, (0, 1), facecolors="darkturquoise", label="Chromosome C")
    plt.legend(loc="upper right")
    plt.xlabel("Position")
    plt.ylabel("k or c")
    plt.title(f"{chrom}_coverage_plot")
    outputgraph = f"brokenbarplots/{barcode}/{chrom}.png"
    plt.savefig(outputgraph)