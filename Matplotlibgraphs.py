import os
import fnmatch
import numpy as np
import pyBigWig as pbw
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.ticker import FuncFormatter

# List and sort all files in the 'bigwigs' directory
inputdir = str(os.environ.get("bigwigfolder"))
barcodelist = os.listdir(inputdir)
barcodelist = fnmatch.filter(barcodelist, "*.bw")
barcodelist.sort()

outputdir = str(os.environ.get("graphoutput"))

# Split the list of files into groups of 4
splitlist = [barcodelist[i:i + 4] for i in range(0, len(barcodelist), 4)]

# Dictionary to hold comparison results per barcode and chromosome
comparison_results = {}


def log1p_formatter(y, pos):
    return f"{np.expm1(y):.0f}"


# Process each group of 4 bigwig files
for sublist in splitlist:
    foldername = ""

    for item in sublist:
        barcode = item[:9]
        foldername += barcode

        bigwigfile = pbw.open(f"bigwigs/{item}")
        chromlist = bigwigfile.chroms().keys()
        chromlist_k = fnmatch.filter(chromlist, "k_chr*")
        chromlist_na = [chrom[2:] for chrom in chromlist_k]

        for chrom in chromlist_na:
            k_chrom = f"k_{chrom}"
            c_chrom = f"c_{chrom}"

            k_chromlen = bigwigfile.chroms(k_chrom)
            c_chromlen = bigwigfile.chroms(c_chrom)

            k_cov = bigwigfile.values(k_chrom, 0, k_chromlen, numpy=False)
            c_cov = bigwigfile.values(c_chrom, 0, c_chromlen, numpy=False)

            k_arr = np.array(k_cov)
            c_arr = np.array(c_cov)

            len_k = len(k_arr)
            len_c = len(c_arr)

            if len_k < len_c:
                k_arr = np.pad(k_arr, (0, len_c - len_k), "constant")
            elif len_c < len_k:
                c_arr = np.pad(c_arr, (0, len_k - len_c), 'constant')

            comparison_k = np.sign(k_arr - c_arr)
            comparison_k[(k_arr == c_arr) & (k_arr == 0)] = 2

            change_points = np.where(np.diff(comparison_k) != 0)[0] + 1
            starts = np.concatenate(([0], change_points))
            ends = np.concatenate((change_points, [len(comparison_k)]))
            widths = ends - starts
            labels = comparison_k[starts]

            key = f"{barcode}_{chrom}"
            comparison_results[key] = {
                "k>c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 1],
                "k<c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == -1],
                "k==c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 0],
                "k==c==0": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 2],
                "k_cov": k_arr,
                "c_cov": c_arr
            }

    max_position = 0
    for item in sublist:
        barcode = item[:9]
        for key in comparison_results:
            if key.startswith(barcode):
                for region_list in comparison_results[key].values():
                    if isinstance(region_list, list):
                        for item in region_list:
                            if isinstance(item, (tuple, list)) and len(item) == 2:
                                start, width = item
                                max_position = max(max_position, start + width)

    outfolder = f"{outputdir}/{foldername}"
    os.makedirs(outfolder, exist_ok=True)

    chrom_set = set("_".join(k.split("_")[1:]) for k in comparison_results.keys() if k.startswith(sublist[0][:9]))

    for chrom in chrom_set:
        fig, axs = plt.subplots(nrows=len(sublist), figsize=(25, 2.5 * len(sublist)), sharex=True)
        if len(sublist) == 1:
            axs = [axs]

        for i, item in enumerate(sublist):
            barcode = item[:9]
            ax = axs[i]
            key = f"{barcode}_{chrom}"

            if key not in comparison_results:
                continue
            data = comparison_results[key]

            x_vals = np.arange(len(data["k_cov"]))
            k_cov_log = np.log1p(data["k_cov"])
            c_cov_log = np.log1p(data["c_cov"])
            max_cov_log = max(k_cov_log.max(), c_cov_log.max())

            data = comparison_results[key]
            ax.broken_barh(data["k>c"], (0, max_cov_log), facecolors="orange", label="k > c")
            ax.broken_barh(data["k<c"], (0, max_cov_log), facecolors="darkturquoise", label="k < c")
            ax.broken_barh(data["k==c"], (0, max_cov_log), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, max_cov_log), facecolors="white", label="k = c = 0")

            ax.plot(x_vals, k_cov_log, color="purple", linewidth=1, label="K coverage")
            ax.plot(x_vals, c_cov_log, color="green", linewidth=1, label="C coverage")

            ax.yaxis.set_major_formatter(FuncFormatter(log1p_formatter))
            ax.set_yticks([max_cov_log])
            ax.set_ylabel(f"{barcode}\nlog(coverage + 1)")
            ax.set_xlim(0, max([s + w for cat in data.values() if isinstance(cat, list) for s, w in cat] + [1]))
            ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

            if i == 0:
                ax.set_title(f"Chromosome {chrom} Comparison")

            if i == len(sublist) - 1:
                ax.set_xlabel("Genomic Position")

        handles, labels = axs[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        fig.legend(by_label.values(), by_label.keys(), loc="upper right")

        plt.tight_layout()
        outputgraph = f"{outfolder}/{chrom}.png"
        plt.savefig(outputgraph)
        plt.close()
        # plt.show()


    def chrom_sort_key(chrom):
        chrom = chrom.lower().replace("chr", "")
        return int(chrom) if chrom.isdigit() else float('inf')


    sorted_chroms = sorted(chrom_set, key=chrom_sort_key)

    max_xlim = 1
    for data in comparison_results.values():
        if isinstance(data, dict):
            xlim = max([s + w for cat in data.values() if isinstance(cat, list) for s, w in cat], default=0)
            max_xlim = max(max_xlim, xlim)

    total_plots = len(sorted_chroms) * len(sublist)
    fig, axs = plt.subplots(nrows=total_plots, figsize=(25, 2.5 * total_plots), sharex=True)
    if total_plots == 1:
        axs = [axs]

    plot_idx = 0
    for chrom in sorted_chroms:
        allgraphname = ""
        for item in sublist:
            barcode = item[:9]
            allgraphname += barcode
            key = f"{barcode}_{chrom}"

            if key not in comparison_results:
                continue
            data = comparison_results[key]

            x_vals = np.arange(len(data["k_cov"]))
            k_cov_log = np.log1p(data["k_cov"])
            c_cov_log = np.log1p(data["c_cov"])
            max_cov_log = max(k_cov_log.max(), c_cov_log.max())

            ax = axs[plot_idx]

            ax.broken_barh(data["k>c"], (0, max_cov_log), facecolors="orange", label="k > c")
            ax.broken_barh(data["k<c"], (0, max_cov_log), facecolors="darkturquoise", label="k < c")
            ax.broken_barh(data["k==c"], (0, max_cov_log), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, max_cov_log), facecolors="white", label="k = c = 0")

            x_vals = np.arange(len(data["k_cov"]))
            ax.plot(x_vals, k_cov_log, color="purple", linewidth=1, label="K coverage")
            ax.plot(x_vals, c_cov_log, color="green", linewidth=1, label="C coverage")

            ax.set_yticks([max_cov_log])
            ax.set_yticklabels([])
            ax.set_ylabel(f"{barcode}", rotation=0, ha='right')
            ax.set_xlim(0, max_xlim)
            ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

            if item == sublist[0]:
                ax.set_title(f"Chromosome {chrom}")

            if item == sublist[-1]:
                ax.set_xlabel("Genomic Position")

            plot_idx += 1

    handles, labels = axs[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="upper right")

    plt.tight_layout(rect=[0, 0, 0.98, 1])
    outputgraph = f"{outfolder}/{allgraphname}_all_chromosomes.png"
    plt.savefig(outputgraph)
    plt.close()
    # plt.show()
