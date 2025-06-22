import numpy as np
import matplotlib.pyplot as plt
import pyBigWig as pbw
import os
import fnmatch

# List and sort all files in the 'bigwigs' directory
barcodelist = os.listdir("bigwigs")
barcodelist = fnmatch.filter(barcodelist, "*.bw")
barcodelist.sort()

# Split the list of files into groups of 4
splitlist = [barcodelist[i:i+4] for i in range(0, len(barcodelist), 4)]

# Dictionary to hold comparison results per barcode and chromosome
comparison_results = {}

# Process each group of 4 bigwig files
for sublist in splitlist:
    foldername = ""

    for item in sublist:
        barcode = item[:9]  # Assume first 9 characters identify the barcode
        foldername += barcode

        # Open the bigwig file
        bigwigfile = pbw.open(f"bigwigs/{item}")

        # Get list of chromosome names
        chromlist = bigwigfile.chroms().keys()

        # Filter chromosome names starting with 'k_chr'
        chromlist_k = fnmatch.filter(chromlist, "k_chr*")

        # Remove the 'k_' prefix to get raw chromosome names
        chromlist_na = [chrom[2:] for chrom in chromlist_k]

        for chrom in chromlist_na:
            k_chrom = f"k_{chrom}"
            c_chrom = f"c_{chrom}"

            # Get chromosome lengths
            k_chromlen = bigwigfile.chroms(k_chrom)
            c_chromlen = bigwigfile.chroms(c_chrom)

            # Extract signal values for each chromosome
            k_cov = bigwigfile.values(k_chrom, 0, k_chromlen, numpy=False)
            c_cov = bigwigfile.values(c_chrom, 0, c_chromlen, numpy=False)

            # Convert to NumPy arrays
            k_arr = np.array(k_cov)
            c_arr = np.array(c_cov)

            # Pad arrays to equal length
            len_k = len(k_arr)
            len_c = len(c_arr)

            if len_k < len_c:
                k_arr = np.pad(k_arr, (0, len_c - len_k), "constant")
            elif len_c < len_k:
                c_arr = np.pad(c_arr, (0, len_k - len_c), 'constant')

            # Compare signal: 1 if k > c, -1 if k < c, 0 if equal
            comparison_k = np.sign(k_arr - c_arr)
            # Mark regions where both k and c are zero with a special code (2)
            comparison_k[(k_arr == c_arr) & (k_arr == 0)] = 2

            # Identify change points where comparison state changes
            change_points = np.where(np.diff(comparison_k) != 0)[0] + 1
            starts = np.concatenate(([0], change_points))
            ends = np.concatenate((change_points, [len(comparison_k)]))
            widths = ends - starts
            labels = comparison_k[starts]

            # Store regions by comparison category
            key = f"{barcode}_{chrom}"
            comparison_results[key] = {
                "k>c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 1],
                "k<c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == -1],
                "k==c": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 0],
                "k==c==0": [(s, w) for s, w, l in zip(starts, widths, labels) if l == 2],  # Special marker
            }

    # Determine the maximum position across all comparisons for x-axis limit
    max_position = 0
    for item in sublist:
        barcode = item[:9]
        for key in comparison_results:
            if key.startswith(barcode):
                for region_list in comparison_results[key].values():
                    for start, width in region_list:
                        max_position = max(max_position, start + width)

    # Create output folder named by concatenated barcodes
    foldername = "".join([item[:9] for item in sublist])
    outfolder = f"mergedbarplots/{foldername}"
    os.makedirs(outfolder, exist_ok=True)

    # Get set of chromosomes (without barcode prefix) from the first file in sublist
    chrom_set = set("_".join(k.split("_")[1:]) for k in comparison_results.keys() if k.startswith(sublist[0][:9]))

    # Create bar plot per chromosome
    for chrom in chrom_set:
        fig, axs = plt.subplots(nrows=len(sublist), figsize=(25, 2.5 * len(sublist)), sharex=True)

        # Ensure axs is iterable
        if len(sublist) == 1:
            axs = [axs]

        for i, item in enumerate(sublist):
            barcode = item[:9]
            ax = axs[i]
            key = f"{barcode}_{chrom}"

            if key not in comparison_results:
                continue

            data = comparison_results[key]

            # Plot different comparison regions
            ax.broken_barh(data["k>c"], (0, 15), facecolors="orange", label="k > c")
            ax.broken_barh(data["k<c"], (0, 15), facecolors="darkturquoise", label="k < c")
            ax.broken_barh(data["k==c"], (0, 15), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, 15), facecolors="white", label="k = c = 0")

            ax.set_yticks([15])
            ax.set_yticklabels([])
            ax.set_ylabel(barcode)
            ax.set_xlim(0, max([s + w for cat in data.values() for s, w in cat] + [1]))

            if i == 0:
                ax.set_title(f"Chromosome {chrom} Comparison")

            if i == len(sublist) - 1:
                ax.set_xlabel("Genomic Position")

        # Only show unique legend labels once
        handles, labels = axs[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        fig.legend(by_label.values(), by_label.keys(), loc="upper right")

        # Save the figure
        outputgraph = f"{outfolder}/{chrom}.png"
        plt.tight_layout()
        plt.savefig(outputgraph)
        plt.close()


    # Sort chromosomes naturally: chr1, chr2, ..., chrX
    def chrom_sort_key(chrom):
        chrom = chrom.lower().replace("chr", "")
        return int(chrom) if chrom.isdigit() else float('inf')


    sorted_chroms = sorted(chrom_set, key=chrom_sort_key)

    # Compute global max x-position across all data
    max_xlim = 1
    for data in comparison_results.values():
        xlim = max([s + w for cat in data.values() for s, w in cat], default=0)
        max_xlim = max(max_xlim, xlim)

    # Prepare subplots
    total_plots = len(sorted_chroms) * len(sublist)
    fig, axs = plt.subplots(
        nrows=total_plots,
        figsize=(25, 2.5 * total_plots),
        sharex=True
    )

    if total_plots == 1:
        axs = [axs]

    plot_idx = 0

    for chrom in sorted_chroms:
        allgraphname = ""
        for item in sublist:
            barcode = item[:9]
            allgraphname = allgraphname + barcode
            key = f"{barcode}_{chrom}"

            if key not in comparison_results:
                continue

            data = comparison_results[key]
            ax = axs[plot_idx]

            ax.broken_barh(data["k>c"], (0, 15), facecolors="orange", label="k > c")
            ax.broken_barh(data["k<c"], (0, 15), facecolors="darkturquoise", label="k < c")
            ax.broken_barh(data["k==c"], (0, 15), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, 15), facecolors="white", label="k = c = 0")

            ax.set_yticks([15])
            ax.set_yticklabels([])
            ax.set_ylabel(f"{barcode}", rotation=0, ha='right')
            ax.set_xlim(0, max_xlim)

            if item == sublist[0]:
                ax.set_title(f"Chromosome {chrom}")

            if item == sublist[-1]:
                ax.set_xlabel("Genomic Position")

            plot_idx += 1

    # Add unified legend
    handles, labels = axs[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc="upper right")

    plt.tight_layout(rect=[0, 0, 0.98, 1])
    outputgraph = f"{outfolder}/{allgraphname}_all_chromosomes.png"
    plt.savefig(outputgraph)
    plt.close()
