import os
import fnmatch
import csv
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

# --- Helpers to work only with numeric region pairs (skip sign_changes) ---
NUMERIC_CATS = ("k>c", "k<c", "k==c", "k==c==0")

def iter_numeric_regions(d, cats=NUMERIC_CATS):
    """Yield (start, width) as ints from the numeric region lists only."""
    for cat in cats:
        for s, w in d.get(cat, []):
            try:
                yield int(s), int(w)
            except (TypeError, ValueError):
                # Skip anything not cleanly coercible to ints
                continue

def max_extent(d, fallback=1):
    """Max of (start+width) over numeric regions, with a fallback."""
    m = fallback
    for s, w in iter_numeric_regions(d):
        m = max(m, s + w)
    return m
################################################################################

def log1p_formatter(y, pos):
    return f"{np.expm1(y):.0f}"


def extract_sign_changes(region_labels, region_starts):
    """Return list of (position, change_str) for +/- transitions only.
    region_labels: array of labels per region (1, -1, 0, or 2)
    region_starts: array of start indices for each region
    We record a change at the boundary i where labels[i-1] -> labels[i]
    Only keep +1<->-1 changes.
    """
    events = []
    for i in range(1, len(region_labels)):
        prev_l = region_labels[i - 1]
        curr_l = region_labels[i]
        if prev_l in (1, -1) and curr_l in (1, -1) and prev_l != curr_l:
            pos = int(region_starts[i])  # boundary index (start of new region)
            if prev_l == 1 and curr_l == -1:
                change = "k>c_to_k<c"
            else:  # prev -1 to +1
                change = "k<c_to_k>c"
            events.append((pos, change))
    return events


# Process each group of 4 bigwig files
for sublist in splitlist:
    foldername = ""

    for item in sublist:
        barcode = item[:9]
        foldername += barcode

        bigwigfile = pbw.open(os.path.join(inputdir, item))
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

            # Regionization
            change_points = np.where(np.diff(comparison_k) != 0)[0] + 1
            starts = np.concatenate(([0], change_points))
            ends = np.concatenate((change_points, [len(comparison_k)]))
            widths = ends - starts
            labels = comparison_k[starts]

            # Extract just the +/- sign change events
            sign_change_events = extract_sign_changes(labels, starts)

            key = f"{barcode}_{chrom}"
            comparison_results[key] = {
                "k>c": [(int(s), int(w)) for s, w, l in zip(starts, widths, labels) if l == 1],
                "k<c": [(int(s), int(w)) for s, w, l in zip(starts, widths, labels) if l == -1],
                "k==c": [(int(s), int(w)) for s, w, l in zip(starts, widths, labels) if l == 0],
                "k==c==0": [(int(s), int(w)) for s, w, l in zip(starts, widths, labels) if l == 2],
                "k_cov": k_arr,
                "c_cov": c_arr,
                "sign_changes": sign_change_events,  # list of (position, change_str)
            }

    # === Compute nearest inverse sign-change across the other three barcodes and write CSV ===
    outfolder = f"{outputdir}/{foldername}"
    os.makedirs(outfolder, exist_ok=True)

    # collect chromosomes present for this group
    chrom_set = set("_".join(k.split("_")[1:]) for k in comparison_results.keys() if k.startswith(sublist[0][:9]))

    csv_path = os.path.join(outfolder, "transition_matches.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "group_id", "chromosome", "barcode", "position", "change",
            "nearest_barcode", "nearest_position", "nearest_change", "distance"
        ])

        for chrom in sorted(chrom_set):
            # build per-barcode lists for this chromosome
            per_bc = {}
            for item in sublist:
                bc = item[:9]
                key = f"{bc}_{chrom}"
                if key in comparison_results:
                    per_bc[bc] = comparison_results[key]["sign_changes"]
                else:
                    per_bc[bc] = []

            # for each barcode's events, find nearest inverse on others
            for bc, events in per_bc.items():
                # flatten other events with barcode tag
                others = []
                for obc, oevents in per_bc.items():
                    if obc == bc:
                        continue
                    for pos, chg in oevents:
                        others.append((obc, pos, chg))

                for pos, chg in events:
                    target = "k<c_to_k>c" if chg == "k>c_to_k<c" else "k>c_to_k<c"
                    best = None  # (distance, obc, opos, ochg)
                    for obc, opos, ochg in others:
                        if ochg != target:
                            continue
                        dist = abs(opos - pos)
                        if best is None or dist < best[0]:
                            best = (dist, obc, opos, ochg)
                    if best is None:
                        writer.writerow([foldername, chrom, bc, pos, chg, "", "", "", ""])
                    else:
                        dist, obc, opos, ochg = best
                        writer.writerow([foldername, chrom, bc, pos, chg, obc, opos, ochg, dist])

    # === Plotting (unchanged) ===
    max_position = 0
    for item in sublist:
        barcode = item[:9]
        for key in comparison_results:
            if key.startswith(barcode) and isinstance(comparison_results[key], dict):
                for s, w in iter_numeric_regions(comparison_results[key]):
                    max_position = max(max_position, s + w)

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

            ax.broken_barh(data["k>c"], (0, max_cov_log), facecolors="yellow", label="k > c")
            ax.broken_barh(data["k<c"], (0, max_cov_log), facecolors="blue", label="k < c")
            ax.broken_barh(data["k==c"], (0, max_cov_log), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, max_cov_log), facecolors="white", label="k = c = 0")

            ax.plot(x_vals, k_cov_log, color="deepskyblue", linewidth=1, label="K coverage")
            ax.plot(x_vals, c_cov_log, color="gold", linewidth=1, label="C coverage")

            ax.yaxis.set_major_formatter(FuncFormatter(log1p_formatter))
            ax.set_yticks([max_cov_log])
            ax.set_ylabel(f"{barcode}\nlog(coverage + 1)")
            ax.set_xlim(0, max_extent(data, fallback=1))
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

    def chrom_sort_key(chrom):
        chrom = chrom.lower().replace("chr", "")
        return int(chrom) if chrom.isdigit() else float('inf')

    sorted_chroms = sorted(chrom_set, key=chrom_sort_key)

    max_xlim = 1
    for data in comparison_results.values():
        if isinstance(data, dict):
            max_xlim = max(max_xlim, max_extent(data, fallback=0))

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

            ax.broken_barh(data["k>c"], (0, max_cov_log), facecolors="yellow", label="k > c")
            ax.broken_barh(data["k<c"], (0, max_cov_log), facecolors="blue", label="k < c")
            ax.broken_barh(data["k==c"], (0, max_cov_log), facecolors="gray", label="k = c")
            ax.broken_barh(data["k==c==0"], (0, max_cov_log), facecolors="white", label="k = c = 0")

            ax.plot(x_vals, k_cov_log, color="deepskyblue", linewidth=1, label="K coverage")
            ax.plot(x_vals, c_cov_log, color="gold", linewidth=1, label="C coverage")

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
