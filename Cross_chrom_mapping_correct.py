import pysam
from collections import defaultdict
import os

infilepath = str(os.environ.get("samtoolssort_output"))
outfilepath = str(os.environ.get("pyoutput1"))

if not infilepath or not os.path.exists(infilepath):
    raise ValueError(f"Invalid input file path: {infilepath}")

print("Reading from:", infilepath)
print("Writing to:", outfilepath)

bamfile = pysam.AlignmentFile(infilepath, "rb")
header = bamfile.header

# Store reads by query name
reads_by_name = defaultdict(list)

# First pass: collect all reads by query name
for read in bamfile.fetch(until_eof=True):
    reads_by_name[read.query_name].append(read)

# Output reads where primary and supplementary align to different chromosomes
outreads = []

for qname, read_group in reads_by_name.items():
    primary_read = None
    supp_reads = []

    for read in read_group:
        if read.is_secondary:
            continue
        if read.is_supplementary:
            supp_reads.append(read)
        else:
            primary_read = read
            
    #if not primary_read and supp_reads:
    #    print(f"[{qname}] has supplementary but no primary")
    #elif primary_read and not supp_reads:
    #    print(f"[{qname}] has primary but no supplementary")
    #elif not primary_read and not supp_reads:
    #    print(f"[{qname}] has neither primary nor supplementary")
    #else:
    #    print(f"[{qname}] has primary AND supplementary")

    if primary_read and supp_reads:
        for supp_read in supp_reads:
            if primary_read.reference_name != supp_read.reference_name:
                outreads.append(primary_read)
                outreads.append(supp_read)

print("Python writing final BAM output to:", outfilepath)
with pysam.AlignmentFile(outfilepath, "wb", header=header) as outbam:
    for read in outreads:
        outbam.write(read)

print(f"Written {len(outreads)} reads to {outfilepath}")
