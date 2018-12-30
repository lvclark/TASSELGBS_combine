#!/usr/bin/env python3
import sys
from bisect import bisect_left

# Directory for TagDigger module (tagdigger_fun.py should be in this folder)
td_dir = "../tagdigger/"
# File containing parameters for the pipeline
paramsFile = "params.txt"

# Import TagDigger functions
import importlib.util
spec = importlib.util.spec_from_file_location("tagdigger_fun", td_dir + "tagdigger_fun.py")
tagdigger_fun = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tagdigger_fun)

# Get SAM file names and other data from parameter file
popnames = []
samfiles = []
ttdfiles = []
min_with_reads = []
min_with_minor = []
with open(paramsFile, "rt") as mycon:
    for line in mycon:
        if line.startswith("Population name:"):
            popnames.append(line[16:].strip())
        if line.startswith("SAM file:"):
            samfiles.append(line[10:].strip())
        if line.startswith("TagTaxaDist file:"):
            ttdfiles.append(line[17:].strip())
        if line.startswith("Min. ind. with reads:"):
            min_with_reads.append(int(line[21:].strip()))
        if line.startswith("Min. ind. with minor allele:"):
            min_with_minor.append(int(line[28:].strip()))
        if line.startswith("Tag database output file:"):
            db_outfile = line[25:].strip()

npops = len(popnames)
if npops < 2:
    sys.exit("Fewer than two populations listed.")
if len(samfiles) != npops or len(ttdfiles) != npops or \
  len(min_with_reads) != npops or len(min_with_minor) != npops:
    sys.exit("Params file not formatted correctly.")

print("Finding tag locations common to {} files:".format(npops))
for s in samfiles:
    print(s)

# Read first SAM file
retained_tags = tagdigger_fun.readTags_TASSELSAM(samfiles[0])
if retained_tags == None:
    sys.exit()

# Read remaining SAM files and pare down to shared cut sites
for sf in samfiles[1:]:
    # organize markers from last round
    retained_markers = tagdigger_fun.extractMarkers(retained_tags[0])
    # sort the retained markers; not necessary but will leave in for now
    markernames_r, alleleindex_r = zip(*sorted(zip(retained_markers[0], retained_markers[1])))
    # read in next SAM
    newtags = tagdigger_fun.readTags_TASSELSAM(sf)
    if newtags == None:
        sys.exit()
    new_markers = tagdigger_fun.extractMarkers(newtags[0])
    markernames_n, alleleindex_n = zip(*sorted(zip(new_markers[0], new_markers[1])))
    # set up lists for new retained tags
    retained_tags_n = [[], []]

    # match retained and new markers
    for mi in range(len(markernames_r)):
        m = markernames_r[mi]
        i = bisect_left(markernames_n, m)
        if markernames_n[i] == m:
            tags_old = [retained_tags[1][j] for j in alleleindex_r[mi][1]]
            tags_new = [newtags[1][j] for j in alleleindex_n[i][1]]
            tags_all = list(set(tags_old + tags_new))
            tags_all.sort()
            ct = tagdigger_fun.compareTags(tags_all)
            allele_names = [m + "_" + "".join([ct1[1][j] for ct1 in ct]) \
                            for j in range(len(tags_all))]
            # eliminate cases of one tag being shorter version of another
            if len(set(allele_names)) < len(allele_names):
                continue
            retained_tags_n[0].extend(allele_names)
            retained_tags_n[1].extend(tags_all)
    retained_tags = retained_tags_n

# filter out monomorphic markers
retained_tags = tagdigger_fun.remove_monomorphic_loci(retained_tags[0], retained_tags[1])

# prepare to export common markers to a Tag Manager database
mtags = tagdigger_fun.mergedTagList(retained_tags)
# extract chromosomes and positions, and format for database
splittemp = [markername.split("-") for markername in mtags[0]]
chrom = [s[0] for s in splittemp]
pos = [int(s[1]) for s in splittemp]
chrompos = dict(zip(mtags[0], list(zip(chrom, pos))))
# write file
tagdigger_fun.writeMarkerDatabase(db_outfile, mtags[0], mtags[1], \
                                  [[["Chromosome", "Posiiton"], chrompos]])
