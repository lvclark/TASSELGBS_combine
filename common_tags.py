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

# Get SAM file names from parameter file
samfiles = []
with open(paramsFile, "rt") as mycon:
    for line in mycon:
        if line.startswith("SAM file:"):
            samfiles.append(line[10:].strip())

if len(samfiles) < 2:
    sys.exit("Fewer than two SAM files listed.")

print("Finding tag locations common to {} files:".format(len(samfiles)))
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
            retained_tags_n[0].extend(allele_names)
            retained_tags_n[1].extend(tags_all)
    retained_tags = retained_tags_n

# filter out monomorphic markers
retained_tags = tagdigger_fun.remove_monomorphic_loci(retained_tags[0], retained_tags[1])

# prepare to export common markers to a Tag Manager database
mtags = tagdigger_fun.mergedTagList(retained_tags)
# extract chromosomes and positions, and format for database
# write file
tagdigger_fun.writeMarkerDatabase()
