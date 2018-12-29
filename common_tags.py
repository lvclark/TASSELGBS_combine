#!/usr/bin/env python3
import sys

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
    markernames, alleleindex = zip(*sorted(zip(retained_markers[0], retained_markers[1])))
    # read in next SAM
    newtags = tagdigger_fun.readTags_TASSELSAM(sf)
    if newtags == None:
        sys.exit()
