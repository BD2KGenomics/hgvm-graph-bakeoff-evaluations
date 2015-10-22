#!/usr/bin/env python2.7
"""
redoStats.py: rerun stats on a GAM file

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit

stats = {
    "total_reads": 0,
    "total_mapped": 0,
    "total_multimapped": 0,
    "primary_scores": collections.Counter(),
    "primary_mismatches": collections.Counter(),
    "secondary_scores": collections.Counter(),
    "secondary_mismatches": collections.Counter(),
    "run_time": None
}

last_alignment = None

for line in sys.stdin:
    # Parse the alignment JSON
    alignment = json.loads(line)
    
    if alignment.has_key("score"):
        # This alignment is aligned.
        # Grab its score
        score = alignment["score"]
    
        # Calculate the mismatches
        length = len(alignment["sequence"])
        matches = 0
        for mapping in alignment.get("path", {}).get("mapping", []):
            for edit in mapping.get("edit", []):
                if (not edit.has_key("sequence") and 
                    edit.get("to_length", None) == edit.get(
                    "from_length", None)):
                    
                    # We found a perfect match edit. Grab its length
                    matches += edit["from_length"]
                    
        # Calculate mismatches as what's left
        mismatches = length - matches
                
    
        if alignment.get("is_secondary", False):
            if (last_alignment is None or 
                last_alignment.get("name") != alignment.get("name") or 
                last_alignment.get("is_secondary", False)):
                
                # This is a secondary alignment without a corresponding primary
                # alignment
                #sys.stderr.write("{} secondary alignment comes after "
                #    "alignment of {} instead of corresponding primary "
                #    "alignment\n".format(alignment.get("name"), 
                #    last_alignment.get("name") if last_alignment is not None 
                #    else "nothing"))
                    
                # Correct it to be primary
                alignment["is_secondary"] = False
                
        
        if alignment.get("is_secondary", False):
            # It's a multimapping. We can have max 1 per read, so it's a
            # multimapped read.
            
            # Log its stats as multimapped
            stats["total_multimapped"] += 1
            stats["secondary_scores"][score] += 1
            stats["secondary_mismatches"][mismatches] += 1
        else:
            # Log its stats as primary. We'll get exactly one of these per
            # read with any mappings.
            stats["total_mapped"] += 1
            stats["primary_scores"][score] += 1
            stats["primary_mismatches"][mismatches] += 1
            
            # We won't see an unaligned primary alignment for this read, so
            # count the read
            stats["total_reads"] += 1
            
        last_alignment = alignment
    
    elif not alignment.get("is_secondary", False):
        # We have an unmapped primary "alignment"
        
        # Count the read by its primary alignment
        stats["total_reads"] += 1
            
# Save the stats as JSON
json.dump(stats, sys.stdout)
