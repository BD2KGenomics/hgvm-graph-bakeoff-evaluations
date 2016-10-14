#!/usr/bin/env python2.7
"""
collateStatistics.py: turn individual .stats files into input files for plotting.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import tempfile
import copy

import numpy
import tsv

from toil.job import Job

from toillib import *

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("in_store", type=IOStore.absolute,
        help="input IOStore to find stats files in (under /stats)")
    parser.add_argument("out_store", type=IOStore.absolute,
        help="output IOStore to put collated plotting files in (under /plots)")
    parser.add_argument("--blacklist", action="append", default=[],
        help="ignore the specified regions, graphs, or region:graph pairs")
    parser.add_argument("--overwrite", action="store_true",
        help="replace cached per-sample statistics with recalculated ones")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def collate_all(job, options):
    """
    Collate all the stats files
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Holds a dict from region name to absolute stats by graph, sample, and stat
    # name.
    stats_by_region = {}
    
    for region in in_store.list_input_directory("stats"):
        # Collate everything in the region, and capture the returned stats dict
        # promise
        stats_by_region[region] = job.addChildJobFn(collate_region, options,
            region, cores=1, memory="1G", disk="10G").rv()
    
    # We also need to merge some overall stats across regions
    job.addFollowOnJobFn(merge_overall, options, stats_by_region,
        cores=1, memory="1G", disk="10G")
    
    
# TODO: use this!        
def concatenate_results(job, options, inputs, output_key):
    """
    Concatenate the given file IDs, and write the result to the given output
    store key.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Make the file to concatenate into
    (handle, path) = tempfile.mkstemp(dir=job.fileStore.getLocalTempDir())
    os.close(handle)     
    
    with open(path, "w") as destination:
        for input_id in inputs:
            with job.file_store.readGlobalFileStream(input_id) as input_stream:
                # Copy over every input
                shutil.copyfileobj(input_stream, destination)
    
    # Save the result
    out_store.write_output_file(path, output_key)

    
    
            
def collate_region(job, options, region):
    """
    Collate all the stats files in a region. Returns a dict from graph and
    sample and stat name to stat value, which may be cached.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # This holds a dict from graph name, then sample name, then stat name to
    # actual stst value.
    stats_cache = collections.defaultdict(lambda: collections.defaultdict(dict))
    
    # This is the cache file for this region, in
    # <graph>\t<sample>\t<stat>\t<value> format
    cache_tsv_key = "plots/cache/{}.tsv".format(region)
    
    # What name will it have locally for us?
    local_filename = os.path.join(job.fileStore.getLocalTempDir(), "temp.tsv")
    
    if out_store.exists(cache_tsv_key) and not options.overwrite:
        # Just read in from that TSV
        
        RealTimeLogger.get().info("Loading cached region {}".format(region))
        
        # Grab the cached results
        out_store.read_input_file(cache_tsv_key, local_filename)
        
        # Read all the pop, value pairs from the TSV
        reader = tsv.TsvReader(open(local_filename))
        
        for graph, sample, stat, value in reader:
            # Read in and place all the values
            stats_cache[graph][sample][stat] = float(value)
        
    else:
    
        # We need to populate the stats cache ourselves.
        for graph in in_store.list_input_directory("stats/{}".format(region)):
            # For each graph
        
            if ("{}:{}".format(region, graph) in options.blacklist or 
                region in options.blacklist or
                graph in options.blacklist):
                
                # We don't want to process this region/graph pair.
                RealTimeLogger.get().info("Skipping {} graph {}".format(region,
                    graph))
                continue
            
            RealTimeLogger.get().info("Processing {} graph {}".format(region,
                graph))
            for result in in_store.list_input_directory("stats/{}/{}".format(
                region, graph)):
            
                # For every sample
                
                # Pull sample name from filename
                sample_name = re.match("(.*)\.json$", result).group(1)
            
                # Grab file
                json_filename = os.path.join(job.fileStore.getLocalTempDir(),
                    "temp.json")
                in_store.read_input_file("stats/{}/{}/{}".format(region, graph,
                    result), json_filename)
                # Read the JSON
                stats = json.load(open(json_filename))
                
                # Compute the answers for this sample
                
                # How many reads are mapped well enough?
                total_mapped_well = sum((
                    stats["primary_matches_per_column"].get(x)
                    for x in stats["primary_matches_per_column"].iterkeys()
                    if float(x) >= 0.95))
                    
                # How many reads are multimapped well enough?
                total_multimapped_well = sum((
                    stats["secondary_matches_per_column"].get(x)
                    for x in stats["secondary_matches_per_column"].iterkeys()
                    if float(x) >= 0.95))
                # How many reads multimapped at all?
                total_multimapped_at_all = sum((
                    stats["secondary_matches_per_column"].get(x)
                    for x in stats["secondary_matches_per_column"].iterkeys()))
                    
                # How many reads are perfect?
                total_perfect = sum((
                    stats["primary_matches_per_column"].get(x)
                    for x in stats["primary_matches_per_column"].iterkeys()
                    if float(x) == 1.0))
                    
                # How many reads have no indels?
                total_no_indels = sum((stats["primary_indels"].get(str(x), 0)
                    for x in xrange(1)))
                
                # How many have any?    
                total_with_indels = sum((count for indels, count in 
                    stats["primary_indels"].iteritems() if count != "0"))
                    
                # How many have one? We'll guess these aren't terrible mappings.
                total_one_indel = stats["primary_indels"].get("1", 0)
                    
                # How many total substitution bases are there?
                substitution_bases = sum((count * float(weight)
                    for weight, count in 
                    stats["primary_substitutions"].iteritems()))
                    
                # How many total mismatches (substitutions + indels) are there?
                mismatch_bases = sum((count * float(weight)
                    for weight, count in 
                    stats["primary_mismatches"].iteritems()))
                    
                # How many total indels are there?
                indel_instances = sum((count * float(weight)
                    for weight, count in 
                    stats["primary_indels"].iteritems()))
                    
                # How many reads are mapped with no substitutions?
                total_no_substitutions = sum(
                    (stats["primary_substitutions"].get(str(x), 0)
                    for x in xrange(1)))
                    
                # What's the total matches per column for primary mappings
                total_matches_per_column = sum((count * float(weight)
                    for weight, count in 
                    stats["primary_matches_per_column"].iteritems()))
                    
                # What's the total score?
                total_score = sum((count * float(weight)
                    for weight, count in 
                    stats["primary_scores"].iteritems()))
                
                # How many reads are there overall for this sample?
                total_reads = stats["total_reads"]
                
                # If no reads got aligned, the sample is broken and we want to
                # skip it (and make the user fix it).
                if total_reads == 0:
                    RealTimeLogger.get().warning(
                        "No reads available for {} {} {}".format(
                        region, graph, sample_name))
                    continue
                
                # How many reads are mapped at all for this sample (not just
                # good enough)?
                total_mapped_at_all = stats["total_mapped"]
                
                # What was the runtime?
                runtime = stats.get("run_time", None)
                if runtime is None:
                    # We need NaN floats if there's no runtime
                    runtime = float("nan")
                else:
                    # Convert to time per read aligned
                    runtime /= total_reads
                    
                # Compute the stats we actually care about and save them
                
                # Get the dict we want to put the computed stats in
                sample_stats = stats_cache[graph][sample_name]
                
                # How many reads are mapped well, single-mapped well, and
                # existing at all? We need these for computing an overall
                # perfect vs unique plot.
                sample_stats["perfect"] = total_perfect
                sample_stats["single_mapped_well"] = (total_mapped_well -
                    total_multimapped_well)
                sample_stats["total_reads"] = total_reads
                
                # What portion have one good mapping?
                sample_stats["portion_single_mapped_well"] = \
                    ((total_mapped_well - total_multimapped_well) / 
                    float(total_reads))
                # What portion only have one mapping at all?
                sample_stats["portion_single_mapped_at_all"] = \
                    ((total_mapped_at_all - total_multimapped_at_all) / 
                    float(total_reads))
                # What portion are mapped well?
                sample_stats["portion_mapped_well"] = (total_mapped_well /
                    float(total_reads))
                # What portion are perfect?
                sample_stats["portion_perfect"] = (total_perfect /
                    float(total_reads))
                # What portion are mapped at all?
                sample_stats["portion_mapped_at_all"] = (total_mapped_at_all / 
                    float(total_reads))
                # What was the portion with no indels?
                sample_stats["portion_no_indels"] = (total_no_indels / 
                    float(total_reads))
                # And the potion mapped with indels (as opposed to unmapped)
                # TODO: what part of this is horrible mappings with tiny scores
                # overall?
                sample_stats["portion_with_indels"] = (total_with_indels / 
                    float(total_reads))
                # How many have one indel exactly?
                sample_stats["portion_one_indel"] = (total_one_indel / 
                    float(total_reads))    
                # And the portion with no substitutions
                sample_stats["portion_no_substitutions"] = \
                    (total_no_substitutions / float(total_reads))
                
                if total_mapped_at_all > 0:
                    # Some things we want to divide by the mapped reads
                
                    # What's the average matches per column for primary alignments
                    # for reads that actually aligned?
                    sample_stats["mean_matches_per_column"] = \
                        (total_matches_per_column / float(total_mapped_at_all))
                    # What's the average score for primary alignments of reads that
                    # actually aligned?
                    sample_stats["mean_score"] = \
                        (total_score / float(total_mapped_at_all))
                else:
                    # If no reads actually mapped we can just put 0
                    
                    sample_stats["mean_matches_per_column"] = 0
                    sample_stats["mean_score"] = 0
                    
                # What was the runtime?
                sample_stats["runtime"] = runtime
                
                try:
                    # See if we have access to these extra stats
                    
                    # How many total bases of reads have primary alignments (to
                    # non-Ns)?
                    total_aligned = sum((count * float(weight)
                        for weight, count in 
                        stats["aligned_lengths"].iteritems()))
                        
                    # The aligned bases are the ones that aren't in indels or
                    # leading/trailing softclips.
                    
                    # Calculate portion of aligned bases that are substitutions.
                    # Doesn't count indels at all.
                    sample_stats["substitution_rate"] = (substitution_bases /
                        float(total_aligned))
                        
                    # Calculate indels per base
                    sample_stats["indel_rate"] = (indel_instances /
                        float(total_aligned))
                        
                except:
                    # Sometimes we just won't have these stats available
                    pass
                    
        
        # Now save all these portion stats we extracted back to the cache
        writer = tsv.TsvWriter(open(local_filename, "w"))
        
        for graph, stats_by_sample in stats_cache.iteritems():
            # For each graph and all the stats for that graph
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                for stat_name, stat_value in stats_by_name.iteritems():
                    # For each stat
                    
                    # Save each stat value
                    writer.line(graph, sample, stat_name, stat_value)
                
        # Close the file and save the results for the next run
        writer.close()
        out_store.write_output_file(local_filename, cache_tsv_key)        
    
    # We want normalized and un-normalized versions of the stats cache
    stats_by_mode = {"absolute": stats_cache}
    
    if stats_cache.has_key("refonly"):
        
        # Deep copy and normalize the stats cache
        normed_stats_cache = copy.deepcopy(stats_cache)
        
        # We want to normalize and the reference exists (i.e. not CENX)
        # Normalize every stat against the reference
        for graph, stats_by_sample in normed_stats_cache.iteritems():
            # For each graph and all the stats for that graph
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                for stat_name in stats_by_name.keys():
                
                    if stats_cache["refonly"].has_key(sample):
                
                        # Get the reference value
                        ref_value = stats_cache["refonly"][sample][stat_name]
                        
                        if ref_value != 0:
                            
                            # Normalize
                            stats_by_name[stat_name] /= ref_value
                        else:
                            stats_by_name[stat_name] = None
                    else:
                        # Nothing to norm against. TODO: maybe complain when
                        # sample sets aren't all the same?
                        stats_by_name[stat_name] = None
                        
        # Register this as a condition
        stats_by_mode["normalized"] = normed_stats_cache
    
    # Now break out the stats from our massive cached dict into files for
    # plotting
        
    for mode, mode_stats_cache in stats_by_mode.iteritems():
        
        # We need some config
        # Where should we route each stat to?
        stat_file_keys = {
            "portion_mapped_well": "plots/{}/mapping.{}.tsv".format(mode,
                region),
            "portion_perfect": "plots/{}/perfect.{}.tsv".format(mode, region),
            "portion_single_mapped_well":
                "plots/{}/singlemapping.{}.tsv".format(mode, region),
            "portion_single_mapped_at_all":
                "plots/{}/singlemappingatall.{}.tsv".format(mode, region),
            "portion_mapped_at_all": "plots/{}/anymapping.{}.tsv".format(mode,
                region),
            "runtime": "plots/{}/runtime.{}.tsv".format(mode, region),
            "portion_no_indels": "plots/{}/noindels.{}.tsv".format(mode,
                region),
            "portion_one_indel": "plots/{}/oneindel.{}.tsv".format(mode,
                region),
            "portion_with_indels": "plots/{}/hasindels.{}.tsv".format(mode,
                region),
            "portion_no_substitutions": "plots/{}/nosubsts.{}.tsv".format(mode,
                region),
            "substitution_rate": "plots/{}/substrate.{}.tsv".format(mode,
                region),
            "indel_rate": "plots/{}/indelrate.{}.tsv".format(mode,
                region),
            "mean_matches_per_column": "plots/{}/matches.{}.tsv".format(mode,
                region),
            "mean_score": "plots/{}/score.{}.tsv".format(mode, region)
        }
        
        # Make a local temp file for each (dict from stat name to file object
        # with a .name).
        stats_file_temps = {name: tempfile.NamedTemporaryFile(
            dir=job.fileStore.getLocalTempDir(), delete=False) 
            for name in stat_file_keys.iterkeys()}
            
        for graph, stats_by_sample in mode_stats_cache.iteritems():
            # For each graph and all the stats for that graph
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                for stat_name, stat_value in stats_by_name.iteritems():
                    # For each stat
                    
                    if not stats_file_temps.has_key(stat_name):
                        # Skip stats that have nowhere to go
                        continue
                    
                    # Write graph and value to the file for the stat, for
                    # plotting
                    stats_file_temps[stat_name].write("{}\t{}\n".format(
                        graph, stat_value))
        
        for stat_name, stat_file in stats_file_temps.iteritems():
            # Flush and close the temp file
            stat_file.flush()
            os.fsync(stat_file.fileno())
            stat_file.close()
            
            # Upload the file
            out_store.write_output_file(stat_file.name,
                stat_file_keys[stat_name])
                
        # Special case: perfect vs unique mapping
        # perfect = precision = y, unique ~= recall = x
        perfect_vs_unique = tempfile.NamedTemporaryFile(
            dir=job.fileStore.getLocalTempDir(), delete=False)
        
        for graph, stats_by_sample in mode_stats_cache.iteritems():
            # For each graph and all the stats for that graph
            
            # We're going to make lists and compute medians
            all_perfect = []
            all_unique = []
            
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                
                # Get the prefect and unique stats
                perfect = stats_by_name["portion_perfect"]
                unique = stats_by_name["portion_single_mapped_well"]
                
                all_perfect.append(perfect)
                all_unique.append(unique)
                
            # Write them in x, y order
            perfect_vs_unique.write("{}\t{}\t{}\n".format(graph,
                numpy.median(all_unique),
                numpy.median(all_perfect)))
        
        perfect_vs_unique.flush()            
        os.fsync(perfect_vs_unique.fileno())
        perfect_vs_unique.close()
           
        # Save the perfect vs unique data
        out_store.write_output_file(perfect_vs_unique.name,
                "plots/{}/perfect_vs_unique.{}.tsv".format(mode, region))
        
    # Return the cached stats. TODO: break out actual collate and upload into
    # its own function that will also do normalization.
    # Make sure to use toillib to strip off the defaultdict lambdas
    return de_defaultdict(stats_cache)
    
def merge_overall(job, options, stats_by_region):
    """
    Given stats by region, graph, sample, and stat name, produce the data files
    needed for overall plots that aren't just concatenations of the individual
    region plots.
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Right now all we need to do is the perfect vs unique plot.
    
    # Make a file to write the points to plot to.
    # perfect = precision = y, unique ~= recall = x
    perfect_vs_unique = tempfile.NamedTemporaryFile(
        dir=job.fileStore.getLocalTempDir(), delete=False)
    
    # For each sample in each graph, we're going to sum total reads, total
    # perfect, and total unique across all regions. Then we're going to compute
    # the perfect and unique values for all the samples, find the median in each
    # dimension, and plot that point. TODO: do a cool cloud instead of a point.
    
    # Holds perfect read count by graph and sample
    perfect_count = collections.defaultdict(collections.Counter)
    # Holds unique read count by graph and sample
    unique_count = collections.defaultdict(collections.Counter)
    # Holds overall read count by graph and sample
    read_count = collections.defaultdict(collections.Counter)
    
    # Keep a set of all observed graph names
    graph_names = set()
    
    for region, stats_by_graph in stats_by_region.iteritems():
        # For each region and its stats for the graphs
        for graph, stats_by_sample in stats_by_graph.iteritems():
            # For each graph and all the stats for that graph
            
            # We have a graph with this name
            graph_names.add(graph)
            
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                
                # Add in perfect reads
                perfect_count[graph][sample] += stats_by_name["perfect"]
                # And (sufficiently) unbique reads
                unique_count[graph][sample] += \
                    stats_by_name["single_mapped_well"]
                # And total reads
                read_count[graph][sample] += stats_by_name["total_reads"]
    
    for graph in graph_names:           
        # For every graph under consideration 
        
        # We'll fill this with perfect and unique portions for all the samples.
        perfect_portions = []
        unique_portions = []
        
        for sample in perfect_count[graph].iterkeys():
            # Now calculate per-sample perfect and unique portions
            perfect_portions.append(float(perfect_count[graph][sample]) /
                read_count[graph][sample])
            unique_portions.append(float(unique_count[graph][sample]) /
                read_count[graph][sample])
        
        # Write them in x, y order
        perfect_vs_unique.write("{}\t{}\t{}\n".format(graph,
            numpy.median(unique_portions),
            numpy.median(perfect_portions)))
    
    perfect_vs_unique.flush()            
    os.fsync(perfect_vs_unique.fileno())
    perfect_vs_unique.close()
       
    # Save the perfect vs unique data
    # TODO: we don't support normalization for this overall plot
    out_store.write_output_file(perfect_vs_unique.name,
            "plots/absolute/perfect_vs_unique.ALL.tsv")
    
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    RealTimeLogger.start_master()
    
    # Make a root job
    root_job = Job.wrapJobFn(collate_all, options,
        cores=1, memory="1G", disk="1G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

