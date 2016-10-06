#!/usr/bin/env python2.7
"""
biasDetector.py: see if there is a population bias in mapping rates for each
graph.

Splits out certain statistics from collateStatistics.py by superpopulation, for
plotting with plotBias.sh.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import tempfile
import urllib2
import copy

import tsv
import numpy
import scipy.stats
import scipy.stats.mstats

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
        help="output IOStore to put collated plotting files in (under /bias)")
    parser.add_argument("--blacklist", action="append", default=[],
        help="ignore the specified region:graph pairs")
    parser.add_argument("--samples", type=argparse.FileType("r"),
        help="limit to samples on this list, one per line")
    parser.add_argument("--index_url", 
        default=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000_genomes_project/1000genomes.sequence.index"), 
        help="URL to index of samples, with SAMPLE_NAME and POPULATION_NAME")
    parser.add_argument("--superpopulation_url",
        default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/"
        "20131219.populations.tsv",
        help="URL to index of superpopulation assignments")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def tsv_reader_with_comments(lines):
    """
    TSV reader iterator that doesn't strip comments.
    """
    
    for line in lines:
        yield line.split("\t")
   
def scan_all(job, options, sample_whitelist):
    """
    Scan all the regions and graphs for bias.
    
    Only looks at samples in the whitelist set, if the whitelist is not None.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Download the superpopulation assignments
    # This holds superpop by pop
    superpopulation_by_population = {}
    
    for parts in tsv.TsvReader(urllib2.urlopen(urllib2.Request(
        options.superpopulation_url))):
        # For each population code (column 1), assign it to the right
        # superpopulation (column 2).
        superpopulation_by_population[parts[1]] = parts[2]
        
    
    RealTimeLogger.get().info("Downloading sample population assignments")    
    
    # Load the 1000 Genomes population assignments.
    # Make a reader that goes through split out lines in the TSV.
    reader = tsv_reader_with_comments(urllib2.urlopen(urllib2.Request(
        options.index_url)))
    
    # Get an iterator over the lines
    lines = iter(reader)
    
    # Grab the headings
    headings = lines.next()
    
    while headings[0].startswith("##"):
        # Skip leading lines that aren't the real header (which starts with #)
        headings = lines.next()
    
    # Which column holds sample names?
    sample_name_column = headings.index("SAMPLE_NAME")
    
    # Which column holds sample populations?
    sample_population_column = headings.index("POPULATION")
    
    # What dict do we fill in? Holds population string by sample name.
    # We now use the superpopulation names for our populations.
    pop_by_sample = {}
    
    # We also want to count samples in each population for debuging
    samples_per_pop = collections.Counter()
    
    for parts in lines:
        # Save superpopulation under sample
        pop_by_sample[parts[sample_name_column]] = \
            superpopulation_by_population[parts[sample_population_column]]
            
        # Count the sample for its population
        samples_per_pop[parts[sample_population_column]] += 1
        
    RealTimeLogger.get().info("Found {} populations:".format(len(
        samples_per_pop))) 
        
    for (pop, count) in samples_per_pop.iteritems():
        RealTimeLogger.get().info("{}: {}".format(pop, count))
    
    for region in in_store.list_input_directory("stats"):
        # Collate everything in the region
        job.addChildJobFn(scan_region, options, region, pop_by_sample,
            sample_whitelist, cores=1, memory="1G", disk="10G")
    
def scan_region(job, options, region, pop_by_sample, sample_whitelist):
    """
    Scan all the graphs in a region for bias.
    
    If sample_whitelist is not None, ignores samples not in that set.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # This holds a dict from graph name, then sample name, then stat name to
    # actual stat value.
    stats_cache = collections.defaultdict(lambda: collections.defaultdict(dict))
    
    # This is the cache file for this region, in
    # <graph>\t<sample>\t<stat>\t<value> format
    cache_tsv_key = "plots/cache/{}.tsv".format(region)
    
    # What name will it have locally for us?
    local_filename = os.path.join(job.fileStore.getLocalTempDir(), "temp.tsv")
    
    if out_store.exists(cache_tsv_key):
        # Just read in from that TSV
        
        RealTimeLogger.get().info("Loading cached region {}".format(region))
        
        # Grab the cached results
        out_store.read_input_file(cache_tsv_key, local_filename)
        
        # Read all the pop, value pairs from the TSV
        reader = tsv.TsvReader(open(local_filename))
        
        # Which samples are going to be skipped?
        skipped_samples = set()
        
        for graph, sample, stat, value in reader:
            # Read every line from the cache and pull out what value for what
            # stat it gives for what sample.
            
            if sample_whitelist is not None and sample not in sample_whitelist:
                # Skip this sample that's not on the list
                skipped_samples.add(sample)
                continue
        
            # Populate our cache dict
            stats_cache[graph][sample][stat] = float(
                value)
                
        RealTimeLogger.get().info("Skipped {} samples".format(
            len(skipped_samples)))
            
    else:
        # Stats haven't been collated
        raise RuntimeError(
            "No graph stats for {}; run collateStatistics.py".format(region))
            
    # We want normalized and un-normalized versions of the stats cache
    stats_by_mode = {"absolute": stats_cache}
    
    if stats_cache.has_key("refonly"):
        
        # Deep copy and normalize the stats cache
        normed_stats_cache = copy.deepcopy(stats_cache)
        
        # We want to normalize and the reference exists (i.e. not CENX)
        # Normalize every stat against the reference, by subtraction
        for graph, stats_by_sample in normed_stats_cache.iteritems():
            # For each graph and all the stats for that graph
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                for stat_name in stats_by_name.keys():
                
                    if stats_cache["refonly"].has_key(sample):
                
                        # Get the reference value
                        ref_value = stats_cache["refonly"][sample][stat_name]
                        
                        
                        # Normalize by subtraction
                        stats_by_name[stat_name] -= ref_value
                        
                    else:
                        # Nothing to norm against. TODO: maybe complain when
                        # sample sets aren't all the same?
                        stats_by_name[stat_name] = None
                        
        # Register this as a condition
        stats_by_mode["normalized"] = normed_stats_cache
            
    
    # Now save stats, parceling out by region and graph
        
    for mode, mode_stats_cache in stats_by_mode.iteritems():
        for graph, stats_by_sample in mode_stats_cache.iteritems():
            # We need some config
            # Where should we route each stat to?
            stat_file_keys = {
                "substitution_rate": "bias/{}/{}/substrate.{}.tsv".format(mode,
                    region, graph),
                "indel_rate": "bias/{}/{}/indelrate.{}.tsv".format(mode,
                    region, graph),
                "portion_perfect": "bias/{}/{}/perfect.{}.tsv".format(mode,
                    region, graph)
            }
        
            # Make a local temp file for each (dict from stat name to file
            # object with a .name).
            stats_file_temps = {name: tempfile.NamedTemporaryFile(
                dir=job.fileStore.getLocalTempDir(), delete=False) 
                for name in stat_file_keys.iterkeys()}
                
            for sample, stats_by_name in stats_by_sample.iteritems():
                # For each sample and all the stats for that sample
                for stat_name, stat_value in stats_by_name.iteritems():
                    # For each stat
                    
                    if not stats_file_temps.has_key(stat_name):
                        # Skip stats that have nowhere to go
                        continue
                    
                    # Write graph and value to the file for the stat, for
                    # plotting, naming it after the pop that the sample is in
                    stats_file_temps[stat_name].write("{}\t{}\n".format(
                        pop_by_sample[sample], stat_value))
            
            for stat_name, stat_file in stats_file_temps.iteritems():
                # Flush and close the temp file
                stat_file.flush()
                os.fsync(stat_file.fileno())
                stat_file.close()
                
                # Upload the file
                out_store.write_output_file(stat_file.name,
                    stat_file_keys[stat_name])
        
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
    
    # Load the sample whitelist, if applicable. Holds a set if we have a
    # whitelist, or None otherwise.
    sample_whitelist = None
    if options.samples is not None:
        # Read all the samples from the file
        sample_whitelist = set([line[0] for line in
            tsv.TsvReader(options.samples)])
    
    RealTimeLogger.start_master()
    
    # Make a root job
    root_job = Job.wrapJobFn(scan_all, options, sample_whitelist,
        cores=1, memory="1G", disk="1G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

