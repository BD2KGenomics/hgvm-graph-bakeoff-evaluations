#!/usr/bin/env python2.7
"""
biasDetector.py: see if there is a population bias in mapping rates for each
graph.

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
    parser.add_argument("in_store",
        help="input IOStore to find stats files in (under /stats)")
    parser.add_argument("out_store",
        help="output IOStore to put collated plotting files in (under /bias)")
    parser.add_argument("--blacklist", action="append", default=[],
        help="ignore the specified region:graph pairs")
    parser.add_argument("--index_url", 
        default=("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000_genomes_project/1000genomes.sequence.index"), 
        help="URL to index of samples, with SAMPLE_NAME and POPULATION_NAME")
    parser.add_argument("--superpopulation_url",
        default="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/"
        "20131219.populations.tsv",
        help="URL to index of superpopulation assignments")
    parser.add_argument("--overwrite", action="store_true",
        help="replace cached per-sample statistics with recalculated ones")
    
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
   
def scan_all(job, options):
    """
    Scan all the regions and graphs for bias.
    
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
            cores=1, memory="1G", disk="10G")
    
def scan_region(job, options, region, pop_by_sample):
    """
    Scan all the graphs in a region for bias.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Will hold stats by population (promises) by graph name
    graph_stats = {}
    
    for graph in in_store.list_input_directory("stats/{}".format(region)):
        # For each
        
        if "{}:{}".format(region, graph) in options.blacklist:
            # We don't want to process this region/graph pair.
            RealTimeLogger.get().info("Skipping {} graph {}".format(region,
                graph))
            continue
        
        # Scan the graph and get its per-population stats
        stats_for_graph = job.addChildJobFn(scan_graph, options, region, graph,
            pop_by_sample, cores=1, memory="1G", disk="10G").rv()
            
        # Save the stats for the graph
        graph_stats[graph] = stats_for_graph
        
    # When stats are all in, look at the, and save the results
    job.addFollowOnJobFn(save_region_stats, options, region, graph_stats,
        cores=1, memory="1G", disk="10G")
        
def scan_graph(job, options, region, graph, pop_by_sample):
    """
    Look at a given graph in a given region and return its mismatch portions by
    population, sorted by sample name.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # This holds lists of portion perfectly-mapped stat values by population
    stats_by_pop = collections.defaultdict(list)
    
    # This is the output TSV it should all land in, in <pop>\t<value> format
    labeled_tsv_key = "bias/distributions/{}/{}.tsv".format(region, graph)
    
    # What name will it have locally for us?
    local_filename = os.path.join(job.fileStore.getLocalTempDir(), "temp.tsv")
    
    if out_store.exists(labeled_tsv_key) and not options.overwrite:
        # Just read in from that TSV
        
        # Grab the cached results
        out_store.read_input_file(labeled_tsv_key, local_filename)
        
        # Read all the pop, value pairs from the TSV
        reader = tsv.TsvReader(open(local_filename))
        
        for pop, stat in reader:
            # Collect all the values and put them in the right populations
            stats_by_pop[pop].append(float(stat))
        
    else:
        # Look in each sample file, compute the per-sample statistic, and save
        # it to the tsv

        RealTimeLogger.get().info("Processing samples for {} graph {}".format(
            region, graph))
    
        for result in in_store.list_input_directory("stats/{}/{}".format(region,
            graph)):
            
            # Pull sample name from filename
            sample_name = re.match("(.*)\.json$", result).group(1)
        
            # Grab file
            json_filename = os.path.join(job.fileStore.getLocalTempDir(),
                "temp.json")
            in_store.read_input_file("stats/{}/{}/{}".format(region, graph,
                result), json_filename)
            # Read the JSON
            stats = json.load(open(json_filename))
            
            # How many reads are mapped perfectly?
            total_perfect = sum((stats["primary_mismatches"].get(str(x), 0)
                    for x in xrange(1)))
                
            # How many reads are there overall for this sample?
            total_reads = stats["total_reads"]
            
            # What portion are mapped perfectly?
            portion_perfect = total_perfect / float(total_reads)
            
            # Add the sample to the distribution, with the name in for sorting
            stats_by_pop[pop_by_sample[sample_name]].append((sample_name,
                portion_perfect))
                
        for pop_name in stats_by_pop.keys():
            # Sort all the distributions by sample name, and throw out the
            # sample names
            stats_by_pop[pop_name] = [value for sample_name, value in
                sorted(stats_by_pop[pop_name])]
            
        # Write all the pop, value pairs to the TSV
        writer = tsv.TsvWriter(open(local_filename, "w"))
        
        for pop, list_of_stats in stats_by_pop.iteritems():
            for stat in list_of_stats:
                # Save each sample for its population
                writer.line(pop, stat)
                
        # Close the file and save the results for the next run
        writer.close()
        out_store.write_output_file(local_filename, labeled_tsv_key)
        
        
    return stats_by_pop
    
def save_region_stats(job, options, region, graph_stats):
    """ 
    Save the biases of the graphs for this region to the output IOStore.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    # Now subtract refOnly out of everything as our normalization.
    
    if graph_stats.has_key("refonly"):
        # Normalize stat against the refonly graph
        
        # Pull out the refonly stats dict, as a deep copy so we don't
        # clobber it. Do it before we manage to normalize refonly itself.
        ref_stats = copy.deepcopy(graph_stats["refonly"])
        
        for graph, stats_by_pop in graph_stats.iteritems():
            
            for pop_name in stats_by_pop.keys():
                # Normalize each population
                
                # Thest teo lists need to correspond
                assert(len(stats_by_pop[pop_name]) == len(ref_stats[pop_name]))
                
                # Zip the two stats lists together and do the division, and
                # replace the non-reference list.
                stats_by_pop[pop_name] = [this_stat - ref_stat for 
                    (this_stat, ref_stat) in itertools.izip(
                    stats_by_pop[pop_name], ref_stats[pop_name])]
                    
            # Save the normalized values
            normed_tsv_key = "bias/normalized_distributions/{}/{}.tsv".format(
                region, graph)
        
            # What name will it have locally for us?
            normed_file = tempfile.NamedTemporaryFile(
                dir=job.fileStore.getLocalTempDir(), delete=False)
                
            # Make a TSV writer
            normed_writer = tsv.TsvWriter(normed_file)
            
            for pop_name, value_list in stats_by_pop.iteritems():
                for value in value_list:
                    # Dump each normalized value with its population
                    normed_writer.line(pop_name, value)
                    
            # Finish up our local temp file
            normed_writer.close()
            
            # Save it to the IOStore. TODO: write stream methods for this!
            out_store.write_output_file(normed_file.name, normed_tsv_key)
                        
    else:
        RealTimeLogger.get().warning("Can't normalize {}".format(region))
        
    # After normalization, do the statistical test
    
    # Grab file to save the overall bias levels for the region in
    local_filename = os.path.join(job.fileStore.getLocalTempDir(),
        "temp.tsv")
        
    # Get a writer to write to it
    writer = tsv.TsvWriter(open(local_filename, "w"))
    
    for graph, stats_by_pop in graph_stats.iteritems():
    
        RealTimeLogger.get().info("Running statistics for {} graph {}".format(
            region, graph))
    
        # Grab all the distributions to compare in a list
        list_of_distributions = stats_by_pop.values()

        RealTimeLogger.get().info("Have {} populations to compare".format(len(
            list_of_distributions)))
        
        try:
            
            # Now we have the data for each graph read in, so we can run stats.
            # Test to see if there is a significant difference in medians among
            # the populations.
            h_statistic, p_value = scipy.stats.mstats.kruskalwallis(
                *list_of_distributions)
                
        except ValueError:
            # We couldn;'t run the test on these particular numbers.
            h_statistic = None
            p_value = None
            
        # Since there probably is, quantify the degree of difference between
        # populations. This is my own metric which may or may not be good.
        
        # Get the median of each distribution
        medians = [numpy.median(numpy.array(l)) for l in list_of_distributions]
        
        # Find the standard deviation and call it the bias level
        bias_level = numpy.std(medians)
    
        # Write this stat to the TSV
        writer.line(graph, h_statistic, p_value, bias_level)
        
    # Finish up the file
    writer.close()
    
    RealTimeLogger.get().info("Uploading results for {} ".format(region))
    
    # Save it to the output store
    out_store.write_output_file(local_filename, "bias/{}.tsv".format(region))
            
    
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
    root_job = Job.wrapJobFn(scan_all, options,
        cores=1, memory="1G", disk="1G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

