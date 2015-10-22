#!/usr/bin/env python2.7
"""
collateStatistics.py: turn individual .stats files into input files for plotting.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import tempfile

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
        help="output IOStore to put collated plotting files in (under /plots)")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
   
def group(iterable, max_count):
    """
    Batch up iterable results. Pads with None.
    
    See <http://stackoverflow.com/a/8290490/402891>
    
    """
    
    # Zip a bunch of copies of the iterable.
    return itertools.izip_longest(*([iter(iterable)] * max_count))

def collate_all(job, options):
    """
    Collate all the stats files
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    
    for region in in_store.list_input_directory("stats"):
        # Collate everything in the region
        job.addChildJobFn(collate_region, options, region, cores=1, memory="1G",
            disk="10G")
    
    
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
            with file_store.readGlobalFileStream(input_id) as input_stream:
                # Copy over every input
                shutil.copyfileobj(input_stream, destination)
    
    # Save the result
    out_store.write_output_file(path, output_key)

    
    
            
def collate_region(job, options, region):
    """
    Collate all the stats files in a region
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    
    # Where will we put our final collated files?
    mapping_key="plots/mapping.{}.tsv".format(region)
    perfect_key="plots/perfect.{}.tsv".format(region)
    one_error_key="plots/oneerror.{}.tsv".format(region)
    single_mapping_key="plots/singlemapping.{}.tsv".format(region)
    runtime_key="plots/runtime.{}.tsv".format(region)
    
    # Where will we build them?
    mapping_temp_name = os.path.join(job.fileStore.getLocalTempDir(),
        "mapping.tsv")
    mapping_temp = open(mapping_temp_name, "w")
    perfect_temp_name = os.path.join(job.fileStore.getLocalTempDir(),
        "perfect.tsv")
    perfect_temp = open(perfect_temp_name, "w")
    one_error_temp_name = os.path.join(job.fileStore.getLocalTempDir(),
        "oneerror.tsv")
    one_error_temp = open(one_error_temp_name, "w")
    single_mapping_temp_name = os.path.join(job.fileStore.getLocalTempDir(),
        "singlemapping.tsv")
    single_mapping_temp = open(single_mapping_temp_name, "w")
    runtime_temp_name = os.path.join(job.fileStore.getLocalTempDir(),
        "runtime.tsv")
    runtime_temp = open(runtime_temp_name, "w")
    
    for graph in in_store.list_input_directory("stats/{}".format(region)):
        # For each 
        RealTimeLogger.get().info("Processing {} graph {}".format(region,
            graph))
        for result in in_store.list_input_directory("stats/{}/{}".format(region,
            graph)):
        
            # Grab file
            local_filename = os.path.join(job.fileStore.getLocalTempDir(),
                "temp.json")
            in_store.read_input_file("stats/{}/{}/{}".format(region, graph,
                result), local_filename)
            # Read the JSON
            stats = json.load(open(local_filename))
            
            # Compute the answers for this sample
            
            # How many reads are mapped well enough?
            total_mapped = sum((stats["primary_mismatches"].get(str(x), 0)
                for x in xrange(3)))
                
            # How many reads are multimapped well enough?
            total_multimapped = sum((stats["secondary_mismatches"].get(str(x),
                0) for x in xrange(3)))
                
            # How many reads are perfect?
            total_perfect = sum((stats["primary_mismatches"].get(str(x), 0)
                for x in xrange(1)))
                
            # How many reads are mapped with <=1 error?
            total_one_error = sum((stats["primary_mismatches"].get(str(x), 0)
                for x in xrange(2)))
                
            # How many reads are there overall for this sample?
            total_reads = stats["total_reads"]
            
            # What was the runtime?
            runtime = stats.get("run_time", None)
            if runtime is None:
                # We need NaN floats if there's no runtime
                runtime = float("nan")
            else:
                # Convert to time per read aligned
                runtime /= total_reads
            
            # What portion are single-mapped?
            portion_single_mapped = (total_mapped - total_multimapped) / float(
                total_reads)
            # What portion are mapped?
            portion_mapped = total_mapped / float(total_reads)
            # What portion are perfect?
            portion_perfect = total_perfect / float(total_reads)
            # What portion are <=1 error?
            portion_one_error = total_one_error / float(total_reads)

            # Append lines to the files            
            mapping_temp.write("{}\t{}\n".format(graph, portion_mapped))
            perfect_temp.write("{}\t{}\n".format(graph, portion_perfect))
            one_error_temp.write("{}\t{}\n".format(graph, portion_one_error))
            single_mapping_temp.write("{}\t{}\n".format(graph,
                portion_single_mapped))
            runtime_temp.write("{}\t{}\n".format(graph, runtime))
            
    # Flush and close the files
    mapping_temp.close()
    perfect_temp.close()
    one_error_temp.close()
    single_mapping_temp.close()
    runtime_temp.close()
    
    # Save them to the final places
    out_store.write_output_file(mapping_temp_name, mapping_key)
    out_store.write_output_file(perfect_temp_name, perfect_key)
    out_store.write_output_file(one_error_temp_name, one_error_key)
    out_store.write_output_file(single_mapping_temp_name, single_mapping_key)
    out_store.write_output_file(runtime_temp_name, runtime_key)
    
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
        
        
        
        
        
        
        
        
        
        

