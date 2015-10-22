#!/usr/bin/env python2.7
"""
parallelAzureDownloader.py: download files from Azure.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import tempfile

from toil.job import Job

from toillib import *

from multiprocessing.pool import ThreadPool

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
        help="input IOStore to download from")
    parser.add_argument("out_store",
        help="output IOStore to put things in")
    parser.add_argument("--overwrite", default=False, action="store_true",
        help="overwrite existing files")
    parser.add_argument("--batch_size", type=int, default=1000,
        help="number of files to copy in a batch")
    
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

def copy_everything(job, options):
    """
    Download the file list and copy all the files.
    
    """
    
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)
    
    batch_count = 0;
    
    for batch in group(in_store.list_input_directory("", recursive=True), options.batch_size):
        # For every batch, strip out the Nones
        batch = [x for x in batch if x is not None]
    
        # Copy everything in that batch
        job.addChildJobFn(copy_batch, options, batch, cores=1, memory="1G",
            disk="10G")
            
        batch_count += 1
        
        if batch_count % 10 == 0:
            
            RealTimeLogger.get().info("Queued {} batches...".format(batch_count))
    
def copy_batch(job, options, batch):
    """
    Copy a batch of files from input to output.
    """
        
    # Set up the IO stores.
    in_store = IOStore.get(options.in_store)
    out_store = IOStore.get(options.out_store)

    # Start some threads
    pool = ThreadPool(10)
    
    
    def download(filename):
        """
        Download each file
        """
        
        if (not options.overwrite) and out_store.exists(filename):
            # Skip existing file
            return
        
        # Make a temp file
        (handle, path) = tempfile.mkstemp(dir=job.fileStore.getLocalTempDir())
        os.close(handle)        
        
        # Download
        in_store.read_input_file(filename, path)
        # Store
        out_store.write_output_file(path, filename)
        
        # Clean up
        os.unlink(path)
        
    # Run all the downloads in parallel
    pool.map(download, batch)
    
        
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
    root_job = Job.wrapJobFn(copy_everything, options,
        cores=1, memory="1G", disk="1G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

