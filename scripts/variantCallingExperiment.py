#!/usr/bin/env python2.7
"""
Create pileups from GAMs, and evaluate different parameter sets for the pileup-
to-glennfile and glenn-to-vcf steps of the variant calling pipeline.

Takes alignments stored like this:
  alignments/brca1/cactus/NA19240.gam

so:
  <alignments IOstore>/<region>/<graph>/<sample>.gam

The graphs must also be accessible:

  graphs/cactus-brca1.vg

so:
  <graphs IOStore>/<graph method>-<region>.vg

note: debruin is exception where -k63 tag gets tacked on at end (treated as special case)


"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
import hashlib
import signal
from threading import Timer
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs, IOStore

class ExperimentCondition:
    """
    Represents a combination of read filtering, graph augmenting, vcf
    generating, and vcf filtering options.
    """
    
    def __init__(self, read_filter_options, pileup_options, call_options, vcf_options):
        """
        Take dicts from option string to option value for read filtering,
        pileup-ing, graph augmenting, and vcf conversion, and produces an
        ExperimentCondition wrapping them. Option values may only be single
        strings, and an option may only occur a single time (since we have a
        dict). TODO: also add VCF filtering
        
        """
        
        self.read_filter_options = read_filter_options
        self.pileup_options = pileup_options
        self.call_options = call_options
        self.vcf_options = vcf_options
        
    def dict_to_string(self, options_dict):
        """
        Convert a dict of option values to a single string.
        """
        
        # Get a nested list of pairs of options and value strings. Use sorting
        # to make sure we get them in a deterministic order.
        nested = ([opt_name, str(opt_value)] for opt_name, opt_value in sorted(options_dict.iteritems()))
        
        # Flatten and join into a string
        return " ".join((item for pair in nested for item in pair))
        
    def string_to_path(self, string):
        """
        Convert a string into a usable path component (directory name).
        """
        
        # Condense spaces to underscores, and remove all but a few safe
        # characters.
        return re.sub('[^a-zA-Z0-9_!.]', "", re.sub('\s+', "_", string).strip("_"))
        
    def get_read_filter_options(self):
        """
        Return the options string for read filtering.
        """
        return self.dict_to_string(self.read_filter_options)
        
    def get_pileup_options(self):
        """
        Return the options string for pileup-ing.
        """
        return self.dict_to_string(self.pileup_options)
        
    def get_call_options(self):
        """
        Return the options string for vg call.
        """
        return self.dict_to_string(self.call_options)
        
    def get_vcf_options(self):
        """
        Return the options string for glenn2vcf.
        """
        return self.dict_to_string(self.vcf_options)
        
    def get_pileup_condition_name(self):
        """
        Return a string that can be part of a filesystem path and which depends
        on all the parameters that affect the pielup.
        """
        
        # Depends on the filter and pileup options
        return self.string_to_path(self.get_read_filter_options() + "!" + self.get_pileup_options())
    
    def get_glennfile_condition_name(self):
        """
        Return a string that can be part of a filesystem path and which depends
        on all the parameters that affect the glenn file.
        """
        
        # Depends on the pileup file and the call options
        return self.string_to_path(self.get_pileup_condition_name() + "!" + self.get_call_options())
        
    def get_vcf_condition_name(self):
        """
        Return a string that can be part of a filesystem path and which depends
        on all the parameters that affect the final VCF file.
        """
        
        # Depends on the gelnnfile and the glenn2vcf options
        return self.string_to_path(self.get_glennfile_condition_name() + "!" + self.get_vcf_options())
        
        
        
def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("in_gams",
        help="input alignment files IOStore")
    parser.add_argument("in_graphs",
        help="input graph files IOStore")
    parser.add_argument("cache",
        help="cache IOStore so we don't have to re-do pileups constantly")
    parser.add_argument("out_dir",
        help="output IOStore, containing experimental results")
    parser.add_argument("--truth",
        help="IOStore to search for truth vcfs for comparison. For each region "
        "we expect <region>/<region>.vcf")
    parser.add_argument("--references",
        help="IOStore to search for reference FASTAs. For each region we "
        "expect <region>/ref.fa")
    parser.add_argument("--regions",
        help="IOStore to search for region BEDs. For each region we "
        "expect <REGION>.bed, with the region in upper-case.")
    parser.add_argument("--blacklist", nargs="+", default=["vglr", "mhc:camel",
        "lrc_kir:camel", "debruijn-k31", "debruijn-k63", "cenx", "simons",
        "curoverse"],
        help="ignore the specified regions, graphs, or region:graph pairs")
    
    args = args[1:]
        
    return parser.parse_args(args)

def run(cmds, stdout = sys.stdout, stderr = sys.stderr, timeout_sec = sys.maxint,
        timeout_dep = None, fail_hard = False):
    """
    Run commands in the given list in the shell, piping each in the next. Throw
    an exception if any of the commands fails or if the whole pipeline times
    out.
    
    
    """
    
    
    def timeout_fail(procs, cmds):
        """
        Called when the given processes, launched from the given commands, have
        run for too long. Kills the processes and logs an error.
        
        """
        
        for proc in procs:
            os.kill(proc.pid, signal.SIGKILL)
            proc.kill()
            
        # we often check to see if some output file exists before running
        # if we timeout, make sure the file exists so rerunning wont timeout again
        # at same place (unless overwrite explicitly desired)
        if timeout_dep is not None and os.path.exists(timeout_dep):
            os.system("rm -rf {}; echo timeout > {}".format(timeout_dep, timeout_dep))
        if fail_hard is True:
            raise RuntimeError("Command: {} timed out".format(" | ".join(cmds)))
        else:
            RealTimeLogger.get().warning("Command: {} timed out".format(" | ".join(cmds)))
    
    RealTimeLogger.get().info("RUN: {}".format(" | ".join(cmds)))

    # We have a list of processes, one per command.
    procs = []
    # We remember the previous process's standard output
    last_stdout = None

    for cmd in cmds[:-1]:
        # All but the last command feed their standard output into pipes
        proc = subprocess.Popen(cmd, shell=True, bufsize=-1, stdin=last_stdout,
                                stdout=subprocess.PIPE, stderr=stderr)
        last_stdout = proc.stdout
        procs.append(proc)
        
    for cmd in cmds[-1:]:
        # The last command, if any, just dumps to our standard output
        proc = subprocess.Popen(cmd, shell=True, bufsize=-1, stdin=last_stdout,
                                stdout=stdout, stderr=stderr)
        procs.append(proc)
    
    
    # We collect the return codes
    statuses = []
        
    # based on a comment in http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    timer = Timer(timeout_sec, timeout_fail, [procs, cmds])
    try:
        timer.start()
        
        for proc, cmd in itertools.izip(procs, cmds):
            sts = proc.wait()
            statuses.append(sts)
        
            if sts != 0:
                message = "Command: {} in pipeline {} exited with non-zero status {}".format(cmd, " | ".join(cmds), sts)
                if fail_hard is True:
                    raise RuntimeError(message)
                else:
                    RealTimeLogger.get().warning(message) 
        
    finally:
        timer.cancel()

    if len(statuses) > 0:
        # Return the max return code (0 if everything worked)
        return max(statuses)
    else:
        # Nothing bad haoppened because nothing happened
        return 0

def alignment_region_tag(alignment_key):
    """
    Given an alignment key formatted like <region>/<graph>/<sample>.gam,
    produce the region name.
    
    """
    
    region = alignment_key.split("/")[-3]
    # TODO: can we not hardcode a list of known regions?
    assert region in ["brca1", "brca2", "cenx", "lrc_kir", "sma", "mhc"]
    return region

def alignment_graph_tag(alignment_key):
    """
    Given an alignment key formatted like <region>/<graph>/<sample>.gam,
    produce the graph name.
    
    """
    return alignment_key.split("/")[-2]
    
def alignment_sample_tag(alignment_key):
    """
    Given an alignment key formatted like <region>/<graph>/<sample>.gam,
    produce the sample name.
    
    """
    return os.path.splitext(os.path.basename(alignment_key))[0]

def graph_key(alignment_key):
    """
    Get the graph key (i.e. relative path) for the garph athat a GAM file was
    aligned against, given a GAM key (i.e. relative path).
    
    """
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    # Special-case the De Bruijn graphs and stick the -kwhatever at the end.
    if graph.find("debruijn") == 0:
        assert graph.find("-") == 8
        tag = graph[8:]
        graph = "debruijn"
    else:
        tag = ""
    return "{}-{}{}.vg".format(graph, region, tag)    
    
def pileup_key(alignment_key, condition):
    """
    Given a GAM file key, get the key under which the pileup for that GAM file
    should be stored, under the given experimental conditions.
    
    """
    
    name = os.path.splitext(os.path.basename(alignment_key))[0]
    name += ".vgpu"
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    # Stratify by region, then graph, then options used for the filtering and pileup, then sample.
    return "/".join([region, graph, condition.get_pileup_condition_name(), name])
    
def glennfile_key(alignment_key, condition):
    """
    Get the key for the Glennfile in the cache, based on the experimental
    condition and the original GAM name.
    
    """
   
    name = os.path.splitext(os.path.basename(alignment_key))[0]
    name += "_call.tsv"
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    # We nest inside the pileup condition so it's easier to find the relevant
    # pileup.
    return "/".join([region, graph, condition.get_pileup_condition_name(),
        condition.get_glennfile_condition_name(), name])
        
def augmented_graph_key(alignment_key, condition):
    """
    Get the key for the augmented graph in the cache, based on the experimental
    condition and the original GAM name.
    
    """
   
    name = os.path.splitext(os.path.basename(alignment_key))[0]
    name += "_augmented.vg"
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    # We nest inside the pileup condition so it's easier to find the relevant
    # pileup.
    return "/".join([region, graph, condition.get_pileup_condition_name(),
        condition.get_glennfile_condition_name(), name])
        
def vcf_key(alignment_key, condition):
    """
    Get the key for the VCF file in the cache, based on the experimental
    condition and the original GAM name.
    
    """
   
    name = os.path.splitext(os.path.basename(alignment_key))[0]
    name += "_sample.vcf"
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    return "/".join([region, graph, condition.get_pileup_condition_name(),
        condition.get_glennfile_condition_name(),
        condition.get_vcf_condition_name(), name])
        
def vcf_log_key(alignment_key, condition):
    """
    Get the key for the glenn2vcf log file in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = os.path.splitext(os.path.basename(alignment_key))[0]
    name += "_sample.vcf.stderr"
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    return "/".join([region, graph, condition.get_pileup_condition_name(),
        condition.get_glennfile_condition_name(),
        condition.get_vcf_condition_name(), name])

def make_pileup(job, gam_key, condition, options):
    """
    Toil job to make a pileup from the given GAM, in the given experimental
    condition.
    
    Loads the GAM from the input GAM IOStore, and saves the pileup to the right
    place in the cache IOStore for the given sample.
    
    Returns nothing.
    
    """
    
    # Make IOStores
    gam_store = IOStore.get(options.in_gams)
    graph_store = IOStore.get(options.in_graphs)
    cache_store = IOStore.get(options.cache)
    
    # Determine output key
    out_pileup_key = pileup_key(gam_key, condition)
    
    if cache_store.exists(out_pileup_key):
        # We already made this file
        return
    
    # Download the GAM
    input_gam = gam_store.get_input_file(job, gam_key)
    
    # And the graph it was aligned to
    input_graph = graph_store.get_input_file(job, graph_key(gam_key))
    
    # Plan where the output pileup goes
    out_pileup_path = "{}/pileup.vgpu".format(job.fileStore.getLocalTempDir())
    
    # Run the filter and pileup steps, and die if they fail
    pipeline = []
    pipeline.append("vg filter {} {}".format(input_gam,
        condition.get_read_filter_options()))
    pipeline.append("vg pileup {} - {} -t {} > {}".format(input_graph,
        condition.get_pileup_options(), job.cores, out_pileup_path))
    run(pipeline, fail_hard = True)
    
    # Upload the pileup to the cache for the current experimental conditions
    cache_store.write_output_file(out_pileup_path, out_pileup_key)
    
    
def make_glennfile_from_pileup(job, gam_key, condition, options):
    """
    Toil job which, assuming that the GAM has already been turned into a pileup
    in the cache for the given experimental condition, produces the augmented
    graph and associated glennfile in the cache for the current condition.
    
    Returns nothing.
    """
    
    # Make IOStores
    graph_store = IOStore.get(options.in_graphs)
    cache_store = IOStore.get(options.cache)
    
    # Determine output filenames
    out_glennfile_key = glennfile_key(gam_key, condition)
    out_augmented_graph_key = augmented_graph_key(gam_key, condition)
    
    if cache_store.exists(out_glennfile_key) and cache_store.exists(out_augmented_graph_key):
        # We already made these files
        return
    
    # Get the non-augmented graph
    input_graph = graph_store.get_input_file(job, graph_key(gam_key))
    
    # Get the pileup from the cache
    input_pileup = cache_store.get_input_file(job, pileup_key(gam_key, condition))
    
    # Plan where the output glennfile goes
    out_glennfile = "{}/sample.glenn".format(job.fileStore.getLocalTempDir())
    
    # And where the output augmented graph goes
    out_augmented_graph = "{}/sample.vg".format(job.fileStore.getLocalTempDir())
    
    # Do the actual vg call-ing
    pipeline = []
    pipeline.append("vg call {} {} {} -l -c {} -t {} > {}".format(
        input_graph, input_pileup, condition.get_call_options(),
        out_glennfile, job.cores, out_augmented_graph))
    run(pipeline, fail_hard = True)
    
    # Save the glennfile and augmented graph back to the cache
    cache_store.write_output_file(out_glennfile, out_glennfile_key)
    cache_store.write_output_file(out_augmented_graph, out_augmented_graph_key)
    
def make_vcf_from_glennfile(job, gam_key, condition, options):
    """
    Toil job which, assuming that the Glennfile and augmented graph have already
    been made for the given sample and experimental condition, produces the VCF.
    
    Needs the regions option to be specified.
    
    Returns nothing.
    """
    
    # Make IOStores
    cache_store = IOStore.get(options.cache)
    region_store = IOStore.get(options.regions)
    
    # Determine output keys
    out_vcf_key = vcf_key(gam_key, condition)
    out_vcf_log_key = vcf_log_key(gam_key, condition)
    
    if cache_store.exists(out_vcf_key) and cache_store.exists(out_vcf_log_key):
        # We already made these files
        return
    
    # Get the augmented graph from the cache
    input_augmented_graph = cache_store.get_input_file(job, augmented_graph_key(gam_key, condition))
    
    # Get the glennfile from the cache
    input_glennfile = cache_store.get_input_file(job, glennfile_key(gam_key, condition))
    
    # Get the BED that tells us where the region is
    region_bed = region_store.get_input_file(job, alignment_region_tag(gam_key).upper() + ".bed")

    with open(region_bed) as f:
        # Read the contig and offset we want our VCF to be in from the BED file.
        contig, offset = f.readline().split()[0:2]
        
    # Get the sample name
    sample_name = alignment_sample_tag(gam_key)
    
    # Plan where to put the output VCF
    out_vcf = "{}/sample.vcf".format(job.fileStore.getLocalTempDir())
    
    # And the glenn2vcf error log (which has bases dropped, etc.)
    out_errlog = "{}/sample.err".format(job.fileStore.getLocalTempDir())
    
    # Do the actual VCF conversion
    pipeline = []
    pipeline.append("glenn2vcf {} {} -o {} -c {} -s {} {} > {} 2> {}".format(
        input_augmented_graph, input_glennfile, offset, contig, sample_name,
        condition.get_vcf_options(), out_vcf, out_errlog))
    run(pipeline, fail_hard = True)
    
    # Save the vcf and its error log back to the cache
    cache_store.write_output_file(out_vcf, out_vcf_key)
    cache_store.write_output_file(out_errlog, out_vcf_log_key)
    
def run_experiment(job, options):
    """
    Toil job to run an experiment on a variety of conditions and compare the
    results.
    """
    
    # Make the IOStore we can search for GAMs
    gam_store = IOStore.get(options.in_gams)
    
    
    for region_dir in gam_store.list_input_directory(""):
        # Within every region we have samples for, look through all the
        # different graphs.
        for graph_dir in gam_store.list_input_directory(region_dir):
            # Within every graph for a region, we have a collection of samples.
            
            if ("{}:{}".format(region_dir, graph_dir) in options.blacklist or
                region_dir in options.blacklist or
                graph_dir in options.blacklist):
                # We don't want to process this region/graph pair.
                RealTimeLogger.get().info("Skipping {} graph {}".format(
                    region_dir, graph_dir))
                continue
                
            for filename in gam_store.list_input_directory("{}/{}".format(
                region_dir, graph_dir)):
                # Look at each potential sample file
                
                # Is this file a sample?
                match = re.match("(.+)\\.gam$", filename)
                
                if not match:
                    # It's not a sample
                    continue
                    
                # Otherwise, compose the full GAM key
                gam_key = "{}/{}/{}".format(region_dir, graph_dir, filename)
                
                # Make some experimental conditions with filter, pileup, call,
                # and glenn2vcf options. TODO: change to long options for
                # readability
                condition = ExperimentCondition({
                    "-r": 0.90,
                    "-d": 0.05,
                    "-e": 0.05,
                    "-a": "",
                    "-f": "",
                    "-u": "",
                    "-s": 10000,
                    "-o": 10
                }, {
                    "-w": 40,
                    "-m": 10,
                    "-q": 10
                }, {
                    "-r": 0.0001,
                    "-b": 0.4,
                    "-f": 0.25,
                    "-d": 11
                }, {
                    "--depth": 10
                })
                
                # Kick off a pipeline to make the variant calls.
                # TODO: assumes all the extra directories we need to read stuff from are set
                pileup_job = job.addChildJobFn(make_pileup, gam_key, condition, options,
                    cores=1, memory="10G", disk="10G")
                glennfile_job = pileup_job.addFollowOnJobFn(make_glennfile_from_pileup, gam_key, condition, options,
                    cores=1, memory="10G", disk="10G")
                vcf_job = glennfile_job.addFollowOnJobFn(make_vcf_from_glennfile, gam_key, condition, options,
                    cores=1, memory="10G", disk="10G")
    
def main(args):
    
    options = parse_args(args) 
    
    RealTimeLogger.start_master()

    # Make a root job
    root_job = Job.wrapJobFn(run_experiment, options,
                             cores=1, memory="2G", disk="2G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

