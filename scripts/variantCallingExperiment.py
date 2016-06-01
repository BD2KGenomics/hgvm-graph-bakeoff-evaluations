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
import tsv
import gzip
import json
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs, IOStore, de_defaultdict

class ExperimentCondition:
    """
    Represents a combination of read filtering, graph augmenting, vcf
    generating, and vcf filtering options.
    """
    
    def __init__(self, read_filter_options, pileup_options, call_options, vcf_options, vcfeval_options):
        """
        Take dicts from option string to option value for read filtering,
        pileup-ing, graph augmenting, vcf conversion, and vcf evaluation, and
        produces an ExperimentCondition wrapping them. Option values may only be
        single strings, and an option may only occur a single time (since we
        have a dict). TODO: also add VCF filtering
        
        """
        
        self.read_filter_options = read_filter_options
        self.pileup_options = pileup_options
        self.call_options = call_options
        self.vcf_options = vcf_options
        self.vcfeval_options = vcfeval_options
        
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
        return re.sub('[^a-zA-Z0-9_.]', "", re.sub('\s+', "_", re.sub("_", "", string)).strip("_"))
        
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
    
    def get_vcfeval_options(self):
        """
        Return the options string for vcfeval.
        """
        return self.dict_to_string(self.vcfeval_options)
        
    def get_pileup_condition_name(self):
        """
        Return a string that can be part of a filesystem path and which depends
        on all the parameters that affect the pielup.
        """
        
        # Depends on the filter and pileup options
        return self.string_to_path(self.get_read_filter_options()) + "/" + self.string_to_path(self.get_pileup_options())
    
    def get_glennfile_condition_name(self):
        """
        Return a string that can be part of a filesystem path (potentially
        including slashes) and which depends on all the parameters that affect
        the glenn file, including those that affect the pileup.
        
        """
        
        # Depends on the pileup file and the call options
        return self.get_pileup_condition_name() + "/" + self.string_to_path(self.get_call_options())
        
    def get_vcf_condition_name(self):
        """
        Return a string that can be part of a filesystem path (potentially
        including slashes) and which depends on all the parameters that affect
        the final VCF file, given the glenn file.
        
        """
        
        # Depends on the gelnnfile and the glenn2vcf options
        return self.get_glennfile_condition_name() + "/" + self.string_to_path(self.get_vcf_options())
        
    def get_vcfeval_condition_name(self):
        """
        Return a string that can be part of a filesystem path (potentially
        including slashes) and which depends on all the parameters that affect
        the vcfeval evaluation.
        
        """
        
        # Depends on the gelnnfile and the glenn2vcf options
        return self.get_vcf_condition_name() + "/" + self.string_to_path(self.get_vcfeval_options())
        
    def __repr__(self):
        """
        Represent this object as a string for debugging.
        """
        
        # This happens to include all the options
        return self.get_vcfeval_condition_name()
        
        
        
        
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
        help="IOStore to search for truth vcfs for comparison. For each sample "
        "and region we expect <sample>/<REGION>.vcf.gz and "
        "<sample>/<REGION>.vcf.gz.tbi")
    parser.add_argument("--regions",
        help="IOStore to search for region BEDs. For each region we "
        "expect <REGION>.bed, with the region in upper-case.")
    parser.add_argument("--sdf", default="data/g1kvcf/chrom.sdf",
        help="SDF-format directory structure for reference assembly for "
        "vcfeval")
    parser.add_argument("--blacklist", nargs="+", default=["vglr", "mhc:camel",
        "lrc_kir:camel", "debruijn-k31", "debruijn-k63", "cenx", "simons",
        "curoverse", "flat1kg", "primary1kg"],
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
   
def cache_key_stem(alignment_key):
    """
    Get the cache key stem (region/graph) for the given GAM key.
    """
    
    region = alignment_region_tag(alignment_key)
    graph = alignment_graph_tag(alignment_key)
    
    return region + "/" + graph
    
    
def pileup_key(alignment_key, condition):
    """
    Given a GAM file key, get the key under which the pileup for that GAM file
    should be stored, under the given experimental conditions.
    
    """
    
    name = alignment_sample_tag(alignment_key)
    name += ".vgpu"
    return "/".join([cache_key_stem(alignment_key), condition.get_pileup_condition_name(), name])
    
def glennfile_key(alignment_key, condition):
    """
    Get the key for the Glennfile in the cache, based on the experimental
    condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_call.tsv"
    return "/".join([cache_key_stem(alignment_key), condition.get_glennfile_condition_name(), name])
        
def augmented_graph_key(alignment_key, condition):
    """
    Get the key for the augmented graph in the cache, based on the experimental
    condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_augmented.vg"
    return "/".join([cache_key_stem(alignment_key), condition.get_glennfile_condition_name(), name])
        
def vcf_compressed_key(alignment_key, condition):
    """
    Get the key for the compressed VCF file in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_sample.vcf.gz"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcf_condition_name(), name])
        
def vcf_index_key(alignment_key, condition):
    """
    Get the key for the index of the compressed VCF file in the cache, based on
    the experimental condition and the original GAM name.
    
    """
   
    return vcf_compressed_key(alignment_key, condition) + ".tbi"
        
def vcf_log_key(alignment_key, condition):
    """
    Get the key for the glenn2vcf log file in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_sample.vcf.stderr"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcf_condition_name(), name])
        
def vcfeval_summary_key(alignment_key, condition):
    """
    Get the key for the vcfeval comparison summary in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_summary.txt"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcfeval_condition_name(), name])
        
def vcfeval_fp_key(alignment_key, condition):
    """
    Get the key for the vcfeval false positives VCF in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_fp.vcf.gz"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcfeval_condition_name(), name])
        
def vcfeval_fn_key(alignment_key, condition):
    """
    Get the key for the vcfeval false negatives VCF in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_fn.vcf.gz"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcfeval_condition_name(), name])
        
def vcfeval_roc_key(alignment_key, condition):
    """
    Get the key for the vcfeval weighted ROC curve in the cache, based on the
    experimental condition and the original GAM name.
    
    """
   
    name = alignment_sample_tag(alignment_key)
    name += "_weighted_roc.tsv.gz"
    return "/".join([cache_key_stem(alignment_key), condition.get_vcfeval_condition_name(), name])
        
def truth_compressed_key(alignment_key):
    """
    Get the key for the compressed truth VCF for this sample in the **truth
    IOstore**, not the cache.
    
    """
   
    sample = os.path.splitext(os.path.basename(alignment_key))[0]
    region = alignment_region_tag(alignment_key)
    return "{}/{}.vcf.gz".format(sample, region.upper())
    
def truth_index_key(alignment_key):
    """
    Get the key for the index of the compressed truth VCF for this sample in the
    **truth IOstore**, not the cache.
    
    """
   
    return truth_compressed_key(alignment_key) + ".tbi"

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
    out_vcf_compressed_key = vcf_compressed_key(gam_key, condition)
    out_vcf_index_key = vcf_index_key(gam_key, condition)
    out_vcf_log_key = vcf_log_key(gam_key, condition)
    
    if (cache_store.exists(out_vcf_compressed_key) and 
        cache_store.exists(out_vcf_index_key) and
        cache_store.exists(out_vcf_log_key)):
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
    
    # And its compressed and indexed versions
    out_vcf_compressed = out_vcf + ".gz"
    out_vcf_index = out_vcf_compressed + ".tbi"
    
    # Plan where to put the intermediate unsorted VCF
    unsorted_vcf = "{}/unsorted.vcf".format(job.fileStore.getLocalTempDir())
    
    # And the glenn2vcf error log (which has bases dropped, etc.)
    out_errlog = "{}/sample.err".format(job.fileStore.getLocalTempDir())
    
    # Do the actual VCF conversion
    pipeline = []
    pipeline.append("glenn2vcf {} {} -o {} -c {} -s {} {} > {} 2> {}".format(
        input_augmented_graph, input_glennfile, offset, contig, sample_name,
        condition.get_vcf_options(), unsorted_vcf, out_errlog))
    run(pipeline, fail_hard = True)
    
    # Sort the VCF
    pipeline = []
    pipeline.append("scripts/vcfsort {} > {}".format(unsorted_vcf, out_vcf))
    run(pipeline, fail_hard = True)
    
    # Compress and index the VCF
    run(["bgzip {} -c > {}".format(out_vcf, out_vcf_compressed)], fail_hard=True)
    # TODO: This is forced to append .tbi as the index name
    run(["tabix -f -p vcf {}".format(out_vcf_compressed)], fail_hard=True)
    
    # Save the compressed VCF, its index, and its error log back to the cache
    cache_store.write_output_file(out_vcf_compressed, out_vcf_compressed_key)
    cache_store.write_output_file(out_vcf_index, out_vcf_index_key)
    cache_store.write_output_file(out_errlog, out_vcf_log_key)
    
def make_vcfeval_from_vcf(job, gam_key, condition, options):
    """
    Compute the performance of the given GAM (aligned to the given graph)
    against the truth set for its region.
    
    Places the vcfeval summary file for the given sample in the cache.
    """
    
    # Make IOStores
    cache_store = IOStore.get(options.cache)
    truth_store = IOStore.get(options.truth)
    
    # Determine where the resuls go
    # Summary
    out_summary_key = vcfeval_summary_key(gam_key, condition)
    # False positives VCF
    out_fp_key = vcfeval_fp_key(gam_key, condition)
    # False negatives VCF
    out_fn_key = vcfeval_fn_key(gam_key, condition)
    # ROC curve data
    out_roc_key = vcfeval_roc_key(gam_key, condition)
    
    if (cache_store.exists(out_summary_key) and
        cache_store.exists(out_fp_key) and
        cache_store.exists(out_fn_key) and
        cache_store.exists(out_roc_key)):
        # We already did this
        return

    # Get the query VCF
    query_vcf_compressed = cache_store.get_input_file(job, vcf_compressed_key(gam_key, condition))
    query_vcf_index = cache_store.get_input_file(job, vcf_index_key(gam_key, condition))
    
    if query_vcf_index != query_vcf_compressed + ".tbi":
        # Hack them over to the right names with symlinks
        new_vcf_name = "{}/sample.vcf.gz".format(job.fileStore.getLocalTempDir())
        os.symlink(query_vcf_compressed, new_vcf_name)
        os.symlink(query_vcf_index, new_vcf_name + ".tbi")
        query_vcf_compressed = new_vcf_name
        query_vcf_index = query_vcf_compressed + ".tbi"
        
    # Find the truth VCF
    truth_vcf_compressed = truth_store.get_input_file(job, truth_compressed_key(gam_key))
    truth_vcf_index = truth_store.get_input_file(job, truth_index_key(gam_key))
    
    if truth_vcf_index != truth_vcf_compressed + ".tbi":
        # Hack them over to the right names with symlinks
        new_vcf_name = "{}/truth.vcf.gz".format(job.fileStore.getLocalTempDir())
        os.symlink(truth_vcf_compressed, new_vcf_name)
        os.symlink(truth_vcf_index, new_vcf_name + ".tbi")
        truth_vcf_compressed = new_vcf_name
        truth_vcf_index = truth_vcf_compressed + ".tbi"
    
    # Decide on an output directory
    out_dir = "{}/vcfeval".format(job.fileStore.getLocalTempDir())
    
    # Do the actual VCF conversion
    pipeline = []
    pipeline.append("rtg vcfeval -b {} -c {} -t {} -o {} {}".format(
        truth_vcf_compressed, query_vcf_compressed, options.sdf, out_dir,
        condition.get_vcfeval_options()))
    run(pipeline, fail_hard=True)
    
    # Save the result files back to the cache
    cache_store.write_output_file(out_dir + "/summary.txt", out_summary_key)
    cache_store.write_output_file(out_dir + "/fp.vcf.gz", out_fp_key)
    cache_store.write_output_file(out_dir + "/fn.vcf.gz", out_fn_key)
    cache_store.write_output_file(out_dir + "/weighted_roc.tsv.gz", out_roc_key)
    
def get_max_f_score(job, gam_key, condition, options):
    """
    Given the GAM file key for a sample that has already had vcfeval run under
    the given conditions, parse the vcfeval roc and return the biggest F score.
    
    """
    
    # Make the IOStore
    cache_store = IOStore.get(options.cache)
    
    # Find the ROC curve
    roc_key = vcfeval_roc_key(gam_key, condition)
    
    # Get the file
    roc_compressed = cache_store.get_input_file(job, roc_key)
    
    # Read it
    reader = tsv.TsvReader(gzip.GzipFile(roc_compressed))
    
    # What's the max F score we found?
    max_f_score = None
    for parts in reader:
        # Parse all the F scores
        f_score = float(parts[6])
        
        if max_f_score is None or f_score > max_f_score:
            # And keep the max
            max_f_score = f_score
            
    # Return the max F score.
    return max_f_score
    
def run_conditions(job, gam_key, conditions, options):
    """
    Run the pipeline for all the conditions given, and put all the results into
    the cache.
    
    Returns a dict from condition to best F score.
    
    Avoids re-doing work in parallel.
    """
    
    # Maps from pileup condition name to pileup job
    pileup_jobs = {}
    # And glennfile condition name to glennfile job
    glennfile_jobs = {}
    # And vcf condition name to glenn2vcf job
    vcf_jobs = {}
    # And vcfeval condition name to vcfeval job
    vcfeval_jobs = {}
    
    # And vcfeval condition name to F score (promise)
    f_score_jobs = {}
    
    # This is the dict from condition to best F score.
    f_scores = {}
    
    for condition in conditions:
        # For each condition, what pileup do we need?
        pileup = condition.get_pileup_condition_name()
        if not pileup_jobs.has_key(pileup):
            # We need to make this pileup
            pileup_jobs[pileup] = job.addChildJobFn(make_pileup, gam_key,
                condition, options, cores=1, memory="10G", disk="10G")
                
        # And what glennfile do we need?
        glennfile = condition.get_glennfile_condition_name()
        if not glennfile_jobs.has_key(glennfile):
            # We need to make this glennfile
            glennfile_jobs[glennfile] = pileup_jobs[pileup].addFollowOnJobFn(
                make_glennfile_from_pileup, gam_key, condition, options,
                cores=1, memory="10G", disk="10G")
                
        # And what vcf do we need?
        vcf = condition.get_vcf_condition_name()
        if not vcf_jobs.has_key(vcf):
            # We need to make this vcf
            vcf_jobs[vcf] = glennfile_jobs[glennfile].addFollowOnJobFn(
                make_vcf_from_glennfile, gam_key, condition, options,
                cores=1, memory="10G", disk="10G")
                
        # And what vcfeval do we need?
        vcfeval = condition.get_vcfeval_condition_name()
        if not vcfeval_jobs.has_key(vcfeval):
            # We need to make this vcfeval result
            vcfeval_jobs[vcfeval] = vcf_jobs[vcf].addFollowOnJobFn(
                make_vcfeval_from_vcf, gam_key, condition, options,
                cores=1, memory="10G", disk="10G")
                
        # And what job do we need to compute the F score?
        if not f_score_jobs.has_key(vcfeval):
            # We need to parse the results and extract the max F score
            f_score_jobs[vcfeval] = vcfeval_jobs[vcfeval].addFollowOnJobFn(
                get_max_f_score, gam_key, condition, options,
                cores=1, memory="2G", disk="2G")
                
            # And save the result in our return dict for this condition
            f_scores[condition] = f_score_jobs[vcfeval].rv()
            
    # Send back f scores by condition.
    return f_scores
                
                
def run_experiment(job, options):
    """
    Toil job to run an experiment on a variety of conditions and compare the
    results.
    """
    
    # Make the IOStore we can search for GAMs
    gam_store = IOStore.get(options.in_gams)
    
    # This will hold best F score by region, graph, sample, and then condition.
    # We stick in dicts by condition.
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    
    
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
                condition = ExperimentCondition(
                    { # vg filter 
                        #"-r": 0.90,
                        "-d": 0.05,
                        "-e": 0.05,
                        "-a": "",
                        "-f": "",
                        "-u": "",
                        "-s": 10000,
                        "-o": 10
                    }, { # vg pileup
                        "-w": 40,
                        "-m": 10,
                        "-q": 10
                    }, { # vg call
                        "-r": 0.0001, # Prior for being heterozygous
                        "-b": 1, # Max strand bias
                        "-f": 0.05, # Min fraction of reads required to support a variant
                        "-d": 2 # Min pileup depth
                    }, { # glenn2vcf
                        "--depth": 10
                    }, { # vcfeval
                        "--all-records": "",
                        "--vcf-score-field": "DP",
                        "--squash-ploidy": ""
                    })
                
                # Kick off a pipeline to make the variant calls.
                # TODO: assumes all the extra directories we need to read stuff from are set
                exp_job = job.addChildJobFn(run_conditions, gam_key, [condition], options,
                    cores=1, memory="2G", disk="10G")
                    
                # Save the best F score by condition under this region, graph, and sample filename
                results[region_dir][graph_dir][filename] = exp_job.rv()
                
    # Give back the results
    # TODO: we run it through JSON to fix the pickle-ability.
    return de_defaultdict(results)
    
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
    
    print("Root return value:")
    print(root_job.rv())
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

