#!/usr/bin/env python2.7
"""
parallelMappingEvaluation.py: Run the mapping evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import fnmatch

import dateutil.parser

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
    parser.add_argument("server_list", type=argparse.FileType("r"),
        help="TSV file continaing <region>\t<url> lines for servers to test")
    parser.add_argument("sample_store", type=IOStore.absolute,
        help="sample input IOStore with <region>/<sample>/<sample>.bam.fq")
    parser.add_argument("out_store", type=IOStore.absolute,
        help="output IOStore to create and fill with alignments and stats")
    parser.add_argument("--server_version", default="v0.6.g",
        help="server version to add to URLs")
    parser.add_argument("--sample_pattern", default="*", 
        help="fnmatch-style pattern for sample names")
    parser.add_argument("--sample_limit", type=int, default=float("inf"), 
        help="number of samples to use")
    parser.add_argument("--edge_max", type=int, default=0, 
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=10, 
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--bin_url",
        default="https://hgvm.blob.core.windows.net/hgvm-bin",
        help="URL to download sg2vg and vg binaries from, without Docker")
    parser.add_argument("--use_path_binaries", action="store_true",
        help="use system vg and sg2vg instead of downloading them")
    parser.add_argument("--overwrite", default=False, action="store_true",
        help="overwrite existing result files")
    parser.add_argument("--restat", default=False, action="store_true",
        help="recompute and overwrite existing stats files")
    parser.add_argument("--reindex", default=False, action="store_true",
        help="don't re-use existing indexed graphs")
    parser.add_argument("--too_old", default=None, type=str,
        help="recompute stats files older than this date")
    parser.add_argument("--index_mode", choices=["rocksdb", "gcsa-kmer",
        "gcsa-mem"], default="gcsa-mem",
        help="type of vg index to use for mapping")
    parser.add_argument("--include_pruned", action="store_true",
        help="use the pruned graph in the index")
    parser.add_argument("--include_primary", action="store_true",
        help="use the primary path in the index")
    
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)
    
    
# Reverse complement needs a global translation table
reverse_complement_translation_table = string.maketrans("ACGTN", "TGCAN")
def reverse_complement(sequence):
    """
    Compute the reverse complement of a DNA sequence.
    
    Follows algorithm from <http://stackoverflow.com/a/26615937>
    """
    
    if isinstance(sequence, unicode):
        # Encode the sequence in ASCII for easy translation
        sequence = sequence.encode("ascii", "replace")
    
    # Translate and then reverse
    return sequence.translate(reverse_complement_translation_table)[::-1]
    
def count_Ns(sequence):
    """
    Return the number of N bases in the given DNA sequence
    """
    
    n_count = 0
    for item in sequence:
        if item == "N":
            n_count += 1
            
    return n_count

def run_all_alignments(job, options):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    if options.use_path_binaries:
        # We don't download any bianries and don't maintain a bin_dir
        bin_dir_id = None
    else:
        # Retrieve binaries we need
        RealTimeLogger.get().info("Retrieving binaries from {}".format(
            options.bin_url))
        bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
        robust_makedirs(bin_dir)
        subprocess.check_call(["wget", "{}/sg2vg".format(options.bin_url),
            "-O", "{}/sg2vg".format(bin_dir)])
        subprocess.check_call(["wget", "{}/vg".format(options.bin_url),
            "-O", "{}/vg".format(bin_dir)])
            
        # Make them executable
        os.chmod("{}/sg2vg".format(bin_dir), 0o744)
        os.chmod("{}/vg".format(bin_dir), 0o744)
        
        # Upload the bin directory to the file store
        bin_dir_id = write_global_directory(job.fileStore, bin_dir,
            cleanup=True)
    
    # Make sure we skip the header
    is_first = True
    
    for line in options.server_list:
        if is_first:
            # This is the header, skip it.
            is_first = False
            continue
        
        # We need to read each non-header line
        
        # Break it into its fields
        parts = line.split("\t")
        
        if parts[0].startswith("#"):
            # Skip comments
            continue
            
        if parts[0].startswith("\n"):
            # Skip newlines
            continue
            
        # Pull out the first 3 fields
        region, url, generator = parts[0:3]
        
        # We cleverly just split the lines out to different nodes
        job.addChildJobFn(run_region_alignments, options, bin_dir_id, region,
            url, cores=16, memory="100G", disk="50G")
            
        # Say what we did
        RealTimeLogger.get().info("Running child for {}".format(parts[1]))
        

def run_region_alignments(job, options, bin_dir_id, region, url):
    """
    For the given region, download, index, and then align to the given graph.
    
    """
    
    RealTimeLogger.get().info("Running on {} for {}".format(url, region))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    if bin_dir_id is not None:
        # Download the binaries
        bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
        read_global_directory(job.fileStore, bin_dir_id, bin_dir)
        # We define a string we can just tack onto the binary name and get
        # either the system or the downloaded version.
        bin_prefix = bin_dir + "/"
    else:
        bin_prefix = ""
    
    # Parse the graph URL. It may be either an http(s) URL to a GA4GH server, or
    # a direct file: URL to a vg graph.
    url_parts = urlparse.urlparse(url, "file")
    
    # Either way, it has to include the graph base name, which we use to
    # identify the graph, and which has to contain the region name (like brca2)
    # and the graph type name (like cactus). For a server it's the last folder
    # in the YURL, with a trailing slash, and for a file it's the name of the
    # .vg file.
    basename = re.match('.*/(.*)(/|\.vg)$', url_parts.path).group(1)
        
    # Get graph name (without region and its associated dash) from basename
    graph_name = basename.replace("-{}".format(region), "").replace(
        "{}-".format(region), "")
    
    # Where do we look for samples for this region in the input?
    region_dir = region.upper()
    
    # What samples do we do? List input sample names up to the given limit.
    input_samples = [n for n in sample_store.list_input_directory(region_dir) \
        if fnmatch.fnmatchcase(n, options.sample_pattern)]
    if len(input_samples) > options.sample_limit:
        input_samples = input_samples[:options.sample_limit]
    
    # Work out the directory for the alignments to be dumped in in the output
    alignment_dir = "alignments/{}/{}".format(region, graph_name)
    
    # Also for statistics
    stats_dir = "stats/{}/{}".format(region, graph_name)
    
    # What smaples have been completed? Map from ID to mtime
    completed_samples = {}
    for filename, mtime in out_store.list_input_directory(stats_dir,
        with_times=True):
        # See if every file is a stats file
        match = re.match("(.*)\.json$", filename)
    
        if match and (options.too_old is None or mtime > options.too_old):
            # We found a sample's stats file, and it's new enough.
            completed_samples[match.group(1)] = mtime
        
        if match and (options.too_old is not None and mtime < options.too_old):
            # Say we hit an mtime thing
            RealTimeLogger.get().info("Need to re-run {} because "
                "{} < {}".format(match.group(1), mtime.ctime(),
                options.too_old.ctime()))
            
    RealTimeLogger.get().info("Already have {}/{} completed samples for {} in "
        "{}".format(len(completed_samples), len(input_samples), basename,
        stats_dir))
    
    # What samples haven't been done yet and need doing
    samples_to_run = []
    
    for sample in input_samples:
        # Split out over each sample
        
        if ((not options.overwrite) and (not options.restat) and 
            completed_samples.has_key(sample)):
            # This is already done.
            RealTimeLogger.get().info("Skipping completed alignment of "
                "{} to {} {}".format(sample, graph_name, region))
                
            continue
        else:
            # We need to run this sample
            samples_to_run.append(sample)
            
    if len(samples_to_run) == 0 and not options.reindex:
        # Don't bother indexing the graph if all the samples are done, and we
        # didn't explicitly ask to do it.
        RealTimeLogger.get().info("Nothing to align to {}".format(basename))
        return
    
    # Where will the indexed graph go in the output
    index_key = "indexes/{}-{}-{}/{}/{}.tar.gz".format(options.index_mode,
        options.kmer_size, options.edge_max, region, graph_name)
    
    if (not options.reindex) and out_store.exists(index_key):
        # See if we have an index already available in the output store from a
        # previous run
        
        RealTimeLogger.get().info("Retrieving indexed {} graph from output "
            "store".format(basename))
            
        # Download the pre-made index directory
        tgz_file = "{}/index.tar.gz".format(job.fileStore.getLocalTempDir())
        out_store.read_input_file(index_key, tgz_file)
        
        # Save it to the global file store and keep around the ID.
        # Will be compatible with read_global_directory
        index_dir_id = job.fileStore.writeGlobalFile(tgz_file, cleanup=True)
        
        RealTimeLogger.get().info("Index for {} retrieved "
            "successfully".format(basename))
        
    else:
        # Download the graph, build the index, and store it in the output store
    
        # Work out where the graph goes
        # it will be graph.vg in here
        graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
        robust_makedirs(graph_dir)
        
        graph_filename = "{}/graph.vg".format(graph_dir)
        
        # Download and fix up the graph with this ugly subprocess pipeline
        # sg2vg "${URL}" -u | vg view -Jv - | vg mod -X 100 - | 
        # vg ids -s - > "graphs/${BASENAME}.vg"
        
        with open(graph_filename, "w") as output_file:
        
            # Hold all the popen objects we need for this
            tasks = []
            
            if url_parts.scheme == "file":
                # Grab the vg graph from a local file
                
                RealTimeLogger.get().info("Reading {} to {}".format(
                    url, graph_filename))
                    
                # Just cat the file. We need a process so we can do the
                # tasks[-1].stdout trick.
                tasks.append(subprocess.Popen(["cat", url_parts.path],
                    stdout=subprocess.PIPE))
            
            elif url.endswith(".vg"):
                # Assume it's a vg file
                
                RealTimeLogger.get().info("Downloading {} to {}".format(
                    url, graph_filename))
                    
                tasks.append(subprocess.Popen(["curl", url],
                    stdout=subprocess.PIPE))
                
            else:
                # Assume it's on a server
                
                # Make the real URL with the version
                versioned_url = url + options.server_version
                
                RealTimeLogger.get().info("Downloading {} to {}".format(
                    versioned_url, graph_filename))
                
                # Do the download
                tasks.append(subprocess.Popen(["{}sg2vg".format(bin_prefix),
                    versioned_url, "-u"], stdout=subprocess.PIPE))
                
                # Convert to vg
                tasks.append(subprocess.Popen(["{}vg".format(bin_prefix),
                    "view", "-Jv", "-"], stdin=tasks[-1].stdout,
                    stdout=subprocess.PIPE))
            
            # And cut nodes
            tasks.append(subprocess.Popen(["{}vg".format(bin_prefix), "mod",
                "-X100", "-"], stdin=tasks[-1].stdout, stdout=subprocess.PIPE))
                
            # And sort ids
            tasks.append(subprocess.Popen(["{}vg".format(bin_prefix), "ids",
                "-s", "-"], stdin=tasks[-1].stdout, stdout=output_file))
                
            # Did we make it through all the tasks OK?
            for task in tasks:
                if task.wait() != 0:
                    raise RuntimeError("Pipeline step returned {}".format(
                        task.returncode))
                     
        # TODO: We sometimes don't see the files written immediately, for some
        # reason. Maybe because we opened them? Anyway, this is a hack to wait
        # for them to be on disk and readable.
        time.sleep(1)
        
        
        
        # Now run the indexer.
        # TODO: support both indexing modes
        RealTimeLogger.get().info("Indexing {}".format(graph_filename))
        
        if options.index_mode == "rocksdb":
            # Make the RocksDB index
            subprocess.check_call(["{}vg".format(bin_prefix), "index", "-s", "-k",
                str(options.kmer_size), "-e", str(options.edge_max),
                "-t", str(job.cores), graph_filename, "-d",
                graph_filename + ".index"])
                
        elif (options.index_mode == "gcsa-kmer" or
            options.index_mode == "gcsa-mem"):
            # We want a GCSA2/xg index. We have to prune the graph ourselves.
            # See <https://github.com/vgteam/vg/issues/286>.
            
            # What will we use as our temp combined graph file (containing only
            # the bits of the graph we want to index, used for deduplication)?
            to_index_filename = "{}/to_index.vg".format(
                job.fileStore.getLocalTempDir())
            
            # Where will we save the kmers?
            kmers_filename = "{}/index.graph".format(
                job.fileStore.getLocalTempDir())
                
            with open(to_index_filename, "w") as to_index_file:
                
                if options.include_pruned:
                
                    RealTimeLogger.get().info("Pruning {} to {}".format(
                        graph_filename, to_index_filename))
                    
                    # Prune out hard bits of the graph
                    tasks = []
                    
                    # Prune out complex regions
                    tasks.append(subprocess.Popen(["{}vg".format(bin_prefix), "mod",
                        "-p", "-l", str(options.kmer_size), "-t", str(job.cores),
                        "-e", str(options.edge_max), graph_filename],
                        stdout=subprocess.PIPE))
                        
                    # Throw out short disconnected chunks
                    tasks.append(subprocess.Popen(["{}vg".format(bin_prefix), "mod",
                        "-S", "-l", str(options.kmer_size * 2),
                        "-t", str(job.cores), "-"], stdin=tasks[-1].stdout,
                        stdout=to_index_file))
                        
                    # Did we make it through all the tasks OK?
                    for task in tasks:
                        if task.wait() != 0:
                            raise RuntimeError("Pipeline step returned {}".format(
                                task.returncode))
                                
                    time.sleep(1)
                    
                if options.include_primary:
                
                    # Then append in the primary path. Since we don't knoiw what
                    # "it's called, we retain "ref" and all the 19", "6", etc paths
                    # "from 1KG.
                    
                    RealTimeLogger.get().info(
                        "Adding primary path to {}".format(to_index_filename))
                    
                    # See
                    # https://github.com/vgteam/vg/issues/318#issuecomment-215102199
                    
                    # Generate all the paths names we might have for primary paths.
                    # It should be "ref" but some graphs don't listen
                    ref_names = (["ref", "x", "X", "y", "Y", "m", "M"] +
                        [str(x) for x in xrange(1, 23)])
                        
                    ref_options = []
                    for name in ref_names:
                        # Put each in a -r option to retain the path
                        ref_options.append("-r")
                        ref_options.append(name)

                    tasks = []

                    # Retain only the specified paths (only one should really exist)
                    tasks.append(subprocess.Popen(
                        ["{}vg".format(bin_prefix), "mod", "-N"] + ref_options + 
                        ["-t", str(job.cores), graph_filename], 
                        stdout=to_index_file))
                        
                    # TODO: if we merged the primary path back on itself, it's
                    # possible for it to braid with itself. Right now we just ignore
                    # this and let those graphs take a super long time to index.
                        
                    # Wait for this second pipeline. We don't parallelize with the
                    # first one so we don't need to use an extra cat step.
                    for task in tasks:
                        if task.wait() != 0:
                            raise RuntimeError("Pipeline step returned {}".format(
                                task.returncode))
                             
                    # Wait to make sure no weird file-not-being-full bugs happen
                    # TODO: how do I wait on child process output?
                    time.sleep(1)
                
            time.sleep(1)
                
            # Now we have the combined to-index graph in one vg file. We'll load
            # it (which deduplicates nodes/edges) and then find kmers.
                
            # Save the problematic file
            out_store.write_output_file(to_index_filename,
                "debug/{}-{}-{}-{}-{}.vg".format(options.index_mode,
                options.kmer_size, options.edge_max, region, graph_name))
                
            with open(kmers_filename, "w") as kmers_file:
            
                tasks = []
                
                RealTimeLogger.get().info("Finding kmers in {} to {}".format(
                    to_index_filename, kmers_filename))
                
                # Deduplicate the graph
                tasks.append(subprocess.Popen(["{}vg".format(bin_prefix),
                    "view", "-v", to_index_filename],
                    stdout=subprocess.PIPE))
                
                # Make the GCSA2 kmers file
                tasks.append(subprocess.Popen(["{}vg".format(bin_prefix),
                    "kmers", "-g", "-B", "-k", str(options.kmer_size),
                    "-H", "1000000000", "-T", "1000000001",
                    "-t", str(job.cores), "-"], stdin=tasks[-1].stdout,
                    stdout=kmers_file))
                    
                # Did we make it through all the tasks OK?
                task_number = 0
                for task in tasks:
                    if task.wait() != 0:
                        raise RuntimeError(
                            "Pipeline step {} returned {}".format(
                            task_number, task.returncode))
                    task_number += 1
                            
                # Wait to make sure no weird file-not-being-full bugs happen
                # TODO: how do I wait on child process output?
                time.sleep(1)
                
            time.sleep(1)
                            
            # Where do we put the GCSA2 index?
            gcsa_filename = graph_filename + ".gcsa"
            
            RealTimeLogger.get().info("GCSA-indexing {} to {}".format(
                    kmers_filename, gcsa_filename))
            
            # Make the gcsa2 index. Make sure to use 3 doubling steps to work
            # around <https://github.com/vgteam/vg/issues/301>
            subprocess.check_call(["{}vg".format(bin_prefix), "index", "-t",
                str(job.cores), "-i", kmers_filename, "-g", gcsa_filename,
                "-X", "3"])
                
            # Where do we put the XG index?
            xg_filename = graph_filename + ".xg"
            
            RealTimeLogger.get().info("XG-indexing {} to {}".format(
                    graph_filename, xg_filename))
                    
            subprocess.check_call(["{}vg".format(bin_prefix), "index", "-t",
                str(job.cores), "-x", xg_filename, graph_filename])
        
        else:
            raise RuntimeError("Invalid indexing mode: " + options.index_mode)
            
        # Define a file to keep the compressed index in, so we can send it to
        # the output store.
        index_dir_tgz = "{}/index.tar.gz".format(
            job.fileStore.getLocalTempDir())
            
        # Now save the indexed graph directory to the file store. It can be
        # cleaned up since only our children use it.
        RealTimeLogger.get().info("Compressing index of {}".format(
            graph_filename))
        index_dir_id = write_global_directory(job.fileStore, graph_dir,
            cleanup=True, tee=index_dir_tgz)
            
        # Save it as output
        RealTimeLogger.get().info("Uploading index of {}".format(
            graph_filename))
        out_store.write_output_file(index_dir_tgz, index_key)
        RealTimeLogger.get().info("Index {} uploaded successfully".format(
            index_key))
            
    RealTimeLogger.get().info("Queueing alignment of {} samples to "
        "{} {}".format(len(samples_to_run), graph_name, region))
            
    job.addChildJobFn(recursively_run_samples, options, bin_dir_id, 
        graph_name, region, index_dir_id, samples_to_run,
        cores=1, memory="4G", disk="4G")
            
    RealTimeLogger.get().info("Done making children for {}".format(basename))
   
def recursively_run_samples(job, options, bin_dir_id, graph_name, region,
    index_dir_id, samples_to_run, num_per_call=10):
    """
    Create child jobs to run a few samples from the samples_to_run list, and a
    recursive child job to create a few more.
    
    This is a hack to deal with the problems produced by having a job with
    thousands of children on the Azure job store: the job graph gets cut up into
    tiny chunks of data and stored as table values, and when you have many table
    store operations one of them is likely to fail and screw up your whole
    serialization process.
    
    We have some logic here to decide how much of the sample needs to be rerun.
    If we get a sample, all we know is that it doesn't have an up to date stats
    file, but it may or may not have an alignment file already.
    
    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Get some samples to run
    samples_to_run_now = samples_to_run[:num_per_call]
    samples_to_run_later = samples_to_run[num_per_call:]
    
    # Work out where samples for this region live
    region_dir = region.upper()    
    
    # Work out the directory for the alignments to be dumped in in the output
    alignment_dir = "alignments/{}/{}".format(region, graph_name)
    
    # Also for statistics
    stats_dir = "stats/{}/{}".format(region, graph_name)
    
    for sample in samples_to_run_now:
        # Split out over each sample that needs to be run
        
        # For each sample, know the FQ name
        sample_fastq = "{}/{}/{}.bam.fq".format(region_dir, sample, sample)
        
        # And know where we're going to put the output
        alignment_file_key = "{}/{}.gam".format(alignment_dir, sample)
        stats_file_key = "{}/{}.json".format(stats_dir, sample)
        
        if (not options.overwrite and out_store.exists(alignment_file_key)):
            # If the alignment exists and is OK, we can skip the alignment
            
            mtime = out_store.get_mtime(stats_file_key)
            if mtime is None or mtime < options.too_old or options.restat:
                # All we need to do for this sample is run stats
                RealTimeLogger.get().info("Queueing stat recalculation"
                    " of {} on {} {}".format(sample, graph_name, region))
                
                job.addFollowOnJobFn(run_stats, options, bin_dir_id,
                    index_dir_id, alignment_file_key, stats_file_key,
                    run_time=None, cores=2, memory="4G", disk="10G")
                    
            else:
                # The stats are up to date and the alignment doesn't need
                # rerunning. This shouldn't happen because this sample shouldn't
                # be on the todo list. But it means we can just skip the sample.
                RealTimeLogger.get().warning("SKIPPING sample "
                    "{} on {} {}".format(sample, graph_name, region))
                    
        else:
        
            # Otherwise we need to do an alignment, and then stats.
            
            RealTimeLogger.get().info("Queueing alignment"
                " of {} to {} {}".format(sample, graph_name, region))
            
            # Go and bang that input fastq against the correct indexed graph.
            # Its output will go to the right place in the output store.
            job.addChildJobFn(run_alignment, options, bin_dir_id, sample,
                graph_name, region, index_dir_id, sample_fastq,
                alignment_file_key, stats_file_key, 
                cores=16, memory="100G", disk="50G")
            
    if len(samples_to_run_later) > 0:
        # We need to recurse and run more later.
        RealTimeLogger.get().debug("Postponing queueing {} samples".format(
            len(samples_to_run_later)))
            
        if len(samples_to_run_later) < num_per_call:
            # Just run them all in one batch
            job.addChildJobFn(recursively_run_samples, options, bin_dir_id,
                    graph_name, region, index_dir_id, samples_to_run_later,
                    num_per_call, cores=1, memory="4G", disk="4G")
        else:
            # Split them up
        
            part_size = len(samples_to_run_later) / num_per_call
            
            RealTimeLogger.get().info("Splitting remainder of {} {} into {} "
                "parts of {}".format(graph_name, region, num_per_call,
                part_size))
            
            for i in xrange(num_per_call + 1):
                # Do 1 more part for any remainder
                
                # Grab this bit of the rest
                part = samples_to_run_later[(i * part_size) :
                    ((i + 1) * part_size)]
                
                if len(part) > 0:
                
                    # Make a job to run it
                    job.addChildJobFn(recursively_run_samples, options,
                        bin_dir_id, graph_name, region, index_dir_id, part,
                        num_per_call, cores=1, memory="4G", disk="4G")
        
        
    
            
def save_indexed_graph(job, options, index_dir_id, output_key):
    """
    Save the index dir tar file in the given output key.
    
    Runs as a child to ensure that the global file store can actually
    produce the file when asked (because within the same job, depending on Toil
    guarantees, it might still be uploading).
    
    """
    
    RealTimeLogger.get().info("Uploading {} to output store...".format(
        output_key))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # Get the tar.gz file
    RealTimeLogger.get().info("Downloading global file {}".format(index_dir_id))
    local_path = job.fileStore.readGlobalFile(index_dir_id)
    
    size = os.path.getsize(local_path)
    
    RealTimeLogger.get().info("Global file {} ({} bytes) for {} read".format(
        index_dir_id, size, output_key))
    
    # Save it as output
    out_store.write_output_file(local_path, output_key)
    
    RealTimeLogger.get().info("Index {} uploaded successfully".format(
        output_key))
    
   
def run_alignment(job, options, bin_dir_id, sample, graph_name, region,
    index_dir_id, sample_fastq_key, alignment_file_key, stats_file_key):
    """
    Align the the given fastq from the input store against the given indexed
    graph (in the file store as a directory) and put the GAM and statistics in
    the given output keys in the output store.
    
    Assumes that the alignment actually needs to be redone.
    
    """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    # How long did the alignment take to run, in seconds?
    run_time = None
    
    if bin_dir_id is not None:
        # Download the binaries
        bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
        read_global_directory(job.fileStore, bin_dir_id, bin_dir)
        # We define a string we can just tack onto the binary name and get
        # either the system or the downloaded version.
        bin_prefix = bin_dir + "/"
    else:
        bin_prefix = ""
    
    # Download the indexed graph to a directory we can use
    graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
    read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)
    
    # Also we need the sample fastq
    fastq_file = "{}/input.fq".format(job.fileStore.getLocalTempDir())
    sample_store.read_input_file(sample_fastq_key, fastq_file)
    
    # And a temp file for our aligner output
    output_file = "{}/output.gam".format(job.fileStore.getLocalTempDir())
    
    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:
    
        # Start the aligner and have it write to the file
        
        # Plan out what to run
        vg_parts = ["{}vg".format(bin_prefix), "map", "-f", fastq_file,
            "-i", "-M2", "-W", "500", "-u", "0", "-U", "-t", str(job.cores), graph_file]
            
        if options.index_mode == "rocksdb":
            vg_parts += ["-d", graph_file + ".index", "-n3", "-k",
                str(options.kmer_size)]
        elif options.index_mode == "gcsa-kmer":
            # Use the new default context size in this case
            vg_parts += ["-x", graph_file + ".xg", "-g", graph_file + ".gcsa",
                "-n5", "-k", str(options.kmer_size)]
        elif options.index_mode == "gcsa-mem":
            # Don't pass the kmer size, so MEM matching is used
            vg_parts += ["-x", graph_file + ".xg", "-g", graph_file + ".gcsa",
                "-n5"]
        else:
            raise RuntimeError("invalid indexing mode: " + options.index_mode)
        
        RealTimeLogger.get().info(
            "Running VG for {} against {} {}: {}".format(sample, graph_name,
            region, " ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        process = subprocess.Popen(vg_parts, stdout=alignment_file)
            
        if process.wait() != 0:
            # Complain if vg dies
            raise RuntimeError("vg died with error {}".format(
                process.returncode))
                
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        
                
    RealTimeLogger.get().info("Aligned {}".format(output_file))
    
    # Upload the alignment
    out_store.write_output_file(output_file, alignment_file_key)
    
    RealTimeLogger.get().info("Need to recompute stats for new "
        "alignment: {}".format(stats_file_key))

    # Add a follow-on to calculate stats. It only needs 2 cores since it's
    # not really prarllel.
    job.addFollowOnJobFn(run_stats, options, bin_dir_id, index_dir_id,
        alignment_file_key, stats_file_key, run_time=run_time,
        cores=2, memory="4G", disk="10G")
            
      
def run_stats(job, options, bin_dir_id, index_dir_id, alignment_file_key,
    stats_file_key, run_time=None):
    """
    If the stats aren't done, or if they need to be re-done, retrieve the
    alignment file from the output store under alignment_file_key and compute the
    stats file, saving it under stats_file_key.
    
    Uses index_dir_id to get the graph, and thus the reference sequence that
    each read is aligned against, for the purpose of discounting Ns.
    
    Can take a run time to put in the stats.

    Assumes that stats actually do need to be computed, and overwrites any old
    stats.

    TODO: go through the proper file store (and cache) for getting alignment
    data.
    
    """
          
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    sample_store = IOStore.get(options.sample_store)
    out_store = IOStore.get(options.out_store)
    
    RealTimeLogger.get().info("Computing stats for {}".format(stats_file_key))
    
    if bin_dir_id is not None:
        # Download the binaries
        bin_dir = "{}/bin".format(job.fileStore.getLocalTempDir())
        read_global_directory(job.fileStore, bin_dir_id, bin_dir)
        # We define a string we can just tack onto the binary name and get either
        # the system or the downloaded version.
        bin_prefix = bin_dir + "/"
    else:
        bin_prefix = ""
        
    # Download the indexed graph to a directory we can use
    graph_dir = "{}/graph".format(job.fileStore.getLocalTempDir())
    read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)
    
    # Load the node sequences into memory. This holds node sequence string by
    # ID.
    node_sequences = {}
    
    # Read the alignments in in JSON-line format
    read_graph = subprocess.Popen(["{}vg".format(bin_prefix), "view", "-j",
        graph_file], stdout=subprocess.PIPE)
        
    for line in read_graph.stdout:
        # Parse the graph chunk JSON
        graph_chunk = json.loads(line)
        
        for node_dict in graph_chunk.get("node", []):
            # For each node, store its sequence under its id. We want to crash
            # if a node exists for which one or the other isn't defined.
            node_sequences[node_dict["id"]] = node_dict["sequence"]
        
    if read_graph.wait() != 0:
        # Complain if vg dies
        raise RuntimeError("vg died with error {}".format(
            read_graph.returncode))
 
    # Declare local files for everything
    stats_file = "{}/stats.json".format(job.fileStore.getLocalTempDir())
    alignment_file = "{}/output.gam".format(job.fileStore.getLocalTempDir())
    
    # Download the alignment
    out_store.read_input_file(alignment_file_key, alignment_file)
           
    # Read the alignments in in JSON-line format
    read_alignment = subprocess.Popen(["{}vg".format(bin_prefix), "view", "-aj",
        alignment_file], stdout=subprocess.PIPE)
       
    # Count up the stats
    stats = {
        "total_reads": 0,
        "total_mapped": 0,
        "total_multimapped": 0,
        "mapped_lengths": collections.Counter(),
        "unmapped_lengths": collections.Counter(),
        "aligned_lengths": collections.Counter(),
        "primary_scores": collections.Counter(),
        "primary_mapqs": collections.Counter(),
        "primary_identities": collections.Counter(),
        "primary_mismatches": collections.Counter(),
        "primary_indels": collections.Counter(),
        "primary_substitutions": collections.Counter(),
        "secondary_scores": collections.Counter(),
        "secondary_mapqs": collections.Counter(),
        "secondary_identities": collections.Counter(),
        "secondary_mismatches": collections.Counter(),
        "secondary_indels": collections.Counter(),
        "secondary_substitutions": collections.Counter(),
        "run_time": run_time
    }
        
    last_alignment = None
        
    for line in read_alignment.stdout:
        # Parse the alignment JSON
        alignment = json.loads(line)
        
        # How long is this read?
        length = len(alignment["sequence"])
        
        if alignment.has_key("score"):
            # This alignment is aligned.
            # Grab its score
            score = alignment["score"]
        
            # Get the mappings
            mappings = alignment.get("path", {}).get("mapping", [])
        
            # Calculate the exact match bases
            matches = 0
            
            # And total up the instances of indels (only counting those where
            # the reference has no Ns, and which aren't leading or trailing soft
            # clips)
            indels = 0
            
            # And total up the number of substitutions (mismatching/alternative
            # bases in edits with equal lengths where the reference has no Ns).
            substitutions = 0
            
            # What should the denominator for substitution rate be for this
            # read? How many bases are in the read and aligned?
            aligned_length = 0
            
            # What's the mapping quality? May not be defined on some reads.
            mapq = alignment.get("mapping_quality", 0.0)
            
            # And the identity?
            identity = alignment["identity"]
            
            for mapping_number, mapping in enumerate(mappings):
                # Figure out what the reference sequence for this mapping should
                # be
                
                position = mapping.get("position", {})
                if position.has_key("node_id"):
                    # We actually are mapped to a reference node
                    ref_sequence = node_sequences[position["node_id"]]
                    
                    # Grab the offset
                    offset = position.get("offset", 0)
                    
                    if mapping.get("is_reverse", False):
                        # We start at the offset base on the reverse strand.

                        # Add 1 to make the offset inclusive as an end poiint                        
                        ref_sequence = reverse_complement(
                            ref_sequence[0:offset + 1])
                    else:
                        # Just clip so we start at the specified offset
                        ref_sequence = ref_sequence[offset:]
                    
                else:
                    # We're aligned against no node, and thus an empty reference
                    # sequence (and thus must be all insertions)
                    ref_sequence = "" 
                    
                # Start at the beginning of the reference sequence for the
                # mapping.
                index_in_ref = 0
                    
                # Pull out the edits
                edits = mapping.get("edit", [])
                    
                for edit_number, edit in enumerate(edits):
                    # An edit may be a soft clip if it's either the first edit
                    # in the first mapping, or the last edit in the last
                    # mapping. This flag stores whether that is the case
                    # (although to actually be a soft clip it also has to be an
                    # insertion, and not either a substitution or a perfect
                    # match as spelled by the aligner).
                    may_be_soft_clip = ((edit_number == 0 and 
                        mapping_number == 0) or 
                        (edit_number == len(edits) - 1 and 
                        mapping_number == len(mappings) - 1))
                        
                    # Count up the Ns in the reference sequence for the edit. We
                    # get the part of the reference string that should belong to
                    # this edit.
                    reference_N_count = count_Ns(ref_sequence[
                        index_in_ref:index_in_ref + edit.get("from_length", 0)])
                        
                    if edit.get("to_length", 0) == edit.get("from_length", 0):
                        # Add in the length of this edit if it's actually
                        # aligned (not an indel or softclip)
                        aligned_length += edit.get("to_length", 0)
                        
                    if (not edit.has_key("sequence") and 
                        edit.get("to_length", 0) == edit.get("from_length", 0)):
                        # The edit has equal from and to lengths, but no
                        # sequence provided.

                        # We found a perfect match edit. Grab its length
                        matches += edit["from_length"]

                        # We don't care about Ns when evaluating perfect
                        # matches. VG already split out any mismatches into non-
                        # perfect matches, and we ignore the N-matched-to-N
                        # case.

                    if not may_be_soft_clip and (edit.get("to_length", 0) !=
                        edit.get("from_length", 0)):
                        # This edit is an indel and isn't on the very end of a
                        # read.
                        if reference_N_count == 0:
                            # Only count the indel if it's not against an N in
                            # the reference
                            indels += 1

                    if (edit.get("to_length", 0) == 
                        edit.get("from_length", 0) and 
                        edit.has_key("sequence")):
                        # The edit has equal from and to lengths, and a provided
                        # sequence. This edit is thus a SNP or MNP. It
                        # represents substitutions.

                        # We take as substituted all the bases except those
                        # opposite reference Ns. Sequence Ns are ignored.
                        substitutions += (edit.get("to_length", 0) -
                            reference_N_count)
                            
                        # Pull those Ns out of the substitution rate denominator
                        # as well.
                        aligned_length -= reference_N_count
                        
                    # We still count query Ns as "aligned" when not in indels
                        
                    # Advance in the reference sequence
                    index_in_ref += edit.get("from_length", 0)

            # Calculate mismatches as what's not perfect matches
            mismatches = length - matches
                    
            if alignment.get("is_secondary", False):
                # It's a multimapping. We can have max 1 per read, so it's a
                # multimapped read.
                
                if (last_alignment is None or 
                    last_alignment.get("name") != alignment.get("name") or 
                    last_alignment.get("is_secondary", False)):
                
                    # This is a secondary alignment without a corresponding primary
                    # alignment (which would have to be right before it given the
                    # way vg dumps buffers
                    raise RuntimeError("{} secondary alignment comes after "
                        "alignment of {} instead of corresponding primary "
                        "alignment\n".format(alignment.get("name"), 
                        last_alignment.get("name") if last_alignment is not None 
                        else "nothing"))
                
                # Log its stats as multimapped
                stats["total_multimapped"] += 1
                stats["secondary_scores"][score] += 1
                stats["secondary_mismatches"][mismatches] += 1
                stats["secondary_indels"][indels] += 1
                stats["secondary_substitutions"][substitutions] += 1
                stats["secondary_mapqs"][mapq] += 1
                stats["secondary_identities"][identity] += 1
            else:
                # Log its stats as primary. We'll get exactly one of these per
                # read with any mappings.
                stats["total_mapped"] += 1
                stats["primary_scores"][score] += 1
                stats["primary_mismatches"][mismatches] += 1
                stats["primary_indels"][indels] += 1
                stats["primary_substitutions"][substitutions] += 1
                stats["primary_mapqs"][mapq] += 1
                stats["primary_identities"][identity] += 1
                
                # Record that a read of this length was mapped
                stats["mapped_lengths"][length] += 1
                
                # And that a read with this many aligned primary bases was found
                stats["aligned_lengths"][aligned_length] += 1
                
                # We won't see an unaligned primary alignment for this read, so
                # count the read
                stats["total_reads"] += 1
        
        elif not alignment.get("is_secondary", False):
            # We have an unmapped primary "alignment"
            
            # Count the read by its primary alignment
            stats["total_reads"] += 1
            
            # Record that an unmapped read has this length
            stats["unmapped_lengths"][length] += 1
            
        # Save the alignment for checking for wayward secondaries
        last_alignment = alignment
                
    with open(stats_file, "w") as stats_handle:
        # Save the stats as JSON
        json.dump(stats, stats_handle)
        
    if read_alignment.wait() != 0:
        # Complain if vg dies
        raise RuntimeError("vg died with error {}".format(
            read_alignment.returncode))
        
    # Now send the stats to the output store where they belong.
    out_store.write_output_file(stats_file, stats_file_key)
    
        
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
    
    if options.too_old is not None:
        # Parse the too-old date
        options.too_old = dateutil.parser.parse(options.too_old)
        assert(options.too_old.tzinfo != None)
    
    RealTimeLogger.start_master()
    
    # Pre-read the input file so we don't try to send file handles over the
    # network.
    options.server_list = list(options.server_list)
    
    # Make a root job
    root_job = Job.wrapJobFn(run_all_alignments, options,
        cores=1, memory="4G", disk="50G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

