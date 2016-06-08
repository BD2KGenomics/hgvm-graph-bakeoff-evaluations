#!/usr/bin/env python2.7
"""
parallelCallingEvaluation.py: Run the mapping and variant calling evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

example run: ./parallelCallingEvaluation.py --batchSystem mesos --mesosMaster 10.0.0.5:5050 --realTimeLogging --logError --logDebug --edge_max 5 --kmer_size 16 --index_mode gcsa-mem --include_primary '/home/cmarkello/hgvmeval-jobstore1' '/home/cmarkello/debug_eval_input/BRCA1.vg' 'ref' '81189' '/home/cmarkello/debug_eval_input/BRCA1/NA12877/NA12877.bam.fq' 'NA12877' '/home/cmarkello/debug_eval_output'

example run: ./parallelCallingEvaluation.py --offset 43044293 --edge_max 5 --kmer_size 16 --index_mode gcsa-mem --include_primary '/home/cmarkello/debug_eval_input/BRCA1.vg' 'ref' '81189' '/home/cmarkello/debug_eval_input/BRCA1/NA12877/NA12877.bam.fq' 'NA12877' '/home/cmarkello/debug_eval_output'
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import ntpath
import pdb

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
    #Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("vg_graph", type=str,
        help="Input vg graph file path")
    parser.add_argument("path_name",
        help="Name of reference path in the graph (eg. ref or 17)")
    parser.add_argument("path_size",
        help="Size of the reference path in the graph")
    parser.add_argument("sample_reads", type=str,
        help="Path to sample reads in fastq format")
    parser.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser.add_argument("out_dir", type=str,
        help="directory where all output will be written")
    parser.add_argument("--offset", type=int,
        help="chromosomal position offset. e.g. 43044293")
    parser.add_argument("--edge_max", type=int, default=5,
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=16,
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--overwrite", default=False, action="store_true",
        help="overwrite existing result files")
    parser.add_argument("--restat", default=False, action="store_true",
        help="recompute and overwrite existing stats files")
    parser.add_argument("--reindex", default=False, action="store_true",
        help="don't re-use existing indexed graphs")
    parser.add_argument("--index_mode", choices=["rocksdb", "gcsa-kmer",
        "gcsa-mem"], default="gcsa-mem",
        help="type of vg index to use for mapping")
    parser.add_argument("--include_pruned", action="store_true",
        help="use the pruned graph in the index")
    parser.add_argument("--include_primary", action="store_true",
        help="use the primary path in the index")
    parser.add_argument("--pileup_opts", type=str,
        default="-w 40 -m 10 -q 10",
        help="options to pass to vg pileup. wrap in \"\"")
    parser.add_argument("--call_opts", type=str,
        default="-r 0.0001 -b 0.4 -f 0.25 -d 10",
        help="options to pass to vg call. wrap in \"\"")

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

def run(cmd, proc_stdout = sys.stdout, proc_stderr = sys.stderr,
        check = True):
    """ run command in shell and throw exception if it doesn't work 
    """
    print(cmd)
    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=proc_stdout, stderr=proc_stderr)
    output, errors = proc.communicate()
    sts = proc.wait()
    if check is True and sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (cmd, sts))
    return output, errors

#def run_indexing(job, options):
def run_indexing(options, job_cores):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    print("Starting indexing...")

    graph_filename = ntpath.basename(options.vg_graph)



    # Now run the indexer.
    # TODO: support both indexing modes
    print("Indexing {}".format(options.vg_graph))

    if options.index_mode == "rocksdb":
        # Make the RocksDB index
        run("vg index -s -k {} -e {} -t {} {} -d {}/{}.index".format(options.kmer_size, options.edge_max, job_cores, options.vg_graph, options.out_dir, graph_filename))

    elif (options.index_mode == "gcsa-kmer" or
        options.index_mode == "gcsa-mem"):
        # We want a GCSA2/xg index. We have to prune the graph ourselves.
        # See <https://github.com/vgteam/vg/issues/286>.

        # What will we use as our temp combined graph file (containing only
        # the bits of the graph we want to index, used for deduplication)?
        to_index_filename = "{}/to_index.vg".format(
            options.out_dir)

        # Where will we save the kmers?
        kmers_filename = "{}/index.graph".format(
            options.out_dir)

        with open(to_index_filename, "w") as to_index_file:

            if options.include_pruned:

                print("Pruning {} to {}".format(
                    graph_filename, to_index_filename))

                # Prune out hard bits of the graph
                tasks = []

                # Prune out complex regions
                tasks.append(subprocess.Popen(["vg", "mod",
                    "-p", "-l", str(options.kmer_size), "-t", str(job_cores),
                    "-e", str(options.edge_max), options.vg_graph],
                    stdout=subprocess.PIPE))

                # Throw out short disconnected chunks
                tasks.append(subprocess.Popen(["vg", "mod",
                    "-S", "-l", str(options.kmer_size * 2),
                    "-t", str(job_cores), "-"], stdin=tasks[-1].stdout,
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

                print(
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
                    ["vg", "mod", "-N"] + ref_options +
                    ["-t", str(job_cores), options.vg_graph],
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

        tasks = []
        with open(kmers_filename, "w") as kmers_file:

            tasks = []

            print("Finding kmers in {} to {}".format(
                to_index_filename, kmers_filename))

            # Deduplicate the graph
            tasks.append(subprocess.Popen(["vg",
                "view", "-v", to_index_filename],
                stdout=subprocess.PIPE))
            
            # Make the GCSA2 kmers file
            tasks.append(subprocess.Popen(["vg",
                "kmers", "-g", "-B", "-k", str(options.kmer_size),
                "-H", "1000000000", "-T", "1000000001",
                "-t", str(job_cores), "-"], stdin=tasks[-1].stdout,
                stdout=kmers_file))
                
            # Did we make it through all the tasks OK?
            for task in tasks:
                if task.wait() != 0:
                    raise RuntimeError("Pipeline step returned {}".format(
                        task.returncode))

            # Wait to make sure no weird file-not-being-full bugs happen
            # TODO: how do I wait on child process output?
            time.sleep(1)

        time.sleep(1)

        # Where do we put the GCSA2 index?
        gcsa_filename = options.out_dir + "/" + graph_filename + ".gcsa"

        print("GCSA-indexing {} to {}".format(
                kmers_filename, gcsa_filename))

        # Make the gcsa2 index. Make sure to use 3 doubling steps to work
        # around <https://github.com/vgteam/vg/issues/301>
        subprocess.check_call(["vg", "index", "-t", str(job_cores), "-i",
            kmers_filename, "-g", gcsa_filename, "-X", "3"])

        # Where do we put the XG index?
        xg_filename = options.out_dir + "/" + graph_filename + ".xg"

        print("XG-indexing {} to {}".format(
                options.vg_graph, xg_filename))

        subprocess.check_call(["vg", "index", "-t", str(job_cores), "-x",
            xg_filename, options.vg_graph])
    
    # Define a file to keep the compressed index in, so we can send it to
    # the output store.
    index_dir_tgz = "{}/index.tar.gz".format(
        options.out_dir)

    
#def run_alignment(job, options, bin_dir_id, sample, graph_name, region,
#    index_dir_id, sample_fastq_key, alignment_file_key, stats_file_key):
def run_alignment(options, job_cores):

    graph_filename = ntpath.basename(options.vg_graph)
    
    graph_dir = options.out_dir
    
    # How long did the alignment take to run, in seconds?
    run_time = None

    # We know what the vg file in there will be named
    graph_file = options.vg_graph

    # Also we need the sample fastq
    fastq_file = options.sample_reads

    # And a temp file for our aligner output
    output_file = "{}/output.gam".format(graph_dir)

    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = ["vg", "map", "-f", fastq_file,
            "-i", "-M2", "-a", "-u", "0", "-U", "-t", str(job_cores), graph_file]

        if options.index_mode == "rocksdb":
            vg_parts += ["-d", graph_dir+"/"+graph_filename+".index", "-n3", "-k",
                str(options.kmer_size)]
        elif options.index_mode == "gcsa-kmer":
            # Use the new default context size in this case
            vg_parts += ["-x", graph_dir+"/"+graph_filename+ ".xg", "-g", graph_dir+"/"+graph_filename + ".gcsa",
                "-n5", "-k", str(options.kmer_size)]
        elif options.index_mode == "gcsa-mem":
            # Don't pass the kmer size, so MEM matching is used
            vg_parts += ["-x", graph_dir+"/"+graph_filename+ ".xg", "-g", graph_dir+"/"+graph_filename + ".gcsa",
                "-n5"]
        else:
            raise RuntimeError("invalid indexing mode: " + options.index_mode)

        print(
            "Running VG for {} against {}: {}".format(options.sample_name, graph_filename,
            " ".join(vg_parts)))

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


    print("Aligned {}. Process took {} seconds.".format(output_file, run_time))


#def run_stats(job, options, bin_dir_id, index_dir_id, alignment_file_key,
#    stats_file_key, run_time=None):
def run_stats(options, job_cores):
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

    print("Computing stats for {}".format(options.sample_name))


    graph_filename = ntpath.basename(options.vg_graph)
    
    graph_dir = options.out_dir
    
    # How long did the alignment take to run, in seconds?
    run_time = None

    # We know what the vg file in there will be named
    graph_file = options.vg_graph

    # Load the node sequences into memory. This holds node sequence string by
    # ID.
    node_sequences = {}

    # Read the alignments in in JSON-line format
    read_graph = subprocess.Popen(["vg", "view", "-j",
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
    stats_file = "{}/stats.json".format(options.out_dir)
    alignment_file = "{}/output.gam".format(options.out_dir)

    # Read the alignments in in JSON-line format
    read_alignment = subprocess.Popen(["vg", "view", "-aj",
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
        "primary_mismatches": collections.Counter(),
        "primary_indels": collections.Counter(),
        "primary_substitutions": collections.Counter(),
        "secondary_scores": collections.Counter(),
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
            else:
                # Log its stats as primary. We'll get exactly one of these per
                # read with any mappings.
                stats["total_mapped"] += 1
                stats["primary_scores"][score] += 1
                stats["primary_mismatches"][mismatches] += 1
                stats["primary_indels"][indels] += 1
                stats["primary_substitutions"][substitutions] += 1

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

def run_calling(options, job_cores):
    
    graph_filename = ntpath.basename(options.vg_graph)
    path_name = options.path_name
    vg_path = options.vg_graph
    gam_path = options.out_dir + "/output.gam"

    # do the pileup.  this is the most resource intensive step,
    # especially in terms of mermory used.
    pu_path = options.out_dir + "/" + options.sample_name + ".pu"
    run("vg pileup {} {} -t {} {} > {}".format(
        vg_path, gam_path, job_cores, options.pileup_opts, pu_path))

    # do the calling.
    tsv_path = options.out_dir + "/" + options.sample_name + "_call.tsv"
    ag_path = options.out_dir + "/" + options.sample_name + "_call.vg"
    run("vg call {} {} -t {} {} -l -c {} > {}".format(
        vg_path, pu_path, job_cores, options.call_opts, tsv_path, ag_path))

    # do the vcf export.
    vcf_path = options.out_dir + "/" + options.sample_name + ".vcf"
    xg_path = options.out_dir + "/" + graph_filename + ".xg"
    offset = xg_path_node_offset(xg_path, options.path_name, options.offset)
    run("glenn2vcf {} {} -o {} -c {} -s {} -l {} > {} 2> {}".format(
        ag_path, tsv_path, options.offset, options.path_name, options.sample_name, options.path_size,
        vcf_path + ".us", vcf_path + ".log"))
    sort_vcf(vcf_path + ".us", vcf_path)
    run("rm {}".format(vcf_path + ".us"))
    run("bgzip {}".format(vcf_path))
    run("tabix -f -p vcf {}".format(vcf_path + ".gz"))

def xg_path_node_id(xg_path, path_name, offset):
    """ use vg find to get the node containing a given path position """
    #NOTE: vg find -p range offsets are 0-based inclusive.  
    stdout, stderr = run("vg find -x {} -p {}:{}-{} | vg mod -o - | vg view -j - | jq .node[0].id".format(
        xg_path, path_name, offset, offset),
                         proc_stdout=subprocess.PIPE)
    return int(stdout)

def xg_path_node_offset(xg_path, path_name, offset):
    """ get the offset of the node containing the given position of a path
    """
    # first we find the node
    node_id = xg_path_node_id(xg_path, path_name, offset)

    # now we find the offset of the beginning of the node
    stdout, stderr = run("vg find -x {} -P {} -n {}".format(
        xg_path, path_name, node_id),
                         proc_stdout=subprocess.PIPE)
    toks = stdout.split()
    # if len > 2 then we have a cyclic path, which we're assuming we don't
    assert len(toks) == 2
    assert toks[0] == str(node_id)
    node_offset = int(toks[1])
    # node_offset must be before
    assert node_offset <= offset
    # sanity check (should really use node size instead of 1000 here)
    assert offset - node_offset < 1000

    return node_offset

def sort_vcf(vcf_path, sorted_vcf_path):
    """ from vcflib """
    run("head -10000 {} | grep \"^#\" > {}".format(
        vcf_path, sorted_vcf_path))
    run("grep -v \"^#\" {} | sort -k1,1d -k2,2n >> {}".format(
        vcf_path, sorted_vcf_path))
    

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
    
    run_indexing(options, job_cores=32)
    
    run_alignment(options, job_cores=32)
    
    run_stats(options, job_cores=2)

    run_calling(options, job_cores=32)

    print("All jobs completed successfully")

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
