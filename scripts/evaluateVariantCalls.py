#!/usr/bin/env python2.7
"""
evaluateVariantCalls.py: download platinum genomes call sets, break them out by
region, and compare them against Glenn-format variant calls on vg graphs, by
using a VCF comparison tool.

"""

import argparse, sys, os, os.path, random, collections, shutil, itertools, glob
import urllib2, urlparse, ftplib, fnmatch, subprocess
import traceback
import re

from toil.job import Job
import tsv

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
    parser.add_argument("graph_store",
        help="IOStore with <graph>-<region>.vg graphs to read")
    parser.add_argument("call_store",
        help="IOStore with <region>/<graph>/<sample>_sample.txt call files")
    parser.add_argument("out_store",
        help="IOStore to save output files in")
    parser.add_argument("--blacklist", action="append", default=[],
        help="ignore the specified regions, graphs, or region:graph pairs")
    parser.add_argument("--reference_metadata", 
        default="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCA_000001405.17_GRCh38.p2/"
        "GCA_000001405.17_GRCh38.p2_assembly_structure/"
        "all_alt_scaffold_placement.txt",
        help="URL to download the reference metadata from")
    parser.add_argument("--reference_fasta", default="ftp://ftp-trace.ncbi.nih."
        "gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/"
        "GRCh38_full_analysis_set_plus_decoy_hla.fa",
        help="URL to full reference FASTA file")
    parser.add_argument("--reference_index", default="ftp://ftp-trace.ncbi.nih."
        "gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/"
        "GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
        help="URL to reference FASTA index")
    parser.add_argument("--reference_dict", default="ftp://ftp-trace.ncbi.nih."
        "gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/"
        "GRCh38_full_analysis_set_plus_decoy_hla.dict",
        help="URL to full reference FASTA dict")
    parser.add_argument("--truth_url", default="ftp://ftp-trace.ncbi.nlm.nih"
        ".gov/giab/ftp/data/NA12878/analysis/"
        "Illumina_PlatinumGenomes_NA12877_NA12878_09162015/hg38/2.0.1/"
        "{0}/{0}.vcf.gz",
        help="URL to download truth VCFs from, with {0} for the sample name")
    parser.add_argument("--gatk_url", default="https://s3-us-west-2.amazonaws"
        ".com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar",
        help="URL at which to find the GATK JAR file")
    parser.add_argument("--regions", nargs="*", 
        default=["BRCA1", "BRCA2", "CENX", "MHC", "SMA", "LRC_KIR"],
        help="region names to extract variants for (upper case)")
    parser.add_argument("--samples", action="append",
        default=["NA12877", "NA12878"],
        help="sample names to process")
    parser.add_argument("--overwrite", action="store_true",
        help="overwrite already downloaded samples")
    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def getRegions(metadata_url):
    """
    Download the assembly metadata file at the given URL, and return a dict from
    upper-case region names to 0-based end-exclusive (contig, start, end)
    tuples. Contig names start with "chr".
    
    """
    
    # Holds the chromosome number for each region?
    region_chromosomes = {}
    # Holds the minimum start position for each region on its chromosome
    region_starts = collections.defaultdict(lambda: float("inf"))
    # Holds the maximum stop position for each region on its chromosome
    region_stops = collections.defaultdict(lambda: float("-inf"))
    
    # Holds the (contig, start, end) tuple for each alt in a given region.
    ranges_by_region = collections.defaultdict(list)
    
    # Hard-code some regions that aren't real alt regions
    ranges_by_region["BRCA1"] = ("chr17", 43044294, 43125482)
    ranges_by_region["BRCA2"] = ("chr13", 32314861, 32399849)
    ranges_by_region["CENX"] = ("chrX", 58605580, 62412542)
    
    # Read the reference database
    database = tsv.TsvReader(urllib2.urlopen(metadata_url))
    
    for parts in database:
        # Parse out all the info for this alt and its parent chromosome
        region_name = parts[7]
        # Grab the chromosome ("1" or "X") that's the parent
        parent_chromosome = parts[5]
        parent_start = int(parts[11])
        parent_stop = int(parts[12])
        alt_contig = parts[3]
        alt_start = int(parts[9])
        alt_stop = int(parts[10])
        
        # Note the region start, stop, and parent chromosome number
        region_chromosomes[region_name] = parent_chromosome
        region_starts[region_name] = min(region_starts[region_name],
            parent_start)
        region_stops[region_name] = max(region_stops[region_name],
            parent_stop)
                
    for region_name in region_chromosomes.iterkeys():
        # Add in the reference ranges that all the alts are alternatives to
        # Make sure to add the chr prefix.
        ranges_by_region[region_name] = ("chr" + 
            region_chromosomes[region_name], region_starts[region_name],
            region_stops[region_name])
            
    # Give back our region info dict
    return ranges_by_region
    
# We want pickleable defaultdicts
def defaultdict_dict():
    return collections.defaultdict(dict)

def defaultdict_set():
    return collections.defaultdict(set)

def compareAllSamples(job, options):
    """
    Root job that kicks off jobs to compare all the samples.
    
    """
    
    # Make the IOStores
    graph_store = IOStore.get(options.graph_store)
    call_store = IOStore.get(options.call_store)
    out_store = IOStore.get(options.out_store)
    
    RealTimeLogger.get().info("Starting evaluation")
    
    
    # Download GATK JAR
    RealTimeLogger.get().info("Downloading GATK...")
    local_gatk = os.path.join(job.fileStore.getLocalTempDir(),
        "GenomeAnalysisTK.jar")
    shutil.copyfileobj(urllib2.urlopen(options.gatk_url), open(local_gatk, "w"))
    
    # Save it to the file store
    gatk_id = job.fileStore.writeGlobalFile(local_gatk, cleanup=True)
    
    # Download reference files
    RealTimeLogger.get().info("Downloading reference...")
    local_reference = os.path.join(job.fileStore.getLocalTempDir(),
        "ref.fa")
    shutil.copyfileobj(urllib2.urlopen(options.reference_fasta), open(local_reference,
        "w"))
    
    # Save it to the file store
    reference_id = job.fileStore.writeGlobalFile(local_reference, cleanup=True)
    
    # And the index file
    local_reference_index = os.path.join(job.fileStore.getLocalTempDir(),
        "ref.fa.fai")
    shutil.copyfileobj(urllib2.urlopen(options.reference_index),
        open(local_reference_index, "w"))
    reference_index_id = job.fileStore.writeGlobalFile(local_reference_index,
        cleanup=True)
        
    # And the reference dict file
    local_reference_dict = os.path.join(job.fileStore.getLocalTempDir(),
        "ref.dict")
    shutil.copyfileobj(urllib2.urlopen(options.reference_dict),
        open(local_reference_dict, "w"))
    reference_dict_id = job.fileStore.writeGlobalFile(local_reference_dict,
        cleanup=True)
    
    
    RealTimeLogger.get().info("Scanning graphs...")
    
    # Work out what regions are where. Note that thgis is keyed by *upper case*
    # region names, while everything else uses lower case names.
    ranges_by_region = getRegions(options.reference_metadata)
    
    # Look through all the graphs and see what graphs we have.
    graph_keys_by_name_and_region = collections.defaultdict(dict)
    for graph_file in graph_store.list_input_directory(""):
        # For each graph, try and pull out its name and region
        match = re.match("^([^\\-]+)-([^\\-]+)-?([^\\-]*).vg$", graph_file)
        
        if match:
            # We found a graph
            
            # Everything around the region is the graph name
            graph_name_parts = [match.group(1)]
            if match.group(3) != "":
                graph_name_parts.append(match.group(3))
            graph_name = "-".join(graph_name_parts)
            
            # The part in the middle (or last if no trailing part of the name)
            # is the region.
            graph_region = match.group(2)
            
            # Save the graph so we can look it up later.
            graph_keys_by_name_and_region[graph_name][graph_region] = graph_file
    
    # We're going to queue up VCF-ing all the samples for all the graphs and
    # downloading the truth for the region in parallel. The stats job for each
    # sample is going to be a follow-on of both the truth download job and the
    # VCF conversion job.
    
    RealTimeLogger.get().info("Scanning samples...")
    
    # This holds all the samples that have any graphs for a region.
    # We'll need to fetch the truth for each entry.
    samples_for_region = collections.defaultdict(set)
    
    # Find the actual sample call files that exist. Stores sets of sample names
    # by region, then graph.
    found_samples = collections.defaultdict(defaultdict_set)
    
    # Loop through all the sample call files
    for region_dir in call_store.list_input_directory(""):
        # Within every region we have samples for, look through all the
        # different graphs.
        for graph_dir in call_store.list_input_directory(region_dir):
            # Within every graph for a region, we have a collection of samples.
            
            if ("{}:{}".format(region_dir, graph_dir) in options.blacklist or
                region_dir in options.blacklist or
                graph_dir in options.blacklist):
                # We don't want to process this region/graph pair.
                RealTimeLogger.get().info("Skipping {} graph {}".format(
                    region_dir, graph_dir))
                continue
                
            for filename in call_store.list_input_directory("{}/{}".format(
                region_dir, graph_dir)):
                # Look at each potential sample file
                
                # Is this file a sample?
                match = re.match("(.+)_sample\\.txt$", filename)
                
                if not match:
                    # It's not a sample
                    continue
                    
                # If we did match, pull out the sample name.
                sample_name = match.group(1)
                
                if (options.samples is not None and
                    sample_name not in options.samples):
                    # Don't do this sample. Only do samples we explicitly asked
                    # for.
                    continue
                
                # Note that this sample has any graphs in this region.
                samples_for_region[region_dir].add(sample_name)
                
                # And that it's under this graph for this region
                found_samples[graph_dir][region_dir].add(sample_name)
    
    # Make all the truth download jobs, by sample
    truth_jobs = collections.defaultdict(dict)
    
    for region_dir, region_samples in samples_for_region.iteritems():
        # For each region
        for sample_name in region_samples:
            # For each sample in the region, download the truth
            
            # Unpack the range parameters
            ref_name, ref_start, ref_end = ranges_by_region[
                region_dir.upper()]
            
            # Do the download, returning a file store ID for a VCF file.
            truth_jobs[region_dir.lower()][sample_name] = job.addChildJobFn(
                downloadTruth, options, sample_name, region_dir, ref_name,
                ref_start, ref_end, cores=1, memory="10G", disk="4G")
            
    
    # Make jobs to convert called samples to VCF. Stores jobs by region, graph,
    # and then sample name.
    conversion_jobs = collections.defaultdict(defaultdict_dict)
        
    # Make another dict for the promises
    conversion_promises = collections.defaultdict(defaultdict_dict)
        
    # Make follow-ons that depend on the conversion and the truth download and
    # actually compare things.
    comparison_jobs = collections.defaultdict(defaultdict_dict)
        
    for graph_dir, regions_for_graph in found_samples.iteritems():
        # For each graph we actually have samples for...
        for region_dir, sample_set in regions_for_graph.iteritems():
            # For each region we actually have samples for on that graph...
            
            # Download and store the graph in the file store, because each
            # sample will need it.
            
            # Download the pre-made index directory
            # What is the graph filename? We assume the vg file exists already.
            graph_basename = graph_keys_by_name_and_region[graph_dir][
                region_dir]
            # Where should it sit locally?
            graph_file = "{}/{}".format(job.fileStore.getLocalTempDir(),
                graph_basename)
            # Download it
            graph_store.read_input_file(graph_basename, graph_file)
            
            # Save it to the global file store and keep around the ID.
            graph_file_id = job.fileStore.writeGlobalFile(graph_file,
                cleanup=True)
            
            RealTimeLogger.get().info("Saved graph {} to {}".format(
                graph_basename, graph_file_id)) 
        
            for sample_name in sample_set:
                # For each sample we actually found...
                
                # Construct the key for the sample's call file
                sample_key = os.path.join(region_dir, graph_dir,
                    "{}_sample.txt".format(sample_name))
                
                # Construct the key where we want the resulting VCF to go
                vcf_key = os.path.join("calls", region_dir, graph_dir,
                    "{}.vcf".format(sample_name))
                
                # Unpack the range parameters
                ref_name, ref_start, _ = ranges_by_region[region_dir.upper()]
                
                # Since ref_start is used as an offset, we can pass a 0-based
                # position in.
                
                # Add a conversion job
                conversion_jobs[region_dir][graph_dir][sample_name] = \
                    job.addChildJobFn(convertGlennToVcf, options, 
                    sample_key, graph_file_id, vcf_key, ref_name, ref_start,
                    cores=1, memory="10G", disk="4G")
                
                # Save the promise for the conversion
                conversion_promises[region_dir][graph_dir][sample_name] = \
                    conversion_jobs[region_dir][graph_dir][sample_name].rv()
                    
                # Add real comparison job
                # Just do it after us, don't do fine-grained dependencies.
                comparison_job = job.addFollowOnJobFn(compareVcfs, options, gatk_id,
                    reference_id, reference_index_id, reference_dict_id,
                    truth_jobs[region_dir.lower()][sample_name].rv(), 
                    conversion_promises[region_dir][graph_dir][sample_name],
                    "stats/gatk/{}/{}/{}.grp".format(region_dir, graph_dir,
                    sample_name))
                    
                
    # Now that the promises dict is populated, add the saving-number-of-bases-
    # dropped job.
    log_dropped_bases_job = None    
    
    for region_name, by_graph in conversion_jobs.iteritems():
        # For each region and the graphs in it
        for graph_name, by_sample in by_graph.iteritems():
            # For each graph and the samples on it
            for sample_name, stats in by_sample.iteritems():
                # For each sample
                
                if log_dropped_bases_job is None:
                    # Make the job to do the logging
                    log_dropped_bases_job = job.addFollowOnJobFn(saveBasesDropped, options, 
                        conversion_promises, "stats/bases_dropped.tsv",
                        cores=1, memory="10G", disk="4G")
                        
    RealTimeLogger.get().info("Done making children")
   
def downloadTruth(job, options, sample_name, region_name, ref_name, ref_start,
    ref_end):
    """
    Download the specified range of the VCF for the specified sample, store it
    in the file store, and return the ID.
    
    """
    
    # Make the IOStores
    graph_store = IOStore.get(options.graph_store)
    call_store = IOStore.get(options.call_store)
    out_store = IOStore.get(options.out_store)
    
    RealTimeLogger.get().info("Get truth for {}:{}-{} for {}".format(ref_name,
        ref_start, ref_end, sample_name))
    
    local_filename = os.path.join(job.fileStore.getLocalTempDir(),
        "sample.vcf")
        
    # Open it for writing
    vcf_file = open(local_filename, "w")
        
    # Download the region to the file. TODO: 1-based or 0-based range?
    subprocess.check_call(["bcftools", "view",
        options.truth_url.format(sample_name), "-r",
        "{}:{}-{}".format(ref_name, ref_start, ref_end)], stdout=vcf_file)
        
    # Close the file
    vcf_file.close()
    
    # Upload it
    file_id = job.fileStore.writeGlobalFile(local_filename)
    
    # Save it in our output directory
    out_store.write_output_file(local_filename, "truth/{}/{}.vcf".format(
        region_name, sample_name))
    
    return file_id
   
def convertGlennToVcf(job, options, glenn_file_key, graph_file_id,
    vcf_file_key, ref_name, ref_start):
    """
    Given the key for a Glenn-format variant file, and the graph file ID in the
    file store, convert the Glenn file to a VCF. Pust variants on a reference
    named ref_name, offset from starting at base 1 by ref_start.
    
    Returns the file store ID of the converted VCF file, which is also saved
    under the given key, and the number of bases of novel material not
    converted.
    
    """
    
    # Make the IOStores
    graph_store = IOStore.get(options.graph_store)
    call_store = IOStore.get(options.call_store)
    out_store = IOStore.get(options.out_store)
    
    RealTimeLogger.get().info("Bang {} against {} into {}".format(
        glenn_file_key, graph_file_id, vcf_file_key))
        
    # Download the input files
    # Glenn file
    local_glenn_file = os.path.join(job.fileStore.getLocalTempDir(),
        "glenn.txt")
    call_store.read_input_file(glenn_file_key, local_glenn_file)
    # And graph
    local_graph_file = job.fileStore.readGlobalFile(graph_file_id)
    
    # Define a local file for the VCF
    local_vcf_file = os.path.join(job.fileStore.getLocalTempDir(),
        "out.vcf")
        
    # Open it for writing
    vcf_file = open(local_vcf_file, "w")
    
    # Run the conversion, send stdout to the file, and give us stderr.
    # "ref_name" is actually the name we *want* to have after conversion, not
    # the name we have now.
    
    # TODO: we subtract 1 off the ref start here because we agree at all with
    # the truth set VCF only if we do that. I'm not sure where exactly we're
    # fixing a coordinate mis-conversion.
    command = ["glenn2vcf", local_graph_file, local_glenn_file, "--contig",
        ref_name, "--offset", str(ref_start - 1)]
    conversion = subprocess.Popen(command,
        stdout=vcf_file, stderr=subprocess.PIPE)
        
    # How many unrepresentable bases were dropped?
    bases_dropped = None
        
    for line in conversion.stderr:
        # Look for the line listing dropped bases
        match = re.match(
            "Had to drop ([0-9]+) bp of unrepresentable variation.", line)
            
        if match:
            # We found the right line, parse out the base count
            bases_dropped = int(match.group(1))
            
        # Pass the line through to our standard error in case something goes
        # wrong.
        sys.stderr.write(line)
        
    if conversion.wait() != 0:
        # Complain if the conversion fails
        
        # Save the files that confused us
        out_store.write_output_file(local_graph_file, "error/graph.vg")
        out_store.write_output_file(local_glenn_file, "error/glenn.txt")
        
        raise RuntimeError("VCF conversion of {} via '{}' returned {}".format(
            glenn_file_key, " ".join(command), conversion.returncode))
            
    vcf_file.close()
    
    # Now upload the VCF
    out_store.write_output_file(local_vcf_file, vcf_file_key)
    
    # Save it to the file store too
    vcf_file_id = job.fileStore.writeGlobalFile(local_vcf_file)
        
    # Return the VCF file ID and the number of bases dropped, as specified.
    return (vcf_file_id, bases_dropped)
    
def saveBasesDropped(job, options, return_value_dict, out_key):
    """
    Given a dict from region name, graph name, and sample name to the retrun
    value tuple of convertGlennToVcf, which is (vcf_file_id, pases_dropped),
    save a TSV in the format <region>\t<graph>\t<sample>\t<bases dropped> to the
    given key.
    
    """
    
    # We need to save the TSV file somewhere
    local_filename = os.path.join(job.fileStore.getLocalTempDir(),
        "stats.tsv")
        
    writer = tsv.TsvWriter(open(local_filename, "w"))
    
    # Make the IOStores
    graph_store = IOStore.get(options.graph_store)
    call_store = IOStore.get(options.call_store)
    out_store = IOStore.get(options.out_store)
    
    for region_name, by_graph in return_value_dict.iteritems():
        # For each region and the graphs in it
        for graph_name, by_sample in by_graph.iteritems():
            # For each graph and the samples on it
            for sample_name, stats in by_sample.iteritems():
                # For each sample and its return value from saveBasesDropped
                
                # Save the bases dropped
                writer.line(region_name, graph_name, sample_name, stats[1])
                
    
    # Save the aggregated output    
    writer.close()  
    out_store.write_output_file(local_filename, out_key)  
    
def compareVcfs(job, options, gatk_id, reference_id, reference_index_id,
    reference_dict_id, truth_id, query_id_and_count, out_key):
    """
    Given the GATK JAR file ID, reference FASTA ID, ID of the index for said
    reference FASTA, ID of the "dict" file for said reference FASTA, and two VCF
    file IDs, compare them.
    
    Save the results of the comparison to the given output key in the output
    IOStore.
    
    """
    
    # Make the IOStores
    graph_store = IOStore.get(options.graph_store)
    call_store = IOStore.get(options.call_store)
    out_store = IOStore.get(options.out_store)
    
    # Download the GTAK
    gatk_jar = os.path.join(job.fileStore.getLocalTempDir(),
        "GenomeAnalysisTK.jar")
    job.fileStore.readGlobalFile(gatk_id, userPath=gatk_jar)
    
    # Download the reference FASTA. We need to specify names manually because
    # GATK explodes unless you feed it a particular folder structure.
    reference_fasta = os.path.join(job.fileStore.getLocalTempDir(),
        "ref.fa")
    job.fileStore.readGlobalFile(reference_id, userPath=reference_fasta)
    
    # And the FASTA index
    reference_index = reference_fasta + ".fai"
    job.fileStore.readGlobalFile(reference_index_id, userPath=reference_index)
    
    # And the FASTA dict
    reference_dict = os.path.splitext(reference_fasta)[0] + ".dict"
    job.fileStore.readGlobalFile(reference_dict_id, userPath=reference_dict)
    
    # Download the VCFs
    truth_filename = job.fileStore.readGlobalFile(truth_id)
    query_filename = job.fileStore.readGlobalFile(query_id_and_count[0])
    
    # Now sort the VCFs
    truth_sorted = os.path.join(job.fileStore.getLocalTempDir(),
        "truth.sorted.vcf")
    query_sorted = os.path.join(job.fileStore.getLocalTempDir(),
        "query.sorted.vcf")
        
    # TODO: make this file come with the script on Mesos
    subprocess.check_call(["scripts/vcfsorter.pl", reference_dict,
        truth_filename], stdout=open(truth_sorted, "w"))
    subprocess.check_call(["scripts/vcfsorter.pl", reference_dict,
        query_filename], stdout=open(query_sorted, "w"))
        
    # Decide on an output filename for GTAK's weird table
    output_filename = os.path.join(job.fileStore.getLocalTempDir(),
        "output.grp")
        
    # Run GATK
    RealTimeLogger.get().info("Running GATK VCF comparison for {}...".format(
        out_key))
    
    
    try:
        subprocess.check_call(["java", "-Xmx8g", "-jar", gatk_jar, "-T",
            "GenotypeConcordance", "-R", reference_fasta, "--eval", query_sorted,
            "--comp", truth_sorted, "-o", output_filename])
    except:
        # If GATK fails, save out input VCFs
        RealTimeLogger.get().error(
            "GATK failed; saving input files for debugging")
        out_store.write_output_file(reference_fasta, "error/ref.fa")
        out_store.write_output_file(query_sorted, "error/query.vcf")
        out_store.write_output_file(truth_sorted, "error/truth.vcf")
        raise
        
    # Upload the answer
    out_store.write_output_file(output_filename, out_key)
    RealTimeLogger.get().info("Generated {}".format(out_key))
    
        
    
def main():
    options = parse_args(sys.argv) # This holds the nicely-parsed options object
    
    RealTimeLogger.start_master()
    
    # Make the root job
    root_job = Job.wrapJobFn(compareAllSamples, options, 
        cores=1, memory="1G", disk="4G")
        
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
        
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()
    
if __name__=="__main__":
    sys.exit(main())
