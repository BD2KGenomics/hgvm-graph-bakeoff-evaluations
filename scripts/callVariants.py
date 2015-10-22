#!/usr/bin/env python2.7
"""
Generates variant calls from read alignments in various ways for comparison. Input will be set of GAM alignemnts to process.  Directory structure must be consistent with Adam's scripts.  ie

alignments stored like this:
  alignments/brca1/cactus/NA19240.gam

so:
  <alignments dir>/<region>/<graph method>/<sample>.gam

The graphs must also be accessible:

  graphs/cactus-brca1.vg

so:
  <graphs dir>/<graph method>-<region>.vg

note: debruin is exception where -k63 tag gets tacked on at end (treated as special case)

all output will be in the variants/ folder.

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string
from toil.job import Job
from toillib import RealTimeLogger, robust_makedirs

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("in_gams", nargs="+",
                        help="input alignment files")
    parser.add_argument("--out_dir", type=str, default="variants",
                        help="output directory")
    parser.add_argument("--graph_dir", type=str, default="graphs",
                        help="name of input graphs directory")
    parser.add_argument("--index_ext", type=str, default=".index",
                        help="extension to find input grpah index")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite files if dont exist")
    parser.add_argument("--fa_path", type=str, default="./altRegions",
                        help="path to search for reference sequences. expects "
                        "references to be in <fa_path>/BRCA1/ref.fa, etc. "
                        " [TODO: support different ref versions?]" )
    parser.add_argument("--vg_cores", type=int, default=1,
                        help="number of cores to give to vg commands")
    parser.add_argument("--vg_only", action="store_true", default=False,
                        help="only do vg call. skip samptools linear call")
    parser.add_argument("--skip", type=str, default="",
                        help="comma-separated list of keywords that "
                        "will cause input gam to be skipped if found in path")
    
    args = args[1:]
        
    return parser.parse_args(args)

def temp_path(options, prefix="tmp", ext="", length=6):
    """ get a temporary file in out_dir/temp
    """
    tempdir = os.path.join(options.out_dir, "temp")
    robust_makedirs(tempdir)
    tag = "".join([random.choice(
        string.ascii_uppercase + string.digits) for i in xrange(length)])
    return os.path.join(tempdir, prefix + tag + ext)


def graph_path(alignment_path, options):
    """ get the graph corresponding to a gam
    """
    region = alignment_region_tag(alignment_path, options)
    graph = alignment_graph_tag(alignment_path, options)
    if graph.find("debruijn") == 0:
        assert graph.find("-") == 8
        tag = graph[8:]
        graph = "debruijn"
    else:
        tag = ""
    return os.path.join(options.graph_dir, "{}-{}{}.vg".format(graph, region, tag))
    
def index_path(graph_path, options):
    """ get index corresponding to a graph
    """
    return "{}{}".format(graph_path, options.index_ext)

def ref_path(alignment_path, options):
    """ get the fasta file corresponding to the reference, assuming region
    name comes after 1st dash in filename
    """
    region = alignment_region_tag(alignment_path, options)
    region = region.upper()
    # NOTE: this logic only supports one reference (ie grchg38), so won't work
    # for simons data.
    return os.path.join(options.fa_path, region, "ref.fa")

def alignment_read_tag(alignment_path, options):
    """ say alignment is bla/gcsa/real/camel-brca1.gam, then return real
    """
    return "real"

def alignment_map_tag(alignment_path, options):
    """ say alignment is bla/gcsa/real/camel-brca1.gam, then return gcsa
    """
    return "kmer"

def alignment_region_tag(alignment_path, options):
    """ say alignment is bla/gcsa/real/camel-brca1.gam, then return brca1
    """
    region = alignment_path.split("/")[-3]
    assert region in ["brca1", "brca2", "cenx", "lrc_kir", "sma", "mhc"]
    return region

def alignment_graph_tag(gam_path, options):
    """ extract the graph method name from gam path
    """
    return gam_path.split("/")[-2]

def alignment_sample_tag(gam_path, options):
    """ get the sample name from gam path
    """
    return os.path.splitext(os.path.basename(gam_path))[0]

def out_dir(alignment_path, options):
    """ get directory to put output corresponding to input file
    """
    return os.path.join(options.out_dir,
                        alignment_region_tag(alignment_path, options),
                        alignment_graph_tag(alignment_path, options))

def pileup_path(alignment_path, options, tag=""):
    """ get output pileup name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}.vgpu".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def sample_vg_path(alignment_path, options, tag=""):
    """ get vg call output name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_sample.vg".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def augmented_vg_path(alignment_path, options, tag=""):
    """ get output augmented variant graph name from vg call -l
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_augmented.vg".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def projected_bam_path(alignment_path, options, tag=""):
    """ get output of surjecting input gam to bam on ref
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_ref.bam".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def linear_vcf_path(alignment_path, options, tag=""):
    """ get output of samtools (linear) variant calling pipelines
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_linear.vcf".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def linear_vg_path(alignment_path, options, tag=""):
    """ get output of vg conversion from linear vcf
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_linear.vg".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)    

def run(cmd, stdout = sys.stdout, stderr = sys.stderr):
    """ run command in shell and barf if it doesn't work
    (copied from system() in sonlib.bioio
    """
    RealTimeLogger.get().info("RUN: {}".format(cmd))

    sts = subprocess.call(cmd, shell=True, bufsize=-1,
                          stdout=stdout, stderr=stderr)
    if sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" %
                           (cmd, sts))
    return sts

def compute_linear_variants(job, input_gam, options):
    """ project to bam, then run samtools to call some variants
    """
    input_graph_path = graph_path(input_gam, options)
    input_index_path = index_path(input_graph_path, options)

    # can only do this if there is a "ref" path in the vg graph
    res_path = temp_path(options)
    run("./vgHasPath.sh {} {} > {}".format(input_graph_path, "ref", res_path))
    has_ref = False
    with open(res_path) as res_file:
        has_ref = res_file.read()[0] == "1"
    run("rm {}".format(res_path))
    
    if has_ref:
        surject_path = projected_bam_path(input_gam, options)
        out_vcf_path = linear_vcf_path(input_gam, options)
        out_vg_path = linear_vg_path(input_gam, options)
        fasta_path = ref_path(input_gam, options)
        do_surject = options.overwrite or not os.path.isfile(surject_path)
        do_vcf = do_surject or not os.path.isfile(out_vcf_path + ".gz")
        do_vg = do_vcf or not os.path.isfile(out_vg_path)

        if do_surject:
            robust_makedirs(os.path.dirname(surject_path))
            prefix_path = temp_path(options, ".prefix")
            # surject to reference path (name hardcoded to ref for now)
            run("vg surject -d {} -p {} -b {} -t {} | samtools sort -o - {}> {}".format(
                input_index_path,
                "ref",
                input_gam,
                options.vg_cores,
                prefix_path,
                surject_path))
            run("rm -f {}".format(prefix_path))

        if do_vcf:
            # todo: we assume that all graphs have same reference fasta, here.
            # this is false for, ex, simons which uses grchg37 instead of 38.

            # create pileup in bcf using samtools
            # http://samtools.sourceforge.net/mpileup.shtml
            assert os.path.isfile(fasta_path)
            robust_makedirs(os.path.dirname(out_vcf_path))
            run("samtools mpileup -I -u -t DP -f {} {} | bcftools call -m -V indels - > {}".format(
                fasta_path,
                surject_path,
                out_vcf_path))

            # make compressed index
            run("bgzip -f {}".format(out_vcf_path))
            run("tabix -f -p vcf {}.gz".format(out_vcf_path))

        if do_vg:
            # and convert back to vg...
            robust_makedirs(os.path.dirname(out_vg_path))
            run("vg construct -v {}.gz -r {} -t {} > {}".format(out_vcf_path, fasta_path,
                                                                options.vg_cores, out_vg_path))
            
    
def compute_vg_variants(job, input_gam, options):
    """ run vg pileup and vg call on the input
    """
    input_graph_path = graph_path(input_gam, options)
    out_pileup_path = pileup_path(input_gam, options)
    out_sample_vg_path = sample_vg_path(input_gam, options)
    out_augmented_vg_path = augmented_vg_path(input_gam, options)

    do_pu = options.overwrite or not os.path.isfile(out_pileup_path)
    do_call = do_pu or not os.path.isfile(out_sample_vg_path)
    do_aug = do_pu or not os.path.isfile(out_augmented_vg_path)

    if do_pu:
        RealTimeLogger.get().info("Computing Variants for {} {}".format(
            input_graph_path,
            input_gam))
        robust_makedirs(os.path.dirname(out_pileup_path))
        run("vg pileup {} {} -t {} > {}".format(input_graph_path,
                                                input_gam,
                                                options.vg_cores,
                                                out_pileup_path))

    if do_call:
        robust_makedirs(os.path.dirname(out_sample_vg_path))
        run("vg call {} {} -r 0.001 -d 50 -e 150 -s 25 -t {} | vg ids -c - | vg ids -s -  > {}".format(input_graph_path,
                                                                                                       out_pileup_path,
                                                                                                       options.vg_cores,
                                                                                                       out_sample_vg_path))

    if do_aug:
        robust_makedirs(os.path.dirname(out_augmented_vg_path))
        run("vg call {} {} -r 0.001 -d 50 -e 150 -s 25 -t {} -l | vg ids -c - | vg ids -s - > {}".format(input_graph_path,
                                                                                                         out_pileup_path,
                                                                                                         options.vg_cores,
                                                                                                         out_augmented_vg_path))

def call_variants(job, options):
    """ run everything (root toil job)
    """
    for input_gam in options.in_gams:
        job.addChildJobFn(compute_vg_variants, input_gam, options,
                          cores=options.vg_cores)
        if not options.vg_only:
            job.addChildJobFn(compute_linear_variants, input_gam, options,
                              cores=options.vg_cores)
    
    
def main(args):
    
    options = parse_args(args) 
    
    RealTimeLogger.start_master()

    filtered_gams = []
    skip_words = options.skip.split(",")
    for gam in options.in_gams:
        skip_gam = False
        for word in skip_words:
            if len(word) > 0 and word in gam:
                skip_gam = True
        if not skip_gam:
            filtered_gams.append(gam)
    options.in_gams = filtered_gams

    for gam in options.in_gams:
        if len(gam.split("/")) < 3 or os.path.splitext(gam)[1] != ".gam":
            raise RuntimeError("Input gam paths must be of the form "
                               ".../<alg>/<reads>/<filename>.gam")

    # Make a root job
    root_job = Job.wrapJobFn(call_variants, options,
        cores=1, memory="2G", disk=0)
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

