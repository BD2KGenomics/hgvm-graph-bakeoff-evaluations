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
from threading import Timer
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
    parser.add_argument("--fa_path", type=str, default="data/altRegions",
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
    parser.add_argument("--g1kvcf_path", type=str, default="data/g1kvcf",
                        help="path to search for 1000 genomes vcf and sequences. expects "
                        "these to be in <g1kvcf_path>/BRCA1/BRCA1.vcf. etc. ")
    parser.add_argument("--platinum_path", type=str, default="data/platinum",
                        help="path to search for platinum genomes vcf. expects "
                        "these to be in <g1kvcf_path>/BRCA1/BRCA1.vcf. etc. ")    
    parser.add_argument("--chrom_fa_path", type=str, default="data/g1kvcf/chrom.fa",
                        help="fasta file with entire chromosome info for all regions")
    parser.add_argument("--call_opts", type=str, default="",
                        help="options to pass to vg call.  wrap in \"\"")
    parser.add_argument("--pileup_opts", type=str, default="",
                        help="options to pass to vg pileup. wrap in \"\"")
    parser.add_argument("--filter_opts", type=str, default="",
                        help="options to pass to vg filter. wrap in \"\"")
    parser.add_argument("--platinum_samples", type=str, default="NA12877,NA12878",
                        help="comma-separated list of sample names that have vcf"
                        " data in platinum folder")
    parser.add_argument("--timeout", type=int, default=sys.maxint,
                        help="timeout in seconds for long jobs (vg surject in this case)")
    parser.add_argument("--skipBaseline", action="store_true",
                        help="dont make vg sample graphs from g1k and platinum vcfs")
    parser.add_argument("--genotype", action="store_true",
                        help="use vg genotype instead of vg call")
    parser.add_argument("--surject", action="store_true",
                        help="attempt to make a surjected bam from each gam")

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

def gam_index_path(alignment_path, options, tag=""):
    """ get output gam name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}.gam.index".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def gam_filter_path(alignment_path, options, tag=""):
    """ get output gam name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}.gam.filter".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def sample_vg_path(alignment_path, options, tag=""):
    """ get vg call output name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_sample.vg".format(tag)
    return os.path.join(out_dir(alignment_path, options), name)

def sample_txt_path(alignment_path, options, tag=""):
    """ get vg call output name from input gam
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_call.tsv".format(tag)
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

def g1k_vcf_path(alignment_path, platinum, filter_indels, options, tag=""):
    """ get (compressed) vcf output of sample filtering of 1000 genomes
    """
    name = os.path.splitext(os.path.basename(alignment_path))[0]
    name += "{}_sample.vcf".format(tag)
    base = "g1kvcf" if platinum is False else "platvcf"
    if filter_indels is True:
        base += "_filter"
    outdir = os.path.join(options.out_dir,
                          alignment_region_tag(alignment_path, options),
                          base)
    return os.path.join(outdir, name)

def g1k_fa_path(alignment_path, platinum, filter_indels, options, tag=""):
    """ get fa output of sample filtering of 1000 genomes
    """
    vcf = g1k_vcf_path(alignment_path, platinum, filter_indels, options, tag)
    return os.path.splitext(vcf)[0] + ".fa"

def g1k_vg_path(alignment_path, platinum, filter_indels, options, tag=""):
    """ get vg path constructed from g1k filtered sample
    """
    vcf = g1k_vcf_path(alignment_path, platinum, filter_indels, options, tag)
    return os.path.splitext(vcf)[0] + ".vg"

def run(cmd, stdout = sys.stdout, stderr = sys.stderr, timeout_sec = sys.maxint,
        timeout_dep = None, fail_hard = False):
    """ run command in shell and barf if it doesn't work or times out 
    """
    RealTimeLogger.get().info("RUN: {}".format(cmd))

    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=stdout, stderr=stderr)
    
    def timeout_fail(proc, cmd):
        os.kill(proc.pid, signal.SIGKILL)
        proc.kill()
        # we often check to see if some output file exists before running
        # if we timeout, make sure the file exists so rerunning wont timeout again
        # at same place (unless overwrite explicitly desired)
        if timeout_dep is not None and os.path.exists(timeout_dep):
            os.system("rm -rf {}; echo timeout > {}".format(timeout_dep, timeout_dep))
        if fail_hard is True:
            raise RuntimeError("Command: {} timed out".format(cmd))
        else:
            RealTimeLogger.get().warning("Command: {} timed out".format(cmd))

    # based on a comment in http://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    timer = Timer(timeout_sec, timeout_fail, [proc, cmd])
    try:
        timer.start()
        stdout, stderr = proc.communicate()
    finally:
        timer.cancel()
    sts = proc.wait()
    
    if sts != 0:
        if fail_hard is True:
            raise RuntimeError("Command: %s exited with non-zero status %i" %
                               (cmd, sts))
        else:
            RealTimeLogger.get().warning("Command: {} exited with non-zero status {}".format(cmd, sts))
    return sts

def compute_linear_variants(job, input_gam, options):
    """ project to bam, then run samtools to call some variants
    """
    input_graph_path = graph_path(input_gam, options)
    input_index_path = index_path(input_graph_path, options)

    # can only do this if there is a "ref" path in the vg graph
    res_path = temp_path(options)
    run("scripts/vgHasPath.sh {} {} > {}".format(input_graph_path, "ref", res_path))
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
                surject_path),
                timeout_sec=options.timeout,
                timeout_dep=surject_path)
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
    out_sample_vcf_path = out_sample_vg_path.replace(".vg", ".vcf")
    out_sample_txt_path = sample_txt_path(input_gam, options)    
    out_augmented_vg_path = augmented_vg_path(input_gam, options)
    out_gam_filter_path = gam_filter_path(input_gam, options)
    out_gam_index_path = gam_index_path(input_gam, options)
    out_bam_path = out_sample_vg_path.replace(".vg", ".bam")
    do_genotype = options.genotype and (options.overwrite or not os.path.isfile(out_sample_vcf_path))
    do_gam_filter= (options.overwrite or not os.path.isfile(out_gam_filter_path))
    do_gam_index = do_genotype and (do_gam_filter or options.overwrite or not os.path.isdir(out_gam_index_path))
    do_pu = not options.genotype and (options.overwrite or not os.path.isfile(out_pileup_path))
    do_call = not options.genotype and (do_pu or not os.path.isfile(out_sample_vcf_path))
    do_surject = options.surject and (options.overwrite or do_gam_filter or not os.path.isfile(out_bam_path))

    # We need an XG here for the mase graph, but I haven't got time to refactor
    # to make it in the right place. So just make it here.
    temp_xg_path = job.fileStore.getLocalTempDir() + "/filter.xg"

    if do_gam_filter or do_pu:
        # Make sure we have the xg index around, which fiulter may need.
        run("vg index -x {} {}".format(temp_xg_path, input_graph_path), fail_hard = True)

    if do_gam_filter:
        robust_makedirs(os.path.dirname(out_pileup_path))
        run("vg filter -x {} {} {} {} > {}".format(temp_xg_path, input_gam,
                                             options.filter_opts, input_graph_path,
                                             out_gam_filter_path),
            fail_hard = True)

    if do_gam_index:
        robust_makedirs(os.path.dirname(out_pileup_path))
        run("rm -rf {} ; vg index {} -N -d {}".format(out_gam_index_path, out_gam_filter_path, out_gam_index_path),
            fail_hard = True)        

    if do_pu:
        robust_makedirs(os.path.dirname(out_pileup_path))
        run("vg pileup {} {} {} -t {} > {}".format(input_graph_path,
                                                   out_gam_filter_path,
                                                   options.pileup_opts,
                                                   options.vg_cores,
                                                   out_pileup_path),
            fail_hard = True)
    ref = None
    bedLength = -1
    if do_call or do_genotype or do_surject:
        robust_makedirs(os.path.dirname(out_sample_vcf_path))
        region = alignment_region_tag(input_gam, options)
        g1kbed_path = os.path.join(options.g1kvcf_path, region.upper() + ".bed")            
        with open(g1kbed_path) as f:
            contig, offset, end = f.readline().split()[0:3]
            bedLength = int(end) - int(offset)
            
        # make the vcf
        # can only do this if there is a "ref" path in the vg graph
        ref = None
        res_path = temp_path(options)
        for ref_name in ["ref", contig]:
            run("scripts/vgHasPath.sh {} {} > {}".format(input_graph_path, ref_name, res_path))
            with open(res_path) as res_file:
                if res_file.read()[0] == "1":
                    ref = ref_name
                    break
        run("rm {}".format(res_path))
                
    if ref is not None:
        if do_genotype:
            run("vg genotype {} {} -S -pv -q -i -C -o {} -r {} -c {} -s {} -t {} > {} 2> {}".format(input_graph_path,
                                                                                                 out_gam_index_path,
                                                                                                 offset,
                                                                                                 ref,
                                                                                                 contig,
                                                                                                 alignment_sample_tag(input_gam, options),
                                                                                                 options.vg_cores,
                                                                                                 out_sample_vcf_path,
                                                                                                 out_sample_vcf_path.replace(".vcf", ".vcf.stderr")),
                fail_hard = True)
        if do_call:
            run("vg call {} {} {} -t {} -o {} -r {} -c {} -S {} -A {} > {} 2> {}".format(input_graph_path,
                                                                                         out_pileup_path,
                                                                                         options.call_opts,
                                                                                         options.vg_cores,
                                                                                         offset,
                                                                                         ref,
                                                                                         contig,
                                                                                         alignment_sample_tag(input_gam, options),
                                                                                         out_augmented_vg_path,
                                                                                         out_sample_vg_path.replace(".vg", ".vcf"),
                                                                                         out_sample_vg_path.replace(".vg", ".vcf.stderr")),
                fail_hard = True)
        if do_surject:
            run("vg index {} -k {} -e {} -s -d {}.index -t {}".format(input_graph_path, 20, 5,
                                                                         os.path.join(os.path.dirname(out_bam_path), "graph"),
                                                                         options.vg_cores),
                fail_hard = True)
            run("vg surject {} -t {} -p {} -b -d {}.index > {}".format(out_gam_filter_path, options.vg_cores, ref,
                                                                       os.path.join(os.path.dirname(out_bam_path), "graph"),
                                                                       out_bam_path),
                fail_hard = True)
            
            # fix up chromosome coordinates so we can display on browser
            if contig[0] != "c":
                contig = "chr{}".format(contig)
            contigLength = {"chr5": 181538259, "chr6": 170805979, "chr13": 114364328,
                            "chr17": 83257441, "chr19": 58617616}
            # in header, change up the contig name and 
            run("samtools view -H {} | sed -e \"s/{}/{}/\" | sed -e \"s/{}/{}/\" > {}.sam".format(out_bam_path,
                                                                                                  ref, contig,
                                                                                                  bedLength, contigLength[contig],
                                                                                                  out_bam_path),
                fail_hard = True)
            # in body, add offset and fix contig, leave in sam for now so we can debug
            run("samtools view -F 256 {} | awk -v OFS=\'\\t\' \'{{$3=\"{}\"; $4=$4+{}; $5=60; $8=$8+{}; print $0}}\' >> {}.sam".format(out_bam_path,
                                                                                                                                     contig,
                                                                                                                                     offset, offset,
                                                                                                                                     out_bam_path),
                fail_hard = True)
            # back to bam
            run("samtools view {}.sam -b -F 4 | samtools sort - --threads {} -o {}".format(out_bam_path, options.vg_cores, out_bam_path),
                fail_hard = True)

            # and index
            run("samtools index -b {}".format(out_bam_path), fail_hard = True)
            
                                                                                                                                 
                                                                              


def compute_snp1000g_baseline(job, input_gam, platinum, filter_indels, options):
    """ make 1000 genomes sample graph by filtering the vcf
    """
    # there is only one g1vcf graph per region per sample
    # this function is also going to get called once for each graph type
    # so we hack here to only run on refonly graphs (arbitrary choice)
    if alignment_graph_tag(input_gam, options) != "refonly":
        return

    sample = alignment_sample_tag(input_gam, options)

    if platinum is True and sample not in options.platinum_samples.split(","):
        return
    
    region = alignment_region_tag(input_gam, options)
    if platinum is False:
        g1kvcf_path = os.path.join(options.g1kvcf_path, region.upper() + ".vcf")
    else:
        g1kvcf_path = os.path.join(options.platinum_path, sample, region.upper() + ".vcf")
    g1kbed_path = os.path.join(options.g1kvcf_path, region.upper() + ".bed")
    filter_vcf_path = g1k_vcf_path(input_gam, platinum, filter_indels, options)
    filter_fa_path = g1k_fa_path(input_gam, platinum, filter_indels, options)
    filter_vg_path = g1k_vg_path(input_gam, platinum, filter_indels, options)
    fasta_path = options.chrom_fa_path

    do_filter = options.overwrite or not os.path.isfile(filter_vcf_path + ".gz")
    do_construct = do_filter or not os.path.isfile(filter_vg_path)

    # make sure we're dealing with a sample that's in the vcf
    if do_filter or do_construct:
        p = subprocess.Popen("grep {} {} | wc -l".format(sample, g1kvcf_path),
                             shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=-1)
        output, _ = p.communicate()
        assert p.wait() == 0
        if int(output) == 0:
            do_filter = False
            do_construct = False

    # make filtered compressed vcf for this sample
    if do_filter:            
        robust_makedirs(os.path.dirname(filter_vcf_path))
        if filter_indels is True:
            filter_input_path = filter_vcf_path + ".in"
            run("scripts/vcfFilterIndels.py {} > {}".format(g1kvcf_path, filter_input_path),
                fail_hard = True)
        else:
            filter_input_path = g1kvcf_path
        run("scripts/vcfFilterSample.py {} {} {} {} {}".format(filter_input_path,
                                                              fasta_path,
                                                              sample,
                                                              filter_vcf_path,
                                                              filter_fa_path),
            fail_hard = True)
        run("scripts/vcfsort {} > {}.sort ; mv {}.sort {}".format(filter_vcf_path,
                                                          filter_vcf_path,
                                                          filter_vcf_path,
                                                          filter_vcf_path))
        run("bgzip -f {}".format(filter_vcf_path), fail_hard = True)
        run("tabix -f -p vcf {}.gz".format(filter_vcf_path), fail_hard = True)

    # load it into a vg graph
    if do_construct:
        with open(g1kbed_path) as bed_file:
            coords = bed_file.readline().split()
            # convert from bed to vcf coordinates by adding one to start
            coords = (coords[0], int(coords[1]) + 1, int(coords[2]))
            run("vg construct -v {}.gz -r {} -t {} -R {}:{}-{} > {}".format(filter_vcf_path, filter_fa_path,
                                                                            options.vg_cores, 
                                                                            coords[0], coords[1], coords[2],
                                                                            filter_vg_path),
                fail_hard = True)    
        
def call_variants(job, options):
    """ run everything (root toil job)
    """
    for input_gam in options.in_gams:
        job.addChildJobFn(compute_vg_variants, input_gam, options,
                          cores=options.vg_cores)
        if not options.skipBaseline:
            # all combinations of [G1KVCF , PLATVCF] 
            job.addChildJobFn(compute_snp1000g_baseline, input_gam, False, False, options,
                              cores=options.vg_cores)
            job.addChildJobFn(compute_snp1000g_baseline, input_gam, True, False, options,
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
                             cores=1, memory="2G", disk="2G")
    
    # Run it and see how many jobs fail
    failed_jobs = Job.Runner.startToil(root_job,  options)
    
    if failed_jobs > 0:
        raise Exception("{} jobs failed!".format(failed_jobs))
                               
    RealTimeLogger.stop_master()
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
        
        
        
        
        
        
        
        
        
        

