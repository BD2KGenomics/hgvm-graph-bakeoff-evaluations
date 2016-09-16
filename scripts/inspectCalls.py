#!/usr/bin/env python2.7
"""
Given a graph and VCF, find the variants in the graph, and 
make little graphviz plots of them.  Can use to, ex, debug
false positive calls. 
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("vg_path", type=str,
                        help="input vg file")
    parser.add_argument("vcf_path", type=str,
                        help="input vcf file")
    parser.add_argument("region", type=str,
                        help="region (for computing offset) [BRCA1, BRCA2, MHC, LRC_KIR, SMA]")
    parser.add_argument("out_dir", type=str,
                        help="output dir to create")
    parser.add_argument("--context", type=int, default=1,
                        help="context expansion for vg find")
    parser.add_argument("--aug_context", type=int, default=10,
                        help="context for expanson for vg find on augmented graph")
    parser.add_argument("--step", type=int, default=1,
                        help="only consider every STEP variants in VCF")
    parser.add_argument("--ref", type=str, default=None,
                        help="use this vg path")
    parser.add_argument("--clip", type=str, default=None,
                        help="clip using this bed file")
    parser.add_argument("--qual", type=float, default=None,
                        help="ignore qualities less than this")
    parser.add_argument("--delta", type=str, default=None,
                        help="subtract this vcf using vcfDelta.py")
    parser.add_argument("--gam", type=str, default=None,
                        help="gam for showing mapped reads as paths")
    parser.add_argument("--filter_opts", type=str,
                        default="-r 0.90 -d 0.00 -e 0.00 -afu -s 10000 -q 15 -o 0 --defray-ends 999",
                        help="options for vg filter")
    parser.add_argument("--gam_view", action="store_true",
                        help="use vg view -A to view gams (instead of embedding as paths)")
    parser.add_argument("--snp_only", action="store_true",
                        help="only snps")
    parser.add_argument("--indel_only", action="store_true",
                        help="only indels")

    args = args[1:]
        
    return parser.parse_args(args)

def run(cmd, proc_stdout = sys.stdout, proc_stderr = sys.stderr,
        check = True):
    """ run command in shell and throw exception if it doesn't work 
    """
    print cmd
    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=proc_stdout, stderr=proc_stderr)
    output, errors = proc.communicate()
    sts = proc.wait()
    if check is True and sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (cmd, sts))
    return output, errors

def index_graph(vg_path, xg_path):
    """ make the xg """
    run("vg index {} -x {}".format(vg_path, xg_path))

def xg_path_node_id(xg_path, path_name, vg_offset):
    """ use vg find to get the node containing a given path position """
    #NOTE: vg find -p range offsets are 0-based inclusive.  
    stdout, stderr = run("vg find -x {} -p {}:{}-{} | vg mod -o - | vg view -j - | jq .node[0].id".format(
        xg_path, path_name, vg_offset, vg_offset),
                         proc_stdout=subprocess.PIPE)
    return int(stdout)                         

def get_offset(region):
    return {"BRCA1" : 43044293, "BRCA2" : 32314860,
            "LRC_KIR" : 54025633, "MHC" : 28510119,
            "SMA" : 69216818}[region.upper()]

def get_path_name(region):
    return {"BRCA1" : 17, "BRCA2" : 13,
            "LRC_KIR" : 19, "MHC" : 6,
            "SMA" : 5}[region.upper()]
    
def main(args):
    
    options = parse_args(args)

    if not os.path.isdir(options.out_dir):
        os.makedirs(options.out_dir)

    assert not options.snp_only or not options.indel_only

    # make the xg index in our output folder
    xg_path = os.path.join(options.out_dir, os.path.basename(options.vg_path) + ".xg")
    index_graph(options.vg_path, xg_path)

    # filter vcf into our output folder
    vcf_path = os.path.join(options.out_dir, os.path.basename(options.vcf_path))
    if options.delta is not None:
        tmp_path = vcf_path + ".pre_delta"
        run("scripts/vcfDelta.py {} {} -p > {}".format(options.vcf_path, options.delta, tmp_path))
        run("bgzip -f {}".format(tmp_path))
        run("tabix -f -p vcf {}.gz".format(tmp_path))
        input_vcf = tmp_path + ".gz"
    else:
        input_vcf = options.vcf_path
    cmd = "bcftools view {}".format(input_vcf)
    if options.clip is not None:
        cmd += " -R {}".format(options.clip)
    if options.snp_only is True:
        cmd += " -v snps,mnps"
    if options.indel_only is True:
        cmd += " -V snps,mnps"        
    if options.qual is not None:
        cmd += " | vcfFilterQuality.py - {}".format(options.qual)
    run(cmd + " > {}".format(vcf_path))

    # read vcf into memory
    vcf_lines = []
    with open(vcf_path) as vcf:
        for line in vcf:
            if line[0] != "#":
                vcf_lines.append(line)

    # chunk gam.  chunk i corresponds to ith vcf line
    if options.gam is not None:
        chunk_path = os.path.join(options.out_dir, "gam_chunks")
        bed_path = os.path.join(chunk_path, "coords.bed")
        chunk_base = os.path.join(chunk_path, "chunk")
        run("mkdir -p {}".format(chunk_path))
        with open(bed_path, "w") as f:
            for i, line in enumerate(vcf_lines[::options.step]):
                toks = line.split()
                vcf_pos = int(toks[1])
                # convert from vcf 1-based position to vg 0-based path offset
                vg_pos = vcf_pos - 1 - get_offset(options.region)
                path_name = get_path_name(options.region) if options.ref is None else options.ref
                f.write("{}\t{}\t{}\n".format(path_name, vg_pos, vg_pos + 1))
        run("vg filter {} -x {} -R {} -B {} {}".format(
            options.gam, xg_path, bed_path, chunk_base, options.filter_opts))

    for i, line in enumerate(vcf_lines[::options.step]):
        # parse vcf line
        toks = line.split()
        vcf_pos = int(toks[1])
        ref = toks[3]
        alt = toks[4]
        filtered=toks[6]
        
        if filtered != "PASS" and filtered != ".":
            # Skip filtered-out variants
            continue
        
        # convert from vcf 1-based position to vg 0-based path offset
        vg_pos = vcf_pos - 1 - get_offset(options.region)
        path_name = get_path_name(options.region) if options.ref is None else options.ref
        
        # figure out the node at this exact location
        ref_node = xg_path_node_id(xg_path, path_name, vg_pos)

        # put a bit of information into the file name
        record_name = "{}_{}_{}_rn_{}".format(vg_pos, ref, alt, ref_node)
        if "XREF" not in line:
            record_name += "_aug"
        
        # write the vcf line to a file
        with open(os.path.join(options.out_dir, record_name + ".vcf"), "w") as f:
            f.write(line)

        # chunk vg
        vg_chunk_path = os.path.join(options.out_dir, record_name + ".vg")
        run("vg find -p {}:{}-{} -x {} -c {} > {}".format(
            path_name, vg_pos, vg_pos, xg_path, options.context, vg_chunk_path))
                    
        if options.gam is not None:
            # make gam chunk
            gam_path = os.path.join(options.out_dir, record_name + ".gam")
            run("mv {}-{}.gam {}".format(chunk_base, i, gam_path))

            # default to embadding alignments as paths since vg view -A
            # is not very readable
            if options.gam_view is False:
                # make a big vg with gam paths
                big_path = os.path.join(options.out_dir, record_name + "_big.vg")
                run("vg mod {} -i {} > {}".format(options.vg_path, gam_path, big_path))

                # re-index
                big_index = os.path.join(options.out_dir, record_name + "_big.xg")
                run("vg index {} -x {}".format(big_path, big_index))

                # figure out the node at this exact location in the augmented index
                aug_ref_node = xg_path_node_id(big_index, path_name, vg_pos)
            
                # chunk augmented vg
                aug_chunk = os.path.join(options.out_dir, record_name + "_big.vg")
                run("vg find -p {}:{}-{} -x {} -c {} > {}".format(
                    path_name, vg_pos, vg_pos, big_index, options.aug_context, aug_chunk))
                    
                    
                # draw the gam-augmented chunk
                aug_png = os.path.join(options.out_dir, record_name.replace(
                    "rn_{}".format(ref_node), "rn_{}".format(aug_ref_node)) + "_gam.png")
                run("vg view -pd {} | dot -Tpng > {}".format(aug_chunk, aug_png))

        # activate gam view if desired
        if options.gam is not None and options.gam_view is True:
            view_opts = "-A {}".format(gam_path)
        else:
            view_opts = ""
            
        # draw the vg chunk
        png_path = os.path.join(options.out_dir, record_name + ".png")
        run("vg view {} -pd {} | dot -Tpng > {}".format(
            view_opts, vg_chunk_path, png_path))

        
        


    return 0
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
