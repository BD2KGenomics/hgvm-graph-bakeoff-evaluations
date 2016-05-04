#!/usr/bin/env python2.7
"""
Mapping high-coverage reads produces a directory structure with
/alignments (gam read alignments, on per sample per graph)
/indexes (indexed vg graphs, compressed). 
This script just crawls the indexes directory, unpacking and 
renaming all the graphs therein, such that the output 
graphs directory can be used by the variant analysis scripts.
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("in_files", nargs="+",
                        help="input archives. ex: high_coverage_alignments/indexes/*/*.tar.gz")
    parser.add_argument("out_dir", type=str, default="variants",
                        help="output directory. ex: high_coverage_graphs")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="overwrite existing files")
    parser.add_argument("--keep_index", action="store_true",
                        help="keep the .index for each graph")
    args = args[1:]
        
    return parser.parse_args(args)

def region_tag(path):
    """ given high_coverage_alignments/indexes/brca1/cactus.tar.gz, return brca1
    """
    return path.split("/")[-2]

def graph_tag(path):
    """ given high_coverage_alignments/indexes/brca1/cactus.tar.gz, return cactus
    """
    return os.path.basename(path).replace(".tar.gz", "")

def out_name(path):
    """ given high_coverage_alignments/indexes/brca1/cactus.tar.gz, return cactus-brca1
    """
    region = region_tag(path)
    graph = graph_tag(path)
    suf = ""
    if graph.find("debruijn") == 0:
        suf = graph[8:]
        graph = "debruijn"
    return "{}-{}{}".format(graph, region, suf)
        

def main(args):
    
    options = parse_args(args)
    out_dir = options.out_dir

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    assert os.path.isdir(out_dir)

    for in_file in options.in_files:
        dest = out_name(in_file)
        if options.overwrite or not os.path.isfile("{}/{}.vg".format(out_dir, dest)):
            # extract into output dir
            os.system("tar zxf {} -C {}".format(in_file, out_dir))
            # expect graph.vg and graph.vg.index
            # rename graph.vg
            os.system("mv {}/graph.vg {}/{}.vg".format(out_dir, out_dir, dest))
            if options.keep_index:
                # rename graph.vg.index
                os.system("mv {}/graph.vg.index {}/{}.vg.index".format(out_dir, out_dir, dest))
            else:
                # delete graph.vg.index 
                os.system("rm -rf {}/graph.vg.index".format(out_dir))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
