#!/usr/bin/env python2.7
"""
Computes concordance between trio of sample graphs, using text output from vg call -c.  This concordance is the proportion of base calls in the child that can be found in at least one parent.  This will only work if the calls were all computed on the same reference graph. prints num_concordant num_discordant pct_concordant

"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit, string

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("child", type=str,
                        help="child calls (from vg call -c)")    
    parser.add_argument("parent1", type=str,
                        help="parent1 calls (from vg call -c)")
    parser.add_argument("parent2", type=str,
                        help="parent2 calls (from vg call -c)")
    
    args = args[1:]
        
    return parser.parse_args(args)

def parse_call(call):
    """ take a comman-separated call, and return lower case tuple
    """
    return map(lambda x : x.lower(), call.split(","))

def is_base(value):
    """ is value a call, ie not a -, which we treat as null
    """
    return value in ["a", "c", "g", "t", "."]

def score_call(child_call, dad_call, mom_call):
    """ return tuple of (concordant, discordant) calls in child
    """
    score = [0, 0]
    def check_call(value, score):
        if is_base(value):
            if value in dad_call or value in mom_call:
                score[0] += 1
            else:
                score[1] += 1
    check_call(child_call[0], score)
    check_call(child_call[1], score)

    return score        
    
def main(args):
    
    options = parse_args(args) 

    tot_right = 0
    tot_wrong = 0

    c_file = open(options.child)
    m_file = open(options.parent1)
    d_file = open(options.parent2)

    c_line = c_file.readline()
    while c_line:
        m_line = m_file.readline()
        d_line = d_file.readline()

        c_toks = c_line.split()
        m_toks = m_line.split()
        d_toks = d_line.split()

        if (len(c_toks) != 4 or len(m_toks) != 4 or len(d_toks) != 4 or
            c_toks[0] != m_toks[0] or c_toks[1] != m_toks[1] or c_toks[2] != m_toks[2] or
            c_toks[0] != d_toks[0] or c_toks[1] != d_toks[1] or c_toks[2] != m_toks[2]):
            raise RuntimeError("Calls cannot be compared because they are invalid " +
                                "or from different positions:" +
                                "\n  child: {}\n  parent1: {}\n  parent2: {}\n".format(c_line,
                                                                                       m_line,
                                                                                       d_line))

        c_call = parse_call(c_toks[3])
        m_call = parse_call(m_toks[3])
        d_call = parse_call(d_toks[3])

        score = score_call(c_call, m_call, d_call)
        tot_right += score[0]
        tot_wrong += score[1]
        
        c_line = c_file.readline()

    c_file.close()
    d_file.close()
    m_file.close()

    print "{}\t{}\t{}".format(tot_right, tot_wrong,
                              float(tot_right) / max(1, (tot_right + tot_wrong)))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
