#!/usr/bin/env python2.7

"""
 Write a line of counts for each unique quality value in the vcf
 counts are: snps indels other, probably a way to do this in bcftools 
 but in rush at the moment
"""


import argparse, sys, os, os.path, random, subprocess, shutil, bisect, math

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("in_vcf", type=str,
                        help="Input vcf file (- for stdin)")
    parser.add_argument("--bed", type=str, default=None,
                        help="Input bed regions")
    parser.add_argument("--balance", type=str, default=None,
                        help="comma separated list of fn,fp,tp tables to balance [debugging]")

                        
    args = args[1:]
    options = parser.parse_args(args)
    return options

# note: can re-use code from vcfFilterQuality here if needed
def get_qual_from_line(line, options):
    # quality 
    return float(line.split("\t")[5])

def vcf_qual_stats(vcf_path, bed_path = None, ignore_keywords = []):
    """ count up snps indels and others for each quality value, and return 
a table with the cumulative results.  note quality is expected to be a number
in 5th column of vcf """

    cmd = "bcftools view {}".format(vcf_path)
    if bed_path is not None:
        assert vcf_path[-3:] == ".gz" # gotta be indexed
        cmd += " -R {}".format(bed_path)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=sys.stderr, bufsize=-1)
    output, _ = p.communicate()
    assert p.wait() == 0
    vcf_file = [line + "\n" for line in output.split("\n")]

    # map quality score --> (snp count, indel count, other count)
    counts = dict()
    
    for line in vcf_file:
        if len(line) > 1 and line[0] != "#" and all(x not in line for x in ignore_keywords):
            toks = line.split("\t")
            ref = toks[3]
            alts = toks[4].split(",")
            qual = get_qual_from_line(line, None)
            if qual not in counts:
                counts[qual] = [0, 0, 0]
            # makes most sense for single allelic , normalized vcfs
            if all(len(alt) == 1 and len(ref) == 1 for alt in alts):
                counts[qual][0] += 1
            elif all(len(alt) != len(ref) for alt in alts):
                counts[qual][1] += 1
            else:
                counts[qual][2] += 1
                    
    # make cumulative table qual, snps, indels, others
    qvals = sorted(counts.keys())
    qvals.reverse()
    table = []
    if len(qvals) > 0:
        table += [[qvals[0], counts[qvals[0]][0], counts[qvals[0]][1], counts[qvals[0]][2]]]
    for i in range(1, len(qvals)):
        table += [[qvals[i]] + [x + y for x, y in zip(table[i-1][1:], counts[qvals[i]])]]

    return table

def balance_tables(fn_table, fp_table, tp_table):
    """ need to make one table """

    # total truth (for inferring false negatives)
    # take sum of last lines of tp and fn
    total = tp_table[-1][1:] if len(tp_table) > 0 else [0, 0, 0]
    if len(fn_table) > 0:
        total = [x + y for x,y in zip(total, fn_table[-1][1:])]

    # all quality values
    quals = set(x[0] for x in fp_table + tp_table)

    # make sure we have entry for every quality, grabbing next value
    # if not
    def extend_table(table):
        table.reverse() # make bisect easier to write if increasing
        for qual in quals:
            pos = bisect.bisect_left(table, [qual, 0, 0, 0])
            if len(table) == 0:
                table.append([qual, 0, 0, 0])
            elif pos == len(table):
                table.append([qual, 0, 0, 0])
            elif qual != table[pos][0]:
                other_row = table[pos]
                table.insert(pos, [qual, other_row[1], other_row[2], other_row[3]])
        assert table == sorted(table)
        assert len(table) == len(quals)
        table.reverse()

    extend_table(fp_table)
    extend_table(tp_table)
 
    # recompute fn_table (becuase we want it in input-qualities)
    # so, we assume fn = total - tp
    del fn_table[:]
    for i in range(len(fp_table)):
        fp_row = fp_table[i]
        tp_row = tp_table[i]
        assert fp_row[0] == tp_row[0]
        fn_table.append([tp_row[0], total[0] - tp_row[1], total[1] - tp_row[2], total[2] - tp_row[3]])
    if len(fn_table) > 0:
        assert all(x >= 0 for x in fn_table[-1][1:])

def main(args):
    options = parse_args(args)

    if options.balance is not None:
        tables = options.balance.split(",")
        assert len(tables) == 3
        parsed_tables = []
        for t in tables:
            with open(t) as f:
                pt = []
                for line in f:
                    toks = line.split("\t")
                    pt.append([float(x) for x in toks])
                parsed_tables.append(pt)

        balance_tables(parsed_tables[0], parsed_tables[1], parsed_tables[2])
        for i in range(len(tables)):
            with open(tables[i] + ".balanced", "w") as f:
                for line in parsed_tables[i]:
                    f.write("\t".join([str(x) for x in line]) + "\n")
        return 0;

    table = vcf_qual_stats(options.in_vcf, options.bed)

    # dump to stdout
    for line in table:
        sys.stdout.write("{}\t{}\t{}\t{}\t\n".format(line[0], line[1], line[2], line[3]))
	 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
