#!/usr/bin/env bash

# make_graphs.sh: make vg graphs for the HGVM pilot regions by progressive alignment

set -e

rm -Rf hgvm.tar.gz
wget http://hgwdev.sdsc.edu/~anovak/hgvm/hgvm.tar.gz
tar -xvzf hgvm.tar.gz
# Skip cenx for now. We can't handle its 15k reads well.
rm -Rf hgvm/CENX
for REGION in `ls hgvm`
do
    # Set up a vg command
    VG_COMMAND=(vg msga -t 32 -k 32 -m 32 -B 640 -D -b ref -f "hgvm/${REGION}/ref.fa")
    
    for FASTA in `ls hgvm/${REGION} | grep -v ref.fa | grep '\.fa$'`
    do
        # Add in all the non-ref FASTAs
        VG_COMMAND=(${VG_COMMAND[@]} "-f" "hgvm/${REGION}/${FASTA}")
    done
    
    # Run the command
    echo ${VG_COMMAND[@]} ">${REGION}.vg"
    time ${VG_COMMAND[@]} >${REGION}.vg || echo "=======FAILED!======"
done
