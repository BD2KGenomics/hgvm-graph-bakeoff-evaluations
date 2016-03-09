#!/usr/bin/env bash 
# mean.sh: take the mean of a column of numbers. Taken from
# <http://www.unix.com/302525444-post2.html?s=40551cead014050bc01956806694ba32>
set -e

awk '{x+=$1}
END {print x/NR}'
