#!/usr/bin/env bash 
# sum.sh: take the sum of a column of numbers. Based on
# <http://www.unix.com/302525444-post2.html?s=40551cead014050bc01956806694ba32>
set -e

awk '{x+=$1}
END {print x}'
