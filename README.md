# HGVM Graph Bakeoff Evaluations

## Background

This repository contains scripts to validate and analyse reference variation graphs submitted to the graph bakeoff [(Details)](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Map-%28HGVM%29-Pilot-Project).  They should be sufficient to reproduce all analysis results. 

## Obtaining Submitted Graphs

Todo: I think we have multiple ways of getting the graphs.  I've been piggybacking on graphs from Adam's scripts.  Sean may be using something else.  Worth moving into one place? 

## Dependencies

Submodules and/or Docker may be the way to go here?
*  [VG](https://github.com/ekg/vg) 
*  [toil](https://github.com/BD2KGenomics/toil)
*  [samtools](https://github.com/samtools/samtools)
*  [bcftools](https://github.com/samtools/bcftools)
*  [htslib](https://github.com/samtools/htslib)

## Graph Properties Comparison

## Graph Mapping Comparison

## Graph Variant Calling Comparison

One evaluation that we perform is to test how well each graph works as a reference for calling variants against.

### Performing Initial Alignments

This evaluation depends on having alignments of 1000 Genomes Project high-coverage reads to each graph.

#### Obtaining High Coverage Reads

The alignments, in turn, are created from reads. Reads are downloaded from the 1000 Genomes FTP using the `scripts/getAltReads.py` script.

For example, to download all of the available high coverage samples to the `high_coverage_reads` directory, running only on the local machine, you can do this:

```
# Remove any previous Toil job store
rm -Rf tree
# Run the download
time scripts/getAltReads.py ./tree high_coverage_reads --batchSystem=singleMachine --sample_ftp_root "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/" --population_pattern="*" --sample_pattern "*" --file_pattern "*.high_coverage.cram"  2>&1 | tee log.txt
```

#### Aligning High Coverage Reads

Currently the easiest way to produce the alignments needed for the variant calling evaluation is to run the read alignment evaluation script on the high coverage samples. If you used the command above to download the reads, this can be accomplished thusly:

```
# Remove any previous Toil job store
rm -Rf tree
# Run the alignment
scripts/parallelMappingEvaluation.py ./tree graph_servers.tsv ./high_coverage_reads ./high_coverage_alignments --batchSystem=singleMachine --kmer_size=27 --edge_max=7 2>&1 | tee log.txt
```

It might be wise to parallelize this step across a cluster; Microsoft Azure Storage is supported for the input and output directories, using a syntax similar to that used for Toil job stores (`azure:<account>:<container>/<prefix>`).

### Expected Alignment Location

The mapping comparison must be run before the variant calling comparison.  In particular the following directory structure it outputs is required (top-level names and locations are generally flexible)

* `graphs/` location containing all original graphs in `.vg` format as well as their rocksdb indexes.  Graph names must be of the form `<tool>-<region>.vg`
* `alignments/<region>/<tool>/sample.gam` directory structure of vg alignments. 
  
The names in these two directories must be consistent.  For example, `alignments/brca1/camel/NA19239.gam` requires the existence of `graphs/camel-brca1.vg` and `graphs/camel-brca1.vg.index` for processing. 

### Calling Variants

To run both the vg and samtools variant calling pipeline on all samples, use the following script

     ./scripts/call_snps.sh <graphs_dir> <alignments_dir> <out_dir>

This will create a structure in `<out_dir>` similar to that in `alignments` (organized by tool and region) with calling output for all samples. **Note**: The number of threads used is hardcoded in this script and should be modified accordingly (see `--maxCores and --vgCores` options in OPTS variable). 

### Analysing Called Variants

To generate basic statistics of the calls, as well as do pairwise kmer comparisons of all output graphs, run

     ./scripts/compare_snps.sh <graphs_dir> <alignments_dir> <calls_dir> <out_dir>
     
The `<calls_dir>` argument here must be the same as the `<out_dir>` argument that was passed to `call_snps.sh`.  Like the latter `compare_snps.sh` has hardcoded threading options.  

The output of this script is organized by region.  In each region directory, will be a variety of heatmaps comparing various output gaphs with each other, as well as some tables in `.tsv` format with some basic counting and size statistics. 


