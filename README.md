# HGVM Graph Bakeoff Evaluations

## Background

This repository contains scripts to validate and analyse reference variation graphs submitted to the graph bakeoff [(Details)](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Map-%28HGVM%29-Pilot-Project).  They should be sufficient to reproduce all analysis results. 

## Obtaining Submitted Graphs

Todo: I think we have multiple ways of getting the graphs.  I've been piggybacking on graphs from Adam's scripts.  Sean may be using something else.  Worth moving into one place? 

## Dependencies

Submodules and/or Docker may be the way to go here?
*  [VG](https://github.com/ekg/vg) Note: vg/vcflib/bin should added to PATH too
*  [toil](https://github.com/BD2KGenomics/toil)
*  [samtools](https://github.com/samtools/samtools)
*  [bcftools](https://github.com/samtools/bcftools)
*  [htslib](https://github.com/samtools/htslib)
*  [corg](https://github.com/adamnovak/corg)
*  Python modules detailed in `requirements.txt`. Installable with `pip install -r requirements.txt`.


## Graph Properties Comparison

## Graph Mapping Comparison

One evaluation that we perform is to test how well each graph works as a read mapping target. This is accomplished by aligning the 1000 Genomes low coverage reads for each region against each graph, and evaluating the resulting alignments.

### Obtaining Low Coverage Reads

To download the low coverage reads for the read mapping evaluation, you can use the `getAltReads.py` script. Note that this script will by default, in single machine mode, download in parallel as many samples as your system has cores, and that overall it wants to download 6 regions for 2,691 samples. To mitigate the load on your system, you can use `--maxCores N` to limit the number of cores used (and thus the number of parallel downloads), and `--sampleLimit N` to limit the number of samples retrieved to the first N found. (Note that this sample will not be representative; it will be biased towards the populations that are iterated over first.)

An example command invocation:

```
# Remove any previous Toil job store
rm -Rf tree
# Run the download
time scripts/getAltReads.py ./tree low_coverage_reads --batchSystem=singleMachine --sample_ftp_root "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/" --population_pattern="*" --sample_pattern "*" --file_pattern "*.low_coverage.cram"  2>&1 | tee log.txt
```

### Running the Mapping Evaluation

To actually run the alignments for the downloaded reads, you should use the `parallelMappingEvaluation.py` script. An example invocation on a single machine is given below:

```
# Remove any previous Toil job store
rm -Rf tree
# Run the alignment
scripts/parallelMappingEvaluation.py ./tree graph_servers.tsv ./low_coverage_reads ./low_coverage_alignments --batchSystem=singleMachine --kmer_size=27 --edge_max=7 2>&1 | tee log.txt
```

Note that if you are running on all 2,691 samples retrieved in the previous step, the command given above will take an impractical amount of time to complete. One option to mitigate this is to only run some samples, with a `--sampleLimit` flag. Another is to use the script on a cluster, by changing the `--batchSystem` argument to any batch system supported by toil (and the job store of `./tree` to a job store that will be accessible from your cluster nodes).

Note that, if your cluster does not have a shared file system, the only storage backend currently supported for the input and output files is Microsoft Azure storage, using a syntax similar to that used for Toil job stores (`azure:<account>:<container>/<prefix>`).

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

The mapping steps above will result make a directory, `./high_coverage_alignments`.  It will contain a mapping for every region for every sample for every tool: `./high_coverage_reads/alignments/<region>/<tool>/sample.gam`  It will also contain tarballs of graphs and indexes used: `./high_coverage_reads/indexes/`.  We begin by extracting these into a new directory

     scripts/extractGraphs.py ./high_coverage_alignments/indexes/*/*.tar.gz ./high_coverage_graphs
	  
### Calling Variants

The variant calling scripts require the original FASTA inputs.  Uncompress them as follows:

     tar xzf data/altRegions.tar.gz

Fasta and vcf data is also required to make some baseline sample graphs for the 1000 Genomes data.  This is prepared with the following commands:

	  tar xzf data/g1kvcf.tar.gz
	  scripts/downloadChromFa.py

To run both the vg and samtools variant calling pipeline on all samples, use the following script

     ./scripts/call_snps.sh ./high_coverage_graphs ./high_coverage_alignments/alignments ./high_coverage_variants

This will create a structure in `./high_coverage_variants` similar to that in `alignments` (organized by tool and region) with calling output for all samples. **Note**: The number of threads used is hardcoded in this script and should be modified accordingly (see `--maxCores and --vgCores` options in OPTS variable). 

### Analysing Called Variants

To generate basic statistics of the calls, as well as do pairwise kmer comparisons of all output graphs, run

     ./scripts/compare_snps.sh ./high_coverage_graphs ./high_coverage_alignments/alignments ./high_coverage_variants ./variants_analysis

The output of this script is organized by region.  In each region directory, will be a variety of heatmaps comparing various output gaphs with each other, as well as some tables in `.tsv` format with some basic counting and size statistics. 


