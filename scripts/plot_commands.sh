#!/usr/bin/env bash

OUTDIR=paperplots

mkdir -p ${OUTDIR}

# Set up the plot parameters
# Include both versions of the 1kg SNPs graph name
PLOT_PARAMS=(
    --categories
    vglr
    trivial
    debruijn-k63
    snp1kg
    camel
    cactus
    prg
    sbg
    refonly
    shifted1kg
    haplo1kg30
    haplo1kg50
    debruijn-k31
    simons
    curoverse
    snp1000g
    level1
    level2
    level3
    --category_labels 
    VGLR
    Unmerged
    'De Bruijn 63'
    1KG
    Camel
    Cactus
    PRG
    7BG
    Primary
    Control
    '1KG Haplo 30'
    '1KG Haplo 50'
    'De Bruijn 31'
    SGDP
    Curoverse
    1KG
    Level1
    Level2
    Level3
    --colors
    '#cab2d6'
    '#b1b300'
    '#ff7f00'
    '#fb9a99'
    '#33a02c'
    '#1f78b4'
    '#6a3d9a'
    '#b15928'
    '#000000'
    '#FF0000'
    '#00FF00'
    '#0000FF'
    '#e31a1c'
    '#b2df8a'
    '#a6cee3'
    '#fb9a99'
    '#FF0000'
    '#00FF00'
    '#0000FF'
)

PLOT_PARAMS_MHC=(
    --categories
    trivial
    prg
    debruijn-k63
    cactus
    snp1kg
    haplo1kg50
    haplo1kg30
    simons
    refonly
    shifted1kg
    debruijn-k31
    sbg
    camel
    --category_labels 
    Unmerged
    PRG
    'De Bruijn 63'
    Cactus
    1KG
    '1KG Haplo 50'
    '1KG Haplo 30'
    SGDP
    Primary
    Control
    'De Bruijn 31'
    7BG
    Camel
    --colors
    '#b1b300'
    '#6a3d9a'
    '#ff7f00'
    '#1f78b4'
    '#fb9a99'
    '#0000FF'
    '#00FF00'
    '#b2df8a'
    '#000000'
    '#FF0000'
    '#e31a1c'
    '#b15928'
    '#33a02c'
)


./scripts/boxplot.py ./low_coverage_alignments//plots/absolute/indelrate.brca2.tsv --title 'Indels per base
in BRCA2 (absolute)' --x_label Graph --y_label 'Indel rate' --save ${OUTDIR}/absolute-indelrate.brca2.pdf --x_sideways --hline_median refonly --best_low --range --sparse_axes --max 0.004 --min 0.002 "${PLOT_PARAMS[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/absolute/indelrate.mhc.tsv --title 'Indels per base
in MHC (absolute)' --x_label Graph --y_label 'Indel rate' --save ${OUTDIR}/absolute-indelrate.mhc.pdf --x_sideways --hline_median refonly --best_low --range --sparse_axes --max 0.015 "${PLOT_PARAMS_MHC[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

scripts/scatter.py ./low_coverage_alignments//plots/absolute/perfect_vs_unique.brca2.tsv --save ${OUTDIR}/absolute-perfect_vs_unique.brca2.pdf --title 'Perfect vs. Unique
Mapping in BRCA2' --x_label 'Portion Uniquely Mapped' --y_label 'Portion Perfectly Mapped' --width 12 --height 9 --sparse_axes --markers o --min_x 0.835 --min_y 0.74 --max_x 0.930 --max_y 0.80 "${PLOT_PARAMS[@]}" --font_size 20 --dpi 90 --no_n --no_legend

scripts/scatter.py ./low_coverage_alignments//plots/absolute/perfect_vs_unique.mhc.tsv --save ${OUTDIR}/absolute-perfect_vs_unique.mhc.pdf --title 'Perfect vs. Unique
Mapping in MHC' --x_label 'Portion Uniquely Mapped' --y_label 'Portion Perfectly Mapped' --width 12 --height 9 --sparse_axes --markers o --min_x 0.80 --min_y 0.60 --max_x 0.9 --max_y 0.75 "${PLOT_PARAMS_MHC[@]}" --font_size 20 --dpi 90 --no_n --no_legend

./scripts/boxplot.py ./low_coverage_alignments//plots/absolute/substrate.brca2.tsv --title 'Substitution rate
in BRCA2 (absolute)' --x_label Graph --y_label 'Substitution rate' --save ${OUTDIR}/absolute-substrate.brca2.pdf --x_sideways --hline_median refonly --max 0.015 --min 0.009 --best_low --range --sparse_axes "${PLOT_PARAMS[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/absolute/substrate.mhc.tsv --title 'Substitution rate
in MHC (absolute)' --x_label Graph --y_label 'Substitution rate' --save ${OUTDIR}/absolute-substrate.mhc.pdf --x_sideways --hline_median refonly --max 0.05 --best_low --range --sparse_axes "${PLOT_PARAMS_MHC[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/normalized/indelrate.brca2.tsv --title 'Indels per base
in BRCA2 (normalized)' --x_label Graph --y_label 'Indel  relative rate' --save ${OUTDIR}/normalized-indelrate.brca2.pdf --x_sideways --hline_median refonly --best_low --range --sparse_axes --max 1.5 --min 0.5 "${PLOT_PARAMS[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/normalized/indelrate.mhc.tsv --title 'Indels per base
in MHC (normalized)' --x_label Graph --y_label 'Indel  relative rate' --save ${OUTDIR}/normalized-indelrate.mhc.pdf --x_sideways --hline_median refonly --best_low --range --sparse_axes --max 1.5 --min 0.5 "${PLOT_PARAMS_MHC[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/normalized/substrate.brca2.tsv --title 'Substitution rate
in BRCA2 (normalized)' --x_label Graph --y_label 'Substitution  relative rate' --save ${OUTDIR}/normalized-substrate.brca2.pdf --x_sideways --hline_median refonly --max 2 --min 0.5 --best_low --range --sparse_axes "${PLOT_PARAMS[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks

./scripts/boxplot.py ./low_coverage_alignments//plots/normalized/substrate.mhc.tsv --title 'Substitution rate
in MHC (normalized)' --x_label Graph --y_label 'Substitution  relative rate' --save ${OUTDIR}/normalized-substrate.mhc.pdf --x_sideways --hline_median refonly --max 2 --min 0.5 --best_low --range --sparse_axes "${PLOT_PARAMS_MHC[@]}" --font_size 20 --dpi 90 --no_n --medians_only --hide_categories refonly --hline_ticks
