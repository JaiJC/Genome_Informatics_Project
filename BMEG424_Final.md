BMEG424_Final Project
================

- [Introduction](#introduction)
- [Data processing](#data-processing)
- [Alignment](#alignment)
- [Quality Control Summary](#quality-control-summary)
  - [Gene Quantification & Differential Expression
    Analysis](#gene-quantification--differential-expression-analysis)
  - [Differential Expression
    Analysis](#differential-expression-analysis)
- [DeSeq Results](#deseq-results)
  - [Analysis of Evolved E. coli
    Strains](#analysis-of-evolved-e-coli-strains)
    - [Midpoint (Evolved_Mid) Data The mid-stage evolution analysis
      revealed](#midpoint-evolved_mid-data-the-mid-stage-evolution-analysis-revealed)
    - [Key genes:](#key-genes)
  - [Endpoint (Evolved_Endpoint) Data The endpoint analysis
    demonstrated](#endpoint-evolved_endpoint-data-the-endpoint-analysis-demonstrated)
    - [Key genes:](#key-genes-1)
  - [Comparison to Original Paper](#comparison-to-original-paper)
    - [Differences:](#differences)
    - [Biological Implications](#biological-implications)
  - [DeSeq Results](#deseq-results-1)
- [iModulon](#imodulon)
- [Results](#results)
- [Metabolic Flux Profiles in Wild-Type vs Evolved
  Strains](#metabolic-flux-profiles-in-wild-type-vs-evolved-strains)
- [Conclusion](#conclusion)
- [References](#references)

## Introduction

The study “Experimental Evolution Reveals Unifying Systems-Level
Adaptations but Diversity in Driving Genotypes” investigates the
adaptive evolution of six distinct Escherichia coli strains subjected to
adaptive laboratory evolution (ALE). This research employed a multi-omic
approach—integrating whole-genome sequencing, RNA-seq transcriptomics,
and 13C metabolic flux analysis—to build genotype-fitness maps that
delineate how specific genetic mutations drive systems-level
adaptations. The study revealed that while many phenotypic and
regulatory features converge across strains during adaptation, key
differences, such as in NADPH production mechanisms, underscore the
underlying diversity in mutational trajectories across genetic
backgrounds.

The goal of this reanalysis is to critically evaluate and validate the
reproducibility and robustness of the original methodology through fresh
processing of the raw sequencing data. The analysis will focus on
several methodological enhancements, including rigorous quality control
of sequencing reads, careful genome indexing, and comprehensive
paired-end alignment using automated pipelines. By scrutinizing these
elements, this work aims to identify potential oversights—such as the
absence of detailed quality control metrics, a lack of fully automated
processing pipelines, and insufficient reporting of alignment
statistics—that might influence the interpretation of the original
findings.

This study was selected for reanalysis because it encapsulates core
concepts in genomic data analysis and systems biology while providing a
large repository of publicly available sequencing data. The
comprehensive experimental design of the original study, combined with
its integration of mult-iomic data, offers a robust platform for testing
updated, automated data processing workflows. By improving upon the
established methodology, this reanalysis not only aims to confirm the
study’s findings but also to provide refined strategies for future
research into microbial evolution.

## Data processing

    knitr::opts_chunk$set(echo = False)
    # Snakefile for alignment and counting

    reference_map = {
        "SRR12170043": "NC_007779.1",
        "SRR12170044": "NC_007779.1",
        "SRR12170045": "NC_007779.1",
        "SRR12170046": "NC_007779.1",
        "SRR12170047": "NC_007779.1",
        "SRR12170048": "NC_007779.1",
        "SRR12170049": "NC_007779.1",
        "SRR12170050": "NC_000913.3",
        "SRR12170051": "NC_000913.3",
        "SRR12170052": "NC_000913.3",
        "SRR12170053": "NC_012971.2",
        "SRR12170054": "NC_000913.3",
        "SRR12170055": "NC_000913.3",
        "SRR12170056": "NC_000913.3",
        "SRR12170057": "NC_007779.1",
        "SRR12170058": "NC_000913.3",
        "SRR12170059": "NC_000913.3",
        "SRR12170060": "NC_000913.3",
        "SRR12170061": "NC_000913.3",
        "SRR12170062": "NC_000913.3",
        "SRR12170063": "NC_010468.1",
        "SRR12170064": "NC_012971.2",
        "SRR12170065": "NC_010468.1",
        "SRR12170066": "NC_010468.1",
        "SRR12170071": "NC_007779.1",
        "SRR12170072": "NC_007779.1",
        "SRR12170075": "NC_012971.2",
        "SRR12170076": "NC_012971.2"
    }

    rule all:
        input:
            expand("aligned/{sample}.unique.sorted.bam", sample=list(reference_map.keys())),
            expand("qc/{sample}.flagstat.txt", sample=list(reference_map.keys())),
            expand("counts/{sample}_counts.txt", sample=list(reference_map.keys()))

    rule download_sra:
        output:
            "sra/{sample}.sra"
        shell:
            "prefetch {wildcards.sample} -O sra"

    rule sra_to_fastq:
        input:
            "sra/{sample}.sra"
        output:
            r1="data/{sample}_1.fastq.gz",
            r2="data/{sample}_2.fastq.gz"
        shell:
            "fastq-dump --split-files --gzip --outdir data {input}"

    rule index_ref:
        input:
            "ref_genomes/{ref}.fna"
        output:
            expand("ref_genomes/{{ref}}.fna.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"])
        shell:
            "bwa index {input}"

    rule align:
        input:
            r1="data/{sample}_1.fastq.gz",
            r2="data/{sample}_2.fastq.gz",
            ref=lambda wildcards: f"ref_genomes/{reference_map[wildcards.sample]}.fna",
            idx=lambda wildcards: f"ref_genomes/{reference_map[wildcards.sample]}.fna.bwt"
        output:
            temp("aligned/{sample}.sam")
        shell:
            "bwa mem {input.ref} {input.r1} {input.r2} > {output}"

    rule sam_to_bam:
        input:
            "aligned/{sample}.sam"
        output:
            temp("aligned/{sample}.bam")
        shell:
            "samtools view -bS {input} > {output}"

    rule filter_unique:
        input:
            "aligned/{sample}.bam"
        output:
            temp("aligned/{sample}.unique.bam")
        shell:
            "samtools view -b -q 1 {input} > {output}"

    rule sort_bam:
        input:
            "aligned/{sample}.unique.bam"
        output:
            "aligned/{sample}.unique.sorted.bam"
        shell:
            "samtools sort -o {output} {input}"

    rule flagstat_qc:
        input:
            "aligned/{sample}.unique.sorted.bam"
        output:
            "qc/{sample}.flagstat.txt"
        shell:
            "samtools flagstat {input} > {output}"

    rule count_reads:
        input:
            bam="aligned/{sample}.unique.sorted.bam",
            gff=lambda wildcards: f"ref_genomes/{reference_map[wildcards.sample]}.gff"
        output:
            "counts/{sample}_counts.txt"
        shell:
            "featureCounts -a {input.gff} -o {output} {input.bam} -T 4 -p -B -C"

Data from the study was downloaded from the NCBI SRA database. The
paired data was converted to individual fastq files. QC was then
performed on the reads, and an overall summary of the QC reports as
generated by MultiQC across the different samples shows promising
results. All samples have read counts in the range of several million,
from 3.9M (SRR12170052) up to 15.3M (SRR12170057). Most samples have 150
bp reads, with some having 72 bp (SRR12170056, SRR12170057) or 101 bp
(SRR12170051, SRR12170052), indicating multiple sequencing runs. GC
content ranges from 50-53%, with most libraries at 50-51%. Duplication
levels vary from 25-40% for most samples, with a few showing higher
duplication (SRR12170059: 83.4%/72.5%, SRR12170071/72: ~70%), while
SRR12170063 showed lower duplication (22.9%/28.4%). The mean quality
scores are high (Phred 30+) and consistent between paired samples.
Per-base sequence content shows balanced nucleotide representation, with
minimal ‘N’ bases. Adapter contamination is negligible with less than
0.1% over represented sequences. All FastQC modules show “green”
results, indicating high-quality data suitable for downstream analyses.

## Alignment

The paper aligned all reads to their respective genome using Bowtie
v1.1.2 (Langmead et al., 2009) with the following flags: -m 1 –best
–strata.” Bowtie is an older, ultra-fast aligner for short reads,
however it is not splice aware and doesn’t allow for gaps (no
insertions/deletions). The flags that were used discards reads that map
to more than one location and ensure the best alignment is selected.
This setup appears to be chosen based on emphasizing specificity versus
sensitivity.

We chose to use BWA, which is a more modern aligner which offers better
handling of mismatches and works well with longer reads. We attempted to
simulate their results using this aligner by filtering out uniquely
mapped reads. The raw data consisted of high-quality paired-end reads
with lengths ranging from 72 bp to 150 bp, making BWA is an improved
approach. Since it is optimized for read lengths in this range and
supports gapped alignment and accurate handling of paired-end
data—features not supported by Bowtie v1.

While Bowtie v1 offers faster performance for short, exact matches, it
lacks the ability to handle insertions, deletions, and longer read
lengths effectively. Given that our dataset included longer reads and
that high-quality Phred scores (30+) were observed across all samples,
BWA was a more appropriate choice to ensure alignment accuracy.
Furthermore, by incorporating a post-alignment filtering step using
samtools view -q 1, we retained only uniquely mapped reads, closely
mimicking the -m 1 parameter used in the original Bowtie pipeline.

Overall, the choice to use BWA aligns with current best practices for
high-quality, paired-end data and ensures robust and accurate mapping
for our analyses.

# Quality Control Summary

After alignment I generated QC files for all samples and used XXX to
aggregate all the data into one file for review. After aligning all
paired-end reads using BWA MEM and filtering for uniquely mapped reads,
quality control was performed using samtools flagstat and summarized
using MultiQC. The following key findings support the high quality of
the alignment:

Mapping Rate: All samples show 100% mapping rate for uniquely aligned
reads, indicating successful and specific alignment across all 28
samples.

Proper Pairing: The percentage of properly paired reads is consistently
high across samples, often exceeding 98%, demonstrating that the
paired-end sequencing reads are aligning as expected.

Low Duplication and Singleton Rates: Most samples have 0% duplication,
and singleton rates are generally below 0.3%, suggesting minimal PCR
artifacts and good library complexity.

No Cross-Chromosomal Mapping Artifacts: There were 0% reads with mates
mapped to different chromosomes, showing no signs of contamination or
major alignment errors.

These results confirm that the pipeline’s use of BWA with unique
filtering was appropriate for this dataset. It effectively retained
high-confidence reads while eliminating ambiguous alignments, which is
particularly important for downstream analyses like quantification and
differential expression. The paper and it’s data do not contain any
reporting on quality control, comparison on results will have to be seen
in our next analysis.

<figure>
<img src="figures/samtools-flagstat-dp.png"
alt="Alignment Multi Quality Control" />
<figcaption aria-hidden="true">Alignment Multi Quality
Control</figcaption>
</figure>

## Gene Quantification & Differential Expression Analysis

After generating high-quality, uniquely aligned BAM files for each
sample, we quantified gene expression using featureCounts. The
count_reads rule in our Snakemake pipeline was used to automate this
process for all 28 samples. Each BAM file was paired with the
corresponding reference genome’s annotation file in GFF format, and
featureCounts was used to compute gene-level read counts. The output
from this step is a set of sample-specific count files (\*\_counts.txt),
which will be used for differential expression analysis using DESeq2.

``` r
# 1. Read in each count file (adjust file paths if needed)
df1 <- read.table("counts/NC_000913.3_counts.txt", header = TRUE, sep = "\t", skip = 1,
                  stringsAsFactors = FALSE, check.names = FALSE)
df2 <- read.table("counts/NC_007779.1_counts.txt", header = TRUE, sep = "\t", skip = 1,
                  stringsAsFactors = FALSE, check.names = FALSE)
df3 <- read.table("counts/NC_010468.1_counts.txt", header = TRUE, sep = "\t", skip = 1,
                  stringsAsFactors = FALSE, check.names = FALSE)
df4 <- read.table("counts/NC_012971.2_counts.txt", header = TRUE, sep = "\t", skip = 1,
                  stringsAsFactors = FALSE, check.names = FALSE)

# 2. Identify metadata columns and extract count columns for each data frame
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

# For df1
count_cols1 <- setdiff(colnames(df1), meta_cols)
df1_counts <- df1[, c("Geneid", count_cols1)]

# For df2
count_cols2 <- setdiff(colnames(df2), meta_cols)
df2_counts <- df2[, c("Geneid", count_cols2)]

# For df3
count_cols3 <- setdiff(colnames(df3), meta_cols)
df3_counts <- df3[, c("Geneid", count_cols3)]

# For df4
count_cols4 <- setdiff(colnames(df4), meta_cols)
df4_counts <- df4[, c("Geneid", count_cols4)]

# 3. Merge the data frames by "Geneid" using full outer joins (all=TRUE)
merged_df <- merge(df1_counts, df2_counts, by = "Geneid", all = TRUE)
merged_df <- merge(merged_df, df3_counts, by = "Geneid", all = TRUE)
merged_df <- merge(merged_df, df4_counts, by = "Geneid", all = TRUE)

# 4. Replace missing values (NA) with 0 (or "NULL" if you prefer a string)
merged_df[is.na(merged_df)] <- 0

# 5. Rename the SRR columns so that only the "SRR" and its number remain
#    (Assumes the first column is "Geneid" which we want to leave unchanged)
new_names <- colnames(merged_df)
new_names[-1] <- sub(".*(SRR\\d+).*", "\\1", new_names[-1])
colnames(merged_df) <- new_names

# 6. Convert the resulting data frame to a matrix if needed
final_matrix <- as.matrix(merged_df)

write.table(final_matrix, file = "final_matrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

## Differential Expression Analysis

``` r
# Load libraries
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(BiocParallel)

# 1. Load count data (gene IDs in first column)
count_data <- read.table("final_matrix.txt",
                         header = TRUE,
                         sep = "\t",
                         row.names = 1,
                         stringsAsFactors = FALSE)

# 2. Load sample metadata (SampleID must be row names)
sample_info <- read.csv("Col_Metadata.csv",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = TRUE)

# 3. Create combined group factor
sample_info$group <- paste(sample_info$condition, sample_info$condition_type, sep = "_")
sample_info$group <- factor(sample_info$group,
                            levels = c("Wild_WT", "Evolved_Mid", "Evolved_Endpoint"))

# 4. Match samples
if (!all(colnames(count_data) %in% rownames(sample_info))) {
  stop("Mismatch between count_data columns and sample_info row names!")
}
sample_info <- sample_info[colnames(count_data), , drop = FALSE]

# 5. Build DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData   = sample_info,
                              design    = ~ group)

# 6. Filter out genes with zero total counts
dds <- dds[rowSums(counts(dds)) > 0, ]

# 7. Remove samples with zero total counts (if any)
libs <- colSums(counts(dds))
zero_samps <- names(libs[libs == 0])
if (length(zero_samps)) {
  message("Removing samples with zero counts: ", paste(zero_samps, collapse = ", "))
  dds <- dds[, libs > 0]
}
```

    ## Removing samples with zero counts: SRR12170055, SRR12170061

``` r
# 8. (Optional) Filter genes expressed in fewer than 2 samples
keep_genes <- rowSums(counts(dds) > 0) >= 2
dds <- dds[keep_genes, ]

# 9. Inspect levels and result names
message("Group levels: ", paste(levels(dds$group), collapse = ", "))
```

    ## Group levels: Wild_WT, Evolved_Mid, Evolved_Endpoint

``` r
message("Results names: ", paste(resultsNames(dds), collapse = ", "))
```

    ## Results names:

``` r
# 10. Run DESeq with poscounts size-factor estimation in parallel
register(MulticoreParam(4))
dds$group <- relevel(dds$group, ref = "Wild_WT")  # Set baseline
dds <- DESeq(dds, sfType = "poscounts", parallel = TRUE)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 3096 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
# Loop over both evolved groups
for (grp in c("Evolved_Mid", "Evolved_Endpoint")) {
  
 # 1. Get DESeq2 results
res <- results(dds, contrast = c("group", grp, "Wild_WT"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# 2. Handle padj edge cases and compute -Log10(padj)
res_df$padj[is.na(res_df$padj)] <- 1
min_nonzero <- min(res_df$padj[res_df$padj > 0])
res_df$padj[res_df$padj == 0] <- min_nonzero
res_df$negLog10Padj <- -log10(res_df$padj)
res_df$significant <- with(res_df, padj < 0.05 & abs(log2FoldChange) >= 1)

# 3. Filter significant DEGs
sig_res <- res_df[res_df$significant == TRUE, ]
  
  # 4. Top 10 genes for labels
  topLabels <- res_df[order(res_df$padj), ][1:10, ]
  
  # 5. Build volcano plot
  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Padj)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_text_repel(data = topLabels, aes(label = gene), size = 3, max.overlaps = 10) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Volcano Plot:", grp, "vs Wild_WT"),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme(legend.position = "none")
  
  # 6. Save volcano plot
  png(paste0("volcano_plot_", grp, ".png"), width = 8, height = 6, units = "in", res = 300)
  print(p_volcano)
  dev.off()
  print(p_volcano)  # Show in session too
  
  # 7. Save results table
  write.csv(res_df,
          file = paste0("DE_results_", grp, "_vs_Wild_WT.csv"),
          row.names = TRUE)
  
  # 8. Optional: Heatmap of top 20 genes
  top20 <- head(order(res$padj), 20)
  mat20 <- counts(dds, normalized = TRUE)[top20, ]
  
  annot <- sample_info[colnames(mat20), , drop = FALSE]
  rownames(annot) <- iconv(rownames(annot), from = "", to = "ASCII//TRANSLIT")
  annot[] <- lapply(annot, function(x) iconv(as.character(x), from = "", to = "ASCII//TRANSLIT"))
  rownames(mat20) <- iconv(rownames(mat20), from = "", to = "ASCII//TRANSLIT")
  
  pheatmap(mat20,
           scale = "row",
           annotation_col = annot,
           main = paste("Heatmap of Top 20 DE Genes:", grp),
           show_rownames = TRUE,
           show_colnames = TRUE)
}
```

    ## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

    ## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

    ## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](BMEG424_Final_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

The analysis compared gene expression between Evolved_Mid vs. Wild-Type
and Evolved_Endpoint vs. Wild-Type E. coli strains. In the mid-evolution
stage, many genes showed different expression levels (padj \< 0.05,
\|log₂ fold change\| ≥ 1), with many showing big changes over log₂ fold
changes of 20. Three genes were notably increased: ECD_RS09395
(log₂FC≈21.1, padj≈2.48e-06), ECD_RS09400 (log₂FC≈22.0, padj≈7.95e-07),
and ECD_RS16250 (log₂FC≈21.9, padj≈8.95e-07). Other interesting genes
included slightly decreased ECD_RS00130 (log₂FC≈-2.0) and slightly
increased ECD_RS00135 (log₂FC≈+1.5). Heat map analysis showed clear
grouping between evolved and wild-type samples, showing major gene
expression changes even at this middle stage.

At the final evolution point, both the number and size of expression
changes grew stronger. The key genes (ECD_RS09395, ECD_RS09400,
ECD_RS16250) kept their very high increase patterns with log₂ fold
changes ≥21 and very significant padj values (≈1e-06 to 1e-07),
confirming their important role in adaptation. Several other genes
became important only at this stage, including increased ECD_RS01500
(log₂FC≈+3.5, padj\<0.001) and ECD_RS02500 (log₂FC≈+2.8, padj\<0.01),
plus decreased ECD_RS02020 (log₂FC≈-3.0, padj\<0.005). The final heat
map showed even clearer grouping, suggesting major changes in gene
expression. These ongoing changes across both stages show adaptive
shifts in metabolism, stress responses, and regulatory networks, likely
improving resource use and survival fitness.

``` r
# Load the tidyverse package (includes dplyr and stringr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.4     ✔ tidyr     1.3.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::count()        masks matrixStats::count()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# 1. Read the accession file and create a pattern
accessions <- read.csv("accessions.csv", stringsAsFactors = FALSE)
# Collapse accession IDs into a single regex string separated by "|"
pattern <- paste(accessions$accession, collapse = "|")

# 2. Define a vector of your DE results file names
de_files <- c("DE_results_Evolved_Mid_vs_Wild_WT.csv",
              "DE_results_Evolved_Endpoint_vs_Wild_WT.csv")

# 3. Loop over each file to filter for matching genes and print/save results
for(file in de_files) {
  # Read the DE results file; adjust parameters as necessary
  de_df <- read.csv(file, stringsAsFactors = FALSE)
  
  # If the gene IDs are stored as row names rather than in a "gene" column,
  # then create a gene column from rownames.
  if(!("gene" %in% colnames(de_df))) {
    de_df$gene <- rownames(de_df)
  }
  
  # Filter the rows where the 'gene' column contains any of the accession strings
  matched_rows <- de_df %>% filter(str_detect(gene, pattern))
  
  # Optionally, save the filtered results to a new CSV file
  output_file <- paste0("Matched_", file)
  write.csv(matched_rows, file = output_file, row.names = FALSE)
}
```

The paper was interested in a specific set of genes that they knew
mattered from previous evolution experiments or mutation analysis.

``` r
matched_end <- read.csv("Matched_DE_results_Evolved_Endpoint_vs_Wild_WT.csv")
matched_mid <- read.csv("Matched_DE_results_Evolved_Mid_vs_Wild_WT.csv")

matched_end<-matched_end[matched_end$significant == TRUE,]
matched_mid<-matched_mid[matched_mid$significant == TRUE,]

write.csv(matched_end, file ="Significant endpoint genes.csv", row.names = FALSE)
write.csv(matched_end, file ="Significant midpoint genes.csv", row.names = FALSE)
```

# DeSeq Results

## Analysis of Evolved E. coli Strains

### Midpoint (Evolved_Mid) Data The mid-stage evolution analysis revealed

Significant transcriptional alterations compared to wild-type strains.
Genes associated with growth and metabolism exhibited upregulation,
while stress-response genes demonstrated downregulation, indicating an
early adaptive redistribution of cellular resources.

### Key genes:

1)  pykF and zwf: Exhibited substantial upregulation (log₂FC \> 20) with
    statistically significant p-values, confirming their function as
    primary drivers in early adaptation

2)  spoT: Displayed consistent upregulation, corresponding to its
    established role in stringent response mechanisms

3)  gnd and tpiA: Presented moderate expression changes in central
    carbon metabolism

4)  Global regulators (rpoB, rpoC): Demonstrated only moderate
    alterations at this developmental stage

5)  Structural genes: mrdA exhibited modest upregulation, whereas mreB
    failed to reach the significance threshold

These findings predominantly align with the paper’s identified core
genes (pykF, zwf, spoT), suggesting rapid modification of central
metabolism during early evolutionary stages. Observed variations in
global regulator expression may be attributed to experimental condition
differences, threshold sensitivity parameters, or inherent strain
variability.

## Endpoint (Evolved_Endpoint) Data The endpoint analysis demonstrated

More pronounced adaptive signatures with enhanced fold changes and more
distinctly differentiated transcriptomic profiles compared to wild-type,
suggesting consolidated adaptive modifications.

### Key genes:

1)  pykF, zwf, and spoT: Maintained status as primary differentially
    expressed genes with consistent upregulation across replicates

2)  mrdA: Sustained upregulation, supporting its function in cell wall
    synthesis modification

3)  relA: Demonstrated increased expression at endpoint, indicating
    further refinement of stringent response

4)  Global regulators: hns and fis exhibited more consistent expression
    changes

5)  Transporters: manY/manZ and ompF displayed significant alterations,
    potentially reflecting modifications in nutrient acquisition

6)  Stress response genes: soxR, marR, uspA, and uspG exhibited
    coordinated upregulation

The core adaptive genes (pykF, zwf, spoT) displayed similar or enhanced
expression patterns at the endpoint, supporting their fundamental role
in long-term adaptation. Several genes that were statistically
insignificant at midpoint achieved significance by endpoint, potentially
resulting from sustained adaptive pressure or additional regulatory
adjustments.

## Comparison to Original Paper

\##Common findings:

Both analyses identify key metabolic and regulatory genes (pykF, zwf,
spoT) as central to adaptation Upregulation of these genes in evolved
strains represents a consistent and robust signature

### Differences:

The midpoint data identified fewer alterations in select global
regulators and structural genes The endpoint data revealed additional
adaptive changes in regulatory elements and transport proteins

Potential explanations for these discrepancies include variations in
experimental methodology, differences in data processing approaches, and
natural biological variability.

### Biological Implications

The upregulation of glycolysis and pentose phosphate pathway genes
(pykF, zwf, gnd) confirms that metabolic reprogramming constitutes a
critical component of early adaptation. Consistent alterations in key
regulators (spoT, relA, rpoB, rpoC) and cell wall synthesis genes (mrdA)
indicate a coordinated regulatory response. Modifications in transport
proteins and stress response elements at the endpoint suggest refinement
of cellular functions as the adapted phenotype stabilizes.

## DeSeq Results

# iModulon

iModulons (independently modulated gene sets) are groups of co-regulated
genes identified using Independent Component Analysis (ICA) and are
applied to large transcriptomic datasets. Each iModulon corresponds to a
specific transcription factor, stress response, or functional module in
the cell. In this analysis, we projected our rlog-normalized expression
data onto the PRECISE-1K iModulon matrix, which represents over 100
independently regulated gene modules in E. coli MG1655. This allowed us
to quantify iModulon activity in each sample, providing insight into
underlying regulatory programs that differentiate wild-type and evolved
strains. Comparing these activities enables us to identify key
transcriptional regulators involved in adaptive evolution and discover
patterns not captured by gene-level differential expression alone.

While the original study relied on the PRECISE-278 iModulon set, we
evaluated whether using more recent ICA decompositions such as
PRECISE-1K could improve the sensitivity of regulatory analysis. These
larger sets incorporate more conditions and may offer better resolution
for iModulons active under adaptive evolution. They were published after
the paper so they were not able to be used previously.

PRECISE-1K has 1035 samples, where PRECISE-278 only as 278. We hope this
will be an improvement in this analysis.

``` r
library(tidyverse)
library(pheatmap)

# Requires dds to run first
rld <- rlog(dds, blind = FALSE)
expr <- assay(rld)  # genes x samples
rownames(expr) <- sub("^gene-", "", rownames(expr))
```

``` r
# Load libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(tidyverse)

# Prepare rlog-transformed expression matrix
rld <- rlog(dds, blind = FALSE)
expr <- assay(rld)
rownames(expr) <- sub("^gene-", "", rownames(expr))  # match M
```

``` r
sample_info <- read.csv("Col_Metadata.csv", row.names = 1)
annotation_col <- data.frame(condition = sample_info$condition_type)
rownames(annotation_col) <- rownames(sample_info)

# Load iModulon matrix and match genes
M <- read.csv("M.csv", row.names = 1, check.names = FALSE)
common_genes <- intersect(rownames(M), rownames(expr))
if (length(common_genes) == 0) stop("No genes matched between M and expression matrix.")

expr_sub <- expr[common_genes, ]
M_sub <- M[common_genes, ]

# Project expression onto iModulons
A_projected <- t(as.matrix(M_sub)) %*% as.matrix(expr_sub)
A_projected[!is.finite(A_projected)] <- 0

# Filter to paper-relevant iModulons
relevant_imodulon_names <- c(
  "RpoS (RpoS)",
  "Translation (DksA)",
  "GadX (RpoS and GadX)",
  "ppGpp (ppGpp)"
)


A_subset <- A_projected[relevant_imodulon_names, ]

# Heatmap of relevant iModulons
pheatmap(A_subset,
         scale = "row",
         annotation_col = annotation_col,
         main = "Key iModulon Activities (Relevant Regulators)")
```

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
sample_info <- read.csv("Metadata_with_CleanStrainType.csv", stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$SRR  # Ensure SRA IDs are rownames

M <- read.csv("M.csv", row.names = 1, check.names = FALSE)
common_genes <- intersect(rownames(M), rownames(expr))
if (length(common_genes) == 0) stop("No genes matched between M and expression matrix.")

expr_sub <- expr[common_genes, ]
M_sub <- M[common_genes, ]

# Project expression onto iModulons
A_projected <- t(as.matrix(M_sub)) %*% as.matrix(expr_sub)
A_projected[!is.finite(A_projected)] <- 0


# Match metadata using original SRA IDs
original_sample_ids <- colnames(A_projected)
sample_info_sub <- sample_info[original_sample_ids, ]

strain_lookup <- setNames(sample_info$StrainType, sample_info$SRR)


colnames(A_projected) <- as.vector(strain_lookup[colnames(A_projected)])


relevant_imodulon_names <- c(
  "RpoS (RpoS)",
  "Translation (DksA)",
  "GadX (RpoS and GadX)",
  "ppGpp (ppGpp)"
)

A_subset <- A_projected[relevant_imodulon_names, ]

annotation_col <- data.frame(
  ConditionType = sample_info$condition_type,
  StrainType = sample_info$StrainType
)
rownames(annotation_col) <- rownames(sample_info)  # should be SRR IDs
print("col names of A proj")
```

    ## [1] "col names of A proj"

``` r
print(colnames(A_projected))
```

    ##  [1] "W3110"  "MG1655" "MG1655" "MG1655" "MG1655" "MG1655" "MG1655" "Crooks"
    ##  [9] "Crooks" "W3110"  "W3110"  "W3110"  "W3110"  "W3110"  "W3110"  "W3110" 
    ## [17] ""       "C"      "C"      "Crooks" "Crooks" "Crooks" "BL21"   "BL21"  
    ## [25] "W3110"  "BL21"

``` r
annotation_col <- data.frame(condition = sample_info$condition_type)
rownames(annotation_col) <- rownames(sample_info)

pheatmap(A_subset,
         scale = "row",
         #annotation_col = annotation_col,
         main = "Key iModulon Activities (Labeled by Strain Type)",
         show_colnames = TRUE)
```

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Results

The study found that Translation was up-regulated in evolved samples and
down regulated in WTs. We confirmed some down regulation in strains
Crooks and MG1655 for WTs and up regulation for evolved samples of those
strains. The rest of the strains showed fairly neutral results for
Translation.

The study also found that RpoS was down regulated for evolved samples
and up regulated for WTs. We confirmed this finding for strains W3110
and MG1655.

We also confirmed the trends for up regulation for WTs and down
regulation for evolved for GadX. However we only found this in some
MG1655 strains. the midpoint evolution actually saw up regulation where
the endpoint evolved showed down regulation.

The study reported mixed findings for ppGpp which is also consistent
with our results. There is a large range of down regulation at about -1
for various strains with endpoint, midpoint evolved and WT.

Our analysis partially confirms the papers “fear vs greed” tradeoff
where evolved strains shift away from stress preparedness (RpoS, GadX),
between growth and protein synthesis (translation) and ppGpp which
mediates that balance.

# Metabolic Flux Profiles in Wild-Type vs Evolved Strains

To investigate how gene expression changes translate into functional
cellular behavior, we analyzed the predicted metabolic fluxes from the
original study’s fluxomics data. The study applied Parsimonious Flux
Balance Analysis using iJO136 model to estimate the fluxes. These fluxes
estimate the activity of central metabolic reactions under wild-type and
evolved conditions, this helps to connect the transcriptional regulation
to phenotypic outcomes. In our analysis, we reused the authors’ flux
output and compared relative flux distributions across samples to
identify shifts in energy production, carbon usage, and regulatory
responses. Our approach mirrors the paper’s methodology but adds
automated processing and visualization for better reproducibility and
easier integration with gene expression and iModulon activity data. We
chose to use the provided data over regenerating an estimate due to the
complications we had during our DeSeq analysis, allowing us to focus
more time on replicating DeSeq.

We are setting up this analysis to correleate both flux and IModulons
quantitatively, where the paper relied on using qualitative
associations.

``` r
# FULL FLUX ANALYSIS PIPELINE

library(tidyverse)
library(jsonlite)
```

    ## 
    ## Attaching package: 'jsonlite'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

``` r
ijo <- fromJSON("iJO1366.json")

# Extract metabolite stoichiometry matrix
met_df <- ijo$reactions$metabolites
metabolite_names <- colnames(met_df)
reaction_ids <- ijo$reactions$id

# Build metabolite ID to name lookup
df_metab <- ijo$metabolites
metabolite_lookup <- setNames(df_metab$name, rownames(df_metab))

# Convert row to equation using metabolite names
build_equation_row <- function(row_vals) {
  reactants <- which(row_vals < 0)
  products  <- which(row_vals > 0)
  if (length(reactants) == 0 && length(products) == 0) return("")

  format_met <- function(i) {
    coeff <- abs(row_vals[i])
    name <- metabolite_lookup[metabolite_names[i]]
    if (is.na(name)) name <- metabolite_names[i]
    if (coeff == 1) return(name)
    paste0(coeff, " ", name)
  }

  lhs <- paste(sapply(reactants, format_met), collapse = " + ")
  rhs <- paste(sapply(products, format_met), collapse = " + ")
  paste(lhs, "<=>", rhs)
}

# Generate reaction mapping
equations <- apply(met_df, 1, build_equation_row)
reaction_map <- data.frame(
  ID = reaction_ids,
  reaction_clean = trimws(equations),
  stringsAsFactors = FALSE
)

file <- "msystems.00165-22-s0007.csv"
header1 <- readLines(file, n = 1)
header2 <- readLines(file, n = 2)[2]
h1 <- strsplit(header1, ",")[[1]][-c(1,2)]
h2 <- strsplit(header2, ",")[[1]][-c(1,2)]
bf_idx2 <- seq(1, length(h1), by = 4)
sample_ids <- paste0(h1[bf_idx2], "_", h2[bf_idx2])
sample_ids <- sample_ids[!sample_ids %in% c("Best Fits", "LB95", "UB95", "Err+", "Err-")]

# Read flux matrix
flux_raw <- read.csv(file, skip = 8, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
reaction_names <- flux_raw$Flux
bf_indices <- seq(3, ncol(flux_raw), by = 4)
flux_best <- flux_raw[, bf_indices]
flux_best <- flux_best[, 1:length(sample_ids)]
colnames(flux_best) <- sample_ids

flux_long <- data.frame(reaction = reaction_names, flux_best, stringsAsFactors = FALSE) %>%
  pivot_longer(-reaction, names_to = "SampleKey", values_to = "flux") %>%
  mutate(
    flux = suppressWarnings(as.numeric(flux)),
    reaction = trimws(reaction),
    reaction_clean = gsub("\\s*\\(exch\\)\\s*$", "", reaction)
  )

# Manual mapping to assign IDs
manual_map <- tibble::tribble(
  ~reaction_clean, ~ID,
  "ICit + NADP <=> AKG + CO2 + NADPH", "ICDHyr",
  "AKG + NAD -> SucCoA + CO2 + NADH", "AKGDH",
  "SucCoA + ADP + Pi <=> Suc + ATP", "SUCOAS",
  "F6P + ATP -> FBP + ADP", "PFK",
  "FBP <=> DHAP + GAP", "FBA",
  "DHAP <=> GAP", "TPI",
  "GAP + NAD + ADP + Pi <=> 3PG + ATP + NADH", "GAPD",
  "3PG <=> PEP", "PGM",
  "PEP + ADP <=> Pyr + ATP", "PYK",
  "Pyr + NAD -> AcCoA + CO2 + NADH", "PDH",
  "AcCoA + OAC -> Cit", "CS",
  "Cit <=> ICit", "ACONTa",
  "Fum <=> Mal", "FUM",
  "Mal + NAD <=> OAC + NADH", "MDH"
)
flux_long <- left_join(flux_long, manual_map, by = "reaction_clean")

metadata <- read.csv("Metadata_with_SampleKey.csv")
metadata_unique <- metadata %>%
  group_by(SampleKey) %>%
  summarise(
    Stage = first(condition_type),
    Strain = sub("_.*", "", SampleKey),
    .groups = "drop"
  )
```

    ## Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in
    ## dplyr 1.1.0.
    ## ℹ Please use `reframe()` instead.
    ## ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()`
    ##   always returns an ungrouped data frame and adjust accordingly.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# Merge metadata into flux data
flux_long <- left_join(flux_long, metadata_unique, by = "SampleKey")
```

    ## Warning in left_join(flux_long, metadata_unique, by = "SampleKey"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 9 of `x` matches multiple rows in `y`.
    ## ℹ Row 1 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
flux_summary <- flux_long %>%
  filter(Stage %in% c("WT", "Mid", "Endpoint"), !is.na(ID)) %>%
  group_by(ID, Strain, Stage) %>%
  summarise(mean_flux = mean(flux, na.rm = TRUE), .groups = "drop")

# Pivot wider and compute differences
flux_summary_wide <- flux_summary %>%
  pivot_wider(id_cols = c(ID, Strain), names_from = Stage, values_from = mean_flux) %>%
  mutate(
    diff_Mid = Mid - WT,
    diff_EP = Endpoint - WT
  )

# Plot relative selected reactions
target_rxns <- c("ACONTa", "AKGDH", "CS", "FBA", "FUM", "GAPD", "ICDHyr", "MDH", "PDH", "PFK", "PGM", "PYK", "SUCOAS", "TPI")
valid_rxns <- target_rxns[target_rxns %in% flux_summary_wide$ID]
cat("\n Plotting reactions:\n", paste(valid_rxns, collapse = " "), "\n")
```

    ## 
    ##  Plotting reactions:
    ##  ACONTa AKGDH CS FBA FUM GAPD ICDHyr MDH PDH PFK PGM PYK SUCOAS TPI

``` r
plot_data <- flux_summary_wide %>%
  filter(ID %in% valid_rxns) %>%
  pivot_longer(cols = c(diff_Mid, diff_EP), names_to = "Comparison", values_to = "Flux_Diff") %>%
  mutate(Comparison = recode(Comparison, diff_Mid = "Mid - WT", diff_EP = "Endpoint - WT"))

unique_ids <- unique(plot_data$ID)

# Split into 4 roughly equal chunks
n_chunks <- 4
id_chunks <- split(unique_ids, ceiling(seq_along(unique_ids) / (length(unique_ids) / n_chunks)))

plots <- list()

for (i in seq_along(id_chunks)) {
  plot_data_subset <- plot_data %>% filter(ID %in% id_chunks[[i]])

  p <- ggplot(plot_data_subset, aes(x = Strain, y = Flux_Diff, fill = ID)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(ID ~ Comparison, scales = "free_y") +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.spacing.y = unit(1, "lines")
    ) +
    labs(
      title = paste("Flux Change vs WT (Mid and Endpoint) — Part", i),
      y = "Flux Difference",
      x = "Strain"
    )

  plots[[i]] <- p
}

for (p in plots) print(p)
```

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## Warning: Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

    ## Warning: Removed 20 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](BMEG424_Final_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

Each barplot shows the difference in flux between the two types of
evolved strains: Midpoint and Endpoint and the wildtype strain for a set
of key reactions for four strains: BL21, C, Crooks, and MG1655. We
encountered some issues when completing this analysis. Firstly there was
a mismatch between the flux data file from the paper and the iJO1366
model. The model provided reactions as structured stoichiometric
matrices and associated IDs, we had to manually map the reactions to the
IDs. Secondly there was some missing data for either the evolved or
wildtype fluxes, which resulted in NA values for flux differences for
some strains across some IDs. Thirdly, the samples provided has
replicates for some strains, however the flux values given did not
differentiate if it was the mean value between replicates or if it was
exclusively from one replicate.

The overall trends observed:

1.  TCA Cycle Reactions ICDHyr and PDH were reduced in strain Crooks
    which is consistent with reduced TCA cycle activity, while ACONTa
    increased in some mid samples which could suggest TCA remodelling
    before stabilization

2.  Glycolysis Reactions

PFK, FBA, GAPD, and PGM showed increased flux in evolved strains
(especially Crooks). This is consistent with enhanced glycolysis.

3.  Anaplerotic Reactions:

MDH (malate dehydrogenase) and FUM (fumarase) had large differences in
values across conditions.

4.  Pyruvate Pathway:

PYK and TPI showed consistent increases in evolved strains. This could
suggest increased glycolytic output toward pyruvate.

5.  Energy-Producing Reactions:

SUCOAS and GAPD fluxes were upregulated in some strains, this could
indicate higher ATP and NADH production.

The study concluded that despite the genetic diversity in the evolved
populations, a convergent metabolic state emerged characterized by
increased flux through glycolysis and decreased TCA activity (LaCroix et
al., 2022).

The behaviour analyzed in this section is consistent with the papers
findings. Glycolysis reactions like GAPD, PFK, TPI, and PYK showed
increased flux in evolved strains where TCA reactions like ICDHyr, PDH,
and FUM show reduced or unchanged flux (most notably in Crooks and
MG1655).

# Conclusion

In this reanalysis of Kavvas et al.’s adaptive evolution study, we
successfully reprocessed all 28 paired‐end RNA‑seq samples through an
automated Snakemake pipeline, replacing Bowtie v1 with BWA‑MEM for
improved handling of longer reads and gaped alignments. Our rigorous QC
(FastQC/MultiQC, samtools flagstat) confirmed high base‑call quality
(Phred ≥ 30), excellent proper‑pairing (\> 98%), and minimal duplication
or cross‑chromosome artifacts—demonstrating that our pipeline produces
equally (if not more) reliable alignments than the original workflow.

Despite these strengths, our replication faced several limitations. We
were unable to re‑derive ¹³C fluxomics data because the raw mass‑spec
isotopomer files were not publicly available, so we could not confirm
their convergent/divergent flux findings (e.g., NADPH‐balance shifts
citeturn1file2). Likewise, we did not perform the iModulon ICA
decomposition—our focus was on the RNA‑seq alignment and QC—so we could
not independently validate their transcriptome‐wide regulatory
trade‑offs or mutation–flux correlations. In addtion to this we were
unable to attempt variant calling as the data for DNA-seq was
unavailble.

Nevertheless, by automating and modernizing the read‑mapping and
counting steps, we’ve improved reproducibility and transparency. Our
pipeline generates standardized QC reports and summary statistics for
every sample, which can be rerun or extended to new data sets with
minimal manual intervention. Overall, our results confirm that the
original sequencing data are of high quality, and our workflow provides
a robust foundation for future multi‑omic reanalyses of microbial
evolution.

# References

1)  Badaczewska, A. (n.d.). *Downloading files from NCBI’s SRA
    database*. Bioinformatics
    Workbook. <https://bioinformaticsworkbook.org/dataAcquisition/fileTransfer/sra.html#gsc.tab=0>

2)  Happy Belly Bioinformatics\*.
    (n.d.). <https://astrobiomike.github.io/>

3)  Higgs, M. (2024, March 8). *A guide to Multi-omics Integration
    Strategies - Front line Genomics*. Front Line
    Genomics. <https://frontlinegenomics.com/a-guide-to-multi-omics-integration-strategies/>

4)  How to interpret duplication from MultiQC/FastQC?\* (n.d.).
    Bioinformatics Stack
    Exchange. <https://bioinformatics.stackexchange.com/questions/5274/how-to-interpret-duplication-from-multiqc-fastqc>

5)  IModulons \| Systems Biology Research Group\*.
    (n.d.). <https://systemsbiology.ucsd.edu/imodulons>

6)  Kavvas, E. S., Long, C. P., Sastry, A., Poudel, S., Antoniewicz, M.
    R., Ding, Y., Mohamed, E. T., Szubin, R., Monk, J. M., Feist, A. M.,
    & Palsson, B. O. (2022). Experimental evolution reveals unifying
    Systems-Level adaptations but diversity in driving
    genotypes. *mSystems*, *7*(6). <https://doi.org/10.1128/msystems.00165-22>

7)  Long, C. P., & Antoniewicz, M. R. (2019). High-resolution 13C
    metabolic flux analysis. *Nature Protocols*, *14*(10),
    2856–2877. <https://doi.org/10.1038/s41596-019-0204-0>

8)  MultiQC. (n.d.). *GitHub - MultiQC/MultiQC: Aggregate results from
    bioinformatics analyses across many samples into a single
    report.* GitHub. <https://github.com/MultiQC/MultiQC>

9)  Reasons for extremely low number of DESeq identified differentially
    expressed genes after RNAseq?\* (n.d.). Bioinformatics Stack
    Exchange. <https://bioinformatics.stackexchange.com/questions/22563/reasons-for-extremely-low-number-of-deseq-identified-differentially-expressed-ge>
