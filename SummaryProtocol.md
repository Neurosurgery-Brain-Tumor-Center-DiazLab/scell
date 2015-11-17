# Single-Cell RNA-seq Analysis Protocols

This document is includes protocols for the pre-processing and analysis of single-cell RNA-seq libraries, including: alignment and read-level QC, quality filtering and outlier library detection, gene panel selection, normalization, dimentionality reduction, clustering and lineage reconstruction. The preprocessing scripts are meant to be run on a high-performance computing cluster running TORQUE.

## 1. Library Pre-processing

#### 1.1. Read Trimming and FastQC

We use the wrapper [**TrimGalore!**](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), which employs [**cutadapt**](https://wiki.gacrc.uga.edu/wiki/Cutadapt) to trim adapter sequences from reads after demultiplexing, as well as to trim low quality bases at the ends of reads.

To trim reads for leftover barcode sequences and low-quality bases, run the shell script [`sub_trim.sh`](https://dl.dropboxusercontent.com/u/9990581/SCell/PreProcessingScripts/sub_trim.sh) from the terminal, specifying the directory of single-cell libraries you want to trim. This usually corresponds to a sequencing plate of single-cell libraries.

```
sh sub_trim.sh sequencingRunDir
```

This script calls and executes the PBS script [`trim.pbs`](https://dl.dropboxusercontent.com/u/9990581/SCell/PreProcessingScripts/trim.pbs), which will submit the trimming jobs for each of the libraries in the specified directory as a job array.

A new subdirectory called `/trimmed` will be created under each individual library's directory, which now contains the quality-and-barcode-trimmed reads in a `fastq` file.

The `/trimmed` subdirectory will also contain the output from [**FastQC**](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for each individual library in your sequencing run. The MATLAB script `summarize_fastqc.m` can be used to summarize FastQC metrics over ensembles of individual cells.

#### 1.2. Read Alignment
The scripts provided here make use of [**TopHat 2**](https://ccb.jhu.edu/software/tophat/index.shtml), a splice junction short read aligner for RNA-Seq data.

To align your single-cell libraries, run the shell script [`sub_align.sh`]() from the terminal, specifying the directory of single-cell libraries you wish to align. This script calls and executes the PBS script [`align.pbs`](), which will submit the alignment jobs for each of the *trimmed* libraries in the specified directory as a job array.

```
sh sub_align.sh /path/to/sequencingRunDir
```

In the `sub_align.sh` script, you must specify:

 • The path to the Bowtie genome index directory containing the necessary `.bt2` files for alignment.

 • The path to the transcriptome index containing the necessary `.bt2` files to be used for the guided alignment.

 • The base name for the `.bt2` files in the Bowtie genome index directory to be used for alignment.

```
export bowtie_idx_path="/usr/refs/Homo_sapiens/NCBI/build37.2/Sequence/Bowtie2Index/"
```

```
export transcriptome_idx_path="/usr/refs/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes"
```
```
export genome_name="genome_name"
```

This will create a new directory called `/tophat_out`, which contains the `.bam` files and alignment summary statistics corresponding to each single-cell library.

Mapping rate statistics contained in the file `accepted_hits.stats` from each library can be summarized by running the script `comp_mapping_rates.m`

#### 1.3. Gene Expression Quantification

[**FeatureCounts**](http://bioinformatics.oxfordjournals.org/content/30/7/923.full.pdf?keytype=ref&ijkey=ZzPz96t2lqzAH6F), part of the [**SubRead**](http://subread.sourceforge.net/) package, assigns aligned reads/fragments to genomic features, such as genes.

The output of this read summarization process is a counts matrix, which contains the number of reads/fragments assigned to each feature for each one of your single-cell libraries.

To obtain a matrix of gene counts from your aligned single-cell libraries, run the shell script
`sub_featureCounts.sh` from the terminal, specifying the name of the single-cell library directory to summarize.

```
sh sub_featureCounts.sh sequencingRunDir
```
This script calls the PBS script `run_featureCounts.pbs`.

You will need a GTF file containing the features you want to summarize aligned reads to. Refer to [this example](https://dl.dropboxusercontent.com/u/9990581/genes_ercc.gtf), or supply your own annotation file.

In the `run_featureCounts.pbs` file, you must specify the path to the GTF file to use for read summarization, e.g.:

```
annot="/usr/refs/Homo_sapiens/NCBI/build37.2/Annotation/genes_ercc.gtf"
```

This will create a new directory named `/featureCounts_out`, which contains a tab-delimited text file with genes as rows, and columns as single-cell libraries.

This matrix is ready to be imported into SCell for quality control and analysis.

#### 1.4. Data Import

Refer to the [Load Data and Metadata](https://github.com/carmensandoval/SCell/blob/master/Manual.md#1-load-data-and-metadata) section of the [SCell Manual](https://github.com/carmensandoval/SCell/blob/master/Manual.md#using-scell) for detailed instructions on how to import your data and start using SCell.

## 2. Single-Cell Data Analysis using SCell

#### 2.1. Library Quality Control and Filtering

SCell can identify low quality libraries by computing their Lorenz statistic. Briefly, it estimates genes expressed at background levels in a given sample, and filters samples whose background fraction is significantly larger than average, via a threshold on the Benjamini-Hochberg corrected q-value.

In our tests, samples that have a small q-value for our Lorenz-statistic have low library complexity, as measured by Gini-Simpson index (Simpson, 1949), and/or they have low coverage, as estimated by the Good-Turing statistic (Good, 1953). Moreover, in our data the Lorenz-statistic correlates with the results of live-dead staining (Pearson-correlation 0.7), marking as outliers cells which appeared red, as well as empty chambers.

Libraries may also be manually filtered and excluded from further analysis based on other criteria shown in the SCell expression profiler, such as live/dead staining image call or number of genes tagged.

Refer to the [Library Quality Control / Outlier Filtering](https://github.com/carmensandoval/SCell/blob/master/Manual.md#2-library-quality-control--outlier-filtering) section of the SCell Manual for detailed instructions on how to perform library quality control and filtering.

#### 2.2 Feature Selection

It is important to select a discriminating panel of genes prior to dimensionality reduction. An ideal gene panel for dimensionality reduction is one that is sampled in a large number of cells, while at the same time exhibiting sufficient inter-cellular variance as to distinguish disparate cell types. In the low-coverage regime, typical of single-cell data, it can be difficult to distinguish when a gene is under-sampled due to technical limitations from when it is lowly expressed.
SCell provides statistics for feature selection. It uses a score statistic derived from a generalized-Poisson model, to test for zero-inflation in each gene’s expression across cells.  

Refer to the [Feature Selection](https://github.com/carmensandoval/SCell/blob/master/Manual.md#feature-selection) section of the SCell Manual for detailed instructions on how to select a meaningful gene set for analysis based on these metrics.

#### 2.3. Library Normalization

SCell can perform count normalization in several ways:

##### a) Normalization by library size (counts per million)

Refer to the [Normalization by Library Size](https://github.com/carmensandoval/SCell/blob/master/Manual.md#3a-feature-selection-and-normalization-by-library-size-cpm) section of the SCell manual.

##### b) Latent variable regression based on cyclin/CDK expression

SCell utilizes canonical-correlation analysis (CCA) on cyclins/CDKs to correlate cell-cycle and gene expression. It will estimate the percentage of genome-wide variance explained by cyclin/CDKs, the specific cyclins/CDKs that explain the highest percentage of variance and the genes that correlate most strongly with cyclin/CDKs.

Additionally, SCell can regress out unwanted variation, using iteratively re-weighted least squares regression. SCell can produce counts normalized by any combination of:

1. Cyclins and cyclin-dependent kinases (a model of cell cycle state)

2. A user supplied count matrix (enabling an arbitrary set of controls).

To normalize for the effect of cell cycle on gene expression, in addition to normalizing libraries by size, refer to the  SCell manual section [Normalization by Human Cyclins and Cyclin-Dependent Kinases (Cell Cycle Regression).](https://github.com/carmensandoval/SCell/blob/master/Manual.md#3b-normalization-by-human-cyclins-and-cyclin-dependent-kinases-cell-cycle-regression)

##### c) Mutual background

#### 2.4. Dimensionality Reduction by PCA

SCell implements PCA for dimensionality reduction, and optionally Varimax-rotation to post-process the PCA. Varimax rotates the PCA axes in order to reduce the number of genes strongly loading two axes. By default, PC1 and PC2 are the axes shown on the PCA score and loading plots, but this can be specified by the user.

To perform PCA on pre-filtered, normalized samples after selecting a meaningful gene panel, refer to the  SCell manual section [Dimensionality Reduction).](https://github.com/carmensandoval/SCell/blob/master/Manual.md#4-dimensionality-reduction-pca)

#### 2.5. Clustering

SCell implements several methods for clustering:

- k-means
- Minkowski-weighted k-means
- Gaussian mixture model
- DBSCAN

#### 2.6 Visualizing Gene Expression

SCell provides options for fitting loess/lowess regressions, as well as several interpolation algorithms (linear, cubic spline, biharmonic and thin-plate spline) to gene expression across samples in PCA space. This allows for the visualization of expression
gradients, and the evaluation of gene expression kinetics along minimum-spanning trees.


#### 2.7 Minimum-Spanning Trees and Lineage Trajectory

Viewing Gene Expression across trees

#### 2.8 Iterative PCA and Clustering

A PCA can be recomputed at any time from the user's sample list, enabling "iterative" PCA learning of population sub-structure.
