#SCell user guide

>**SCell** is an integrated software tool for quality filtering, normalization, feature selection, iterative dimensionality reduction, clustering and the estimation of gene-expression gradients from large ensembles of single-cell RNA-seq datasets. SCell is open source, and implemented with an intuitive graphical interface.
>
>Binary executables for Windows, MacOS and Linux are available at http://sourceforge.net/projects/scell, source code and pre-processing scripts are available from https://github.com/diazlab/SCell. If you use this software for your research, please cite:

## Table of contents
- [Installation](#installation)
- [Load data and metadata](#load_data)
- [Quality control & outlier filtering](#qc)
- [Gene panel selection and normalization](#norm)
- [Dimensionality reduction and PCA](#pca)
- [Visualize gene expression across cells](#visualization)
- [Clustering](#clustering)
- [Fit a minimum-spanning tree or minimum-cost path](#fit_mst)
- [Iterative PCA](#iter_pca)
- [Save and load sessions](#save)

## 1. <a id="installation"></a>Installation
SCell runs under 64bit Mac OSX, Windows, and Linux. Start by downloading the appropriate installation package
from:
[https://sourceforge.net/projects/scell/](https://sourceforge.net/projects/scell/) 

SCell is released under the GNU General Public License:
[http://www.gnu.org/licenses/](http://www.gnu.org/licenses/). The SCell source code can
be obtained from
[https://github.com/diazlab/scell](https://github.com/songlab/scell).

## 2. <a id="load_data"></a>Load data and metadata

SCell accepts as input a matrix of raw gene counts with **genes as rows and cells as columns.** The table should have the format shown below.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/CountsSample.png" width=550>

To import your gene expression dataset into SCell, select **Load data** in the main window, then choose your raw gene counts file from the dialog box.

![MainWindow](Images/main_window1.png?raw=true)


SCell also accepts a metadata table associated with the libraries in your dataset. The metadata matrix can include:

- a live-dead image call from your single-cell workflow
- number of mapped reads
- number of unmapped reads

The metadata file must have the format shown below.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/MetadataSample.png" width="300">

To import your metadata into SCell, select **Load meta-data** in the main window, then choose the corresponding file from the dialog box.

The library names and associated metadata for the imported cells is displayed in the interactive expression profiler in the main window.

## 3. <a id="qc"></a>Library Quality Control / Outlier Filtering

SCell can identify low quality libraries by computing their Lorenz statistic. Briefly, it estimates genes expressed at background levels in a given sample, and filters samples whose background fraction is significantly larger than average, via a threshold on the (Benjamini-Hochberg corrected) q-value.

In our tests, samples that have a small q-value for our Lorenz-statistic have low library complexity, as measured by Gini-Simpson index (Simpson, 1949), and they have low coverage, as estimated by the Good-Turing statistic (Good, 1953). Moreover, in our data the Lorenz-statistic correlates with the results of live-dead staining (Pearson-correlation 0.7).

To perform QC on the loaded and selected single-cell RNA-seq libraries, select **QC selected libraries** in the lower panel of the main window.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/QCButton.png" width="350">

SCell will estimate library complexity, Lorenz statistic and PRESEQ score for the selected libraries, and will generate the plots shown below.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/libraryComplexityProgress.png" width="300">

___Lorenz Curves___

[Add figure legend.]

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/Lorenz.png" width="300">

___Per-Cell Gene Expression Distribution___

[Add figure legend.]

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/GeneExpressionQuantiles.png" width="300">

___Simpson Diversity and PRESEQ Coverage Boxplots___

[Add figure legend.]

SCell reports library coverage and marginal return estimates based on the PRESEQ (Daley and Smith, 2014) extrapolation curve.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/QCBoxplots.png" width="300">

SCell displays these quality metrics, Gini-Simpson index, along with the user’s metadata, in the interactive expression profiler.

Once these metrics have been estimated, **filter low-complexity or outlier cells** by selecting **Lorenz** or **PRESEQ** from the dropdown menu in the main window lower panel, and setting a threshold to filter by (default=0.05, Lorenz).

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/FIlterByButton.png" width="300">

Then, click **Filter by** to unckeck cells that did not meet the established threshold in the interactive expression profiler. This will cause low quality/outlier libraries to be excluded from all downstream analysis.

You may also de-select libraries that you wish to exclude based on other criteria shown in the expression profiler, such as live/dead staining image call or number of genes tagged.

## 4. <a id="norm"></a>Gene Panel Selection and Library Normalization

Once your single-cell data has been filtered for low-quality samples, SCell can perform normalization of your libraries in several ways:

a) Normalization by library size (counts per million)

b) Latent variable regression based on cyclin/CDK expression

c) Mutual background latent variable regression

### Feature Selection
It is also important to select a discriminating panel of genes useful for dimensionality reduction.

An ideal gene for dimensionality reduction is one that is sampled in a large number of cells, while at the same time exhibiting sufficient inter-cellular variance as to distinguish disparate cell types. In the low-coverage regime, typical of single-cell data, it can be difficult to distinguish when a gene is under-sampled due to technical limitations from when it is lowly expressed.

SCell provides statistics for feature selection. It uses a score statistic derived from a generalized-Poisson model, to test for zero-inflation in each gene’s expression across cells.  

#### 4a. Feature Selection and Normalization by Library size (CPM)

To normalize samples by library size and select meaningful genes for analysis, follow these steps:

- In the main panel, select the **Normalize selected libraries** button to launch the **normalization tool** window.

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/NormalizeButton.png" width="300">

- Once in the normalization tool, click on **Select Genes**.

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/normTool.png" width="250">

- SCell will perform gene variance and zero-inflation analysis, and a **Gene Selection Tool** window will be launched.

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/geneVarianceProgress.png" width="250">

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/ZIprogress.png" width="250">

- In the **Gene Selection Tool** window, set an Index of Dispersion percentile threshold for genes, as well as a threshold for the fraction of cells expressing a given gene. You can also set a zero-inflation power threshold. On the lower panel, SCell displays a list of the genes in your dataset and their values for these metrics. This list, as well as the displayed plot, can be exported by selecting the **Export plot** or **Export gene list** buttons.

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/GeneSelection.png" width="350">

- Once the thresholds for a gene panel have been chosen, click **Use these genes**.
Only the genes that meet your criteria will be used for downstream analysis.

- You will be taken back to the **Normalization Tool** window. Click **Done**.

    <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/normToolDone.png" width="250">

Your libraries are now normalized by library size (counts per million) .You will be taken back to the Main Panel, where you may proceed to dimensionality reduction and clustering of the filtered, normalized libraries.

#### 4b. Normalization by Human Cyclins and Cyclin-Dependent Kinases (Cell Cycle Regression)

To normalize for the effect of cell cycle on gene expression (remove unwanted variation due to cell cycle state), in addition to normalizing libraries by size, follow these steps:

- Follow steps 1 to 5 from the **Normalization by Library size (CPM)** section above, in order to select a meaningful gene panel to include in the analysis.

- On the **Normalization Tool** window, check the **Human cyclins/CDKs** box, then click **Normalize samples**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/cyclinNormalization.png" width="250">

SCell will estimate the percentage of genome-wide variance explained by cyclin/CDKs, the specific cyclins/CDKs that explain the highest percentage of variance and the genes that correlate most strongly with cyclin/CDKs.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/cyclinNormProgress.png" width="250">

A barchart will be displayed, showing Cyclin/CDK contribution to the observed variance in gene expression.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/cyclinCorrelations.png" width="400">

You will be prompted to specify the location to store the file for cyclin/CDK correlations, then for cyclin/CDK-gene correlations.

Choose a location in the dialog box, then click **Save**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/saveCyclinCorr.png" width="400">
<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/saveGeneCorrToCyclins.png" width="400">

SCell utilizes canonical-correlation analysis (CCA) on cyclins/CDKs to correlate cell-cycle and gene expression.

___Cyclin-CDK Correlations File___

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/cyclinCorrTSV.png" width="250">

___Cyclin-CDK/Gene Correlations File___

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/geneCorrWithCyclinsTSV.png" width="350">

After specifying a location to save the gene-cyclin correlation and cyclin-cyclin correlation TSV files, normalization will take place. Please allow a few minutes for its completion.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/normalizationCyclins_Progress.png" width="300">

Once the libraries have been normalized, you will be taken back to the Main Panel, where you may proceed to dimensionality reduction and clustering of the filtered, normalized libraries.

Note: SCell can produce counts normalized by any combination of:

- Cyclins and cyclin-dependent kinases (a model of cell cycle state)

- A user supplied count matrix (enabling an arbitrary set of controls).

## 5. <a id="pca"></a>Dimensionality Reduction (PCA)

SCell implements PCA for dimensionality reduction, and optionally Varimax-rotation.

On the lower panel of the Main Window, select **Analyze selected libraries** to perform dimensionality reduction via PCA on the QC-filtered and normalized libraries.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/AnalyzeButton.png" width="350">

Two interactive plots will be displayed, which allow the user to explore samples in PCA space, with gene-level and sample-level metadata displayed upon mouse-over:


  ○ ___Interactive PCA Sample Scores Plot___

View samples in PCA space.

Mouse-over samples to view their identity and associated metrics, or click on a sample to mark and add it to the **working sample list**.
By default, PC1 and PC2 are the axes shown on the plot, but this can be modified (see below).

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PCAScores.png" width="300"> <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/CellAnnotations.png" width="300">

  ○ ___Interactive Gene Loadings Plot___

Explore strongly loading genes in your dataset.

In this dataset, **many canonical neuronal** genes are loading negative PC1.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PCAGenesNeuronal.png" width="300">
<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PCAGenesNeuronal2.png" width="300">

Many **canonical radial glia (neocortical neural stem cells)** marker genes strongly load positive PC1.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PCAGenesRG.png" width="300">
<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PCAGenesRG2.png" width="300">

Mouse over genes in the **Gene Loading plot** to view their gene symbol, mean expression, IOD and the percentage of cells expressing the gene in the **Gene Annotations** panel.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/GeneAnnotations.png" width="300">

Click on genes on the **gene loading plot** to add them to the **Working gene list**. Alternatively, to find a gene in the loading plot, type its gene symbol in the **Working gene list** panel search box, then click **Add/Find gene**. The gene will be highlighted on the loading plot, and added to the working list.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/workingGeneList_Search.png" width="300">

The current working gene list can be exported to a file for downstream analysis, by selecting **Save gene list to file** in the **Working gene list** panel.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/workingGeneListSave.png" width="300">

You will be prompted to choose a location to store the file.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/saveLoadingGeneList.png" width="300">

#### Refresh PCA Axes

 By default, PC1 and PC2 are the axes shown on the PCA score and loading plots, but this can be modified. Enter the axes you would like to plot, then select **Refresh PCA axes** to generate new Sample Scores and Gene Loadings plots.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/refreshPCAaxes.png" width=350>

#### Varimax Rotation

You can choose to apply Varimax Rotation to post-process the PCA. This rotates the PCA axes in order to reduce the number of genes strongly loading two axes.

To apply Varimax rotation, simply click **Varimax rotate**. The plots will be updated accordingly.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/varimaxRotate.png" width=350>


## 6. <a id="visualization"></a>Visualize Gene Expression across Cells.

Enter the name of the gene whose expression you want to visualize.

Select **surface** or **contour**, and a regression method from the drop-down menu, then click on **Plot expression** to generate a plot of gene expression across cells in PCA space.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/VisGeneExpressionButton.png" width="300">

In this dataset, we will plot the expression of _VIM_ and _DCX_, a radial glia and newborn neuron marker, respectively.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/VIMExpressionContourPlot.png" width="300"> <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/VIMExpressionSurface.png" width="297">

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/DCXExpressionContourPlot.png" width="300">
<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/DCXExpressionSurfacePlot.png" width="297">

The height of the points in the surface plot corresponds to their normalized gene expression value for the selected gene.

## 7. <a id="clustering"></a>Clustering

SCell can implement 4 different clustering  algorithms:

○ k-means

○ Gaussian mixture

○ Minkowski weighted k-means

○ DBSCAN

In the **Clustering panel**, select one of the clustering algorithms from the drop-down menu, then click  **Cluster cells**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/ClusterMenu.png" width="300">

A new window will prompt you to enter the desired number of clusters, as well as the number of replicates and distance metric to use. Enter these parameters, then click **Done**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/kmeansParameters.png" width="302">

###### Example: K-means Cluster Assignments

A plot will be generated with samples color-coded according to their assigned cluster.

Upon mouse-over, sample information is displayed in the **Cell Annotation** panel.

Upon selection, samples are added to the **Working sample list**, which can be exported to a file, or used to run a new PCA iteration (see below).

<img src= "https://dl.dropboxusercontent.com/u/9990581/Scell/SCell_Screenshots/kmeansClusters.png" width="600">


## 8. <a id="mst_fit"></a>Fitting a Minimum Spanning Tree or minimum cost path

Once cells have been clustered, SCell can compute a Minimum Spanning Tree to predict a lineage trajectory across clusters.

Select a **root cluster** from the drop-down menu in the **Fit MST** panel, then click **Fit Tree**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/MSTButton.png" width="300">

In this working example, we have selected **Cluster 2** as the root cluster, given that its cells express high levels of canonical radial glia marker genes.

The MST trajectory passes through Intermediate Progenitor Cells (IPCs) **(Cluster 4)** and continues on to neurons **(Cluster 1)**.

Note that a branching event also occurs at **Cluster 4**. The cells in **Cluster 3** have the transcriptional profile corresponding to interneurons, a neuronal population which stems from a distinct lineage and region in the brain, and which migrates to the neocortex during development.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/MST_withVarimax.png" width=300>

#### Visualize Average Gene Expression Across MST Clusters

To view a profile of the transcriptional dynamics of a gene across the estimated lineage trajectory, repeat the steps outlined above in **5. Visualize Gene Expression across Cells**, after fitting a Minimum-Spanning Tree.

A surface / contour plot will be generated showing the gene's expression levels across cells in PCA space.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PAX6Surface.png" width=350>

Additionally, the following plots will be displayed, describing the expression level dynamics of the gene across the defined clusters, following the estimated trajectory.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PAX6expressionClusters.png" width="300"> <img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/PAX6ExpressionClusters2.png" width="295">

## 9. <a id="iter_pca"></a>Iterative PCA

SCell can perform PCA on a subset of libraries selected from the first round of analysis.

Select the cells from the cluster you wish to further analyze on the displayed sample scatterplot.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/iterativePCAselect.png" width="400">

The selected samples will appear on the **Working Sample List** in the main panel.

Click **Refresh PCA using samples** to perform a new iteration of PCA on the selected libraries only.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/iterativePCArefresh.png" width="400">

A new **sample scores** plot will be displayed containing only the selected libraries in PCA space. You can choose to apply Varimax rotation as described above, or to change the principal components plotted on each axis.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/iterativePCAsampleScores.png" width="600">

A new **gene loading** plot will be displayed as well, showing the new contribution of genes to each principal component shown.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/iterativePCAplot.png" width="400"><img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/iterativePCAsampleList.png" width="400">

The expression of selected genes across this subset of samples can be displayed as described above.

Samples can be clustered and a new minimum-spanning tree can be fit to further characterize the subset of cells being analyzed.

## 10. <a id="save"></a>Saving and Loading Sessions

#### Save Current Session

To save the current session, including library QC and filtering, feature selection and normalization, select **Save Session** on the main window.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/saveSessionButton.png" width="350">

You will be prompted to choose a location to store the **.mat** session. Choose a name for the session, and click **Save**.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/saveSession.png" width="350">

#### Load Session from File

To restore a previously saved session into SCell, select **Load data** on the main window.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/mainPanelLoadSession.png" width="350">

You will be prompted to select the .mat file you wish to load. Once the file has been selected, select **Open**. Your data, selected libraries, feature selection and normalization, will be restored.

<img src= "https://dl.dropboxusercontent.com/u/9990581/SCell/SCell_Screenshots/loadSessionLocation.png" width="350">
