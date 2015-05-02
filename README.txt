1. Implement a 3 window GUI, 2 scatterplots with PCA scores and feature loadings respectively, and one readout window. The readout window is partially implemented as pca_tool.fig
a. Information about the samples is displayed in textboxes as the user mouses over a sample in the plot. 
b. Buttons or menu items on the readout window allow the user to cluster the points by k-means, Gaussian mixture model, Minkowski weighted k-means, or a user defined list. The first two algorithms are implemented in MATLAB and the Minkowski weighted k-means I have an implementation of. So you wouldn't need to implement those algorithms, just collect the parameters from a gui interface (either as a pop-up or in a preferences menu item), invoke the clustering routine and refresh the PCA plots, color coded so that samples within the same cluster have the same color.
c. The user should be able to refresh the plot with different PCA axes (2 & 3 instead of 1 & 2 for example)
d. The user should be able to select samples from the PCA score plot, and regenerate a new PCA with just those samples.
e. The PCA itself is generated from the data matrix d.nrmC. Rows are features, columns are samples. See below.
2. The code in part 1 should interface with code already written that does pre-processing.
There are additional modules that I will work on myself, and after the goals I describe above are met, I may have additional work involving testing, compiling the code for major platforms and compiling the underlying algorithms into command-line tools and C++ libraries.
3. The user should be able to curate feature (genes) and sample (cells) lists by clicking on the PCA plot windows. Additionally, the user should be able to add all features whose loading on a given PCA is above a threshold. This is sketched out in the bottom half of pca_tool.fig.The "Compute Gene ontology" button doesn't need to be implemented right now.   

I've implemneted alot of this functionality already in a previous academic code that can be resused for this project, it can be found at https://github.com/diazlab/HiTSelect. The application itself is very different, and you may have better ways to implement these things. 

pca_tool.m will be invoked from main.m It will take as input a data object d, with the following fields:

d = 

       slbls: {1x107 cell} - cell array of strings, sample labels
      counts: [14625x107 double] - (# features x # samples) matrix of feature counts
       nrmC : [14625X107 double] - (# features x # samples) matrix of normalized feature counts**These are the values to use for PCA**
         cpm: [14625x107 double] - (# features x # samples) matrix normalized counts (by sequenicng depth) "counts per million"
         ent: [1x14625 double] - positive integer feature meta-data, Entrez ID 
       gsymb: {14625x1 cell} - feature label, gene name
      length: [14625x1 double] - positive integer feature meta-data, gene length
        cidx: [1x107 double] - sample meta data, not relevant here
          qc: 1 - flag, not relevant here
          sf: [1x107 double] - sample meta-data, not relevant here
        gidx: [1x14625 double] - feature meta-data, not relevant here
         iod: [14625x1 double] - "Index of dispersion", feature meta-data displayed in data window
         pnz: [14625x1 double] - "Percent non-zeros", feature meta-data displayed in data window
     iod_fdr: [14625x1 double] - "Index of dispersion false discovery rate", feature meta-data not relevant here
    zinf_fdr: [14625x1 double] - feature meta-data, not relevant here
       niche: {1x107 cell} - sample meta-data, not relevant here
      source: {1x107 cell} - sample meta-data, not relevant here
       brain: [1x107 double] - sample meta-data, not relevant here
