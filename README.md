# scRNAseq-preprocessing
Documentation of and code for analyses performed to assess the effect of different methods for cleanup and pre-processing of raw scRNAseq data on downtream analyses.

## Data
In order to assess different tools and options for pre-processing scRNA-seq data, we downloaded publicly available data sets that were generated with 10x Chromium chemistry and sequenced on various contemporary Illumina sequencing instruments. We focused on datasets generated for mouse (*Mus musculus*) as the genome assembly and annotation are of exceptionally high quality, and there are no issues regarding patient anonymity that restrict data access as in the case human data. Below are the datasets we analyzed.

| Species | Strain | Sample ID | Tissue | Sequencing Chemistry | Platform | Estimated Cells | Mean Reads Per Cell | Median Genes Per Cell | Data Source | Fastq Link(s) | Notes |
|---------|--------|-----------|--------|----------------------|----------|-----------------|---------------------|-----------------------|-------------|---------------|-------|
|  Mouse  |   NA   | neuron_10k_v3 |cortex, hippocampus and sub ventricular zone | Universal 3' Gene Expression v3 | NovaSeq | 11,831 | 30,184 | 3,684 | [10x](https://www.10xgenomics.com/datasets/10-k-brain-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0) | same as data source | E18 developmental stage |
