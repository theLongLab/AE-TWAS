# Users’ Manual of AE-TWAS
## Overview
Transcriptome-wide association study (TWAS) is an emerging model leveraging gene expressions to direct genotype-phenotype association mapping. A key component in TWAS is the prediction of gene expressions; and many statistical approaches have been developed along this line. However, a problem is that many genes have low expression heritability, limiting the performance of any predictive model. In this work, hypothesizing that appropriate denoising may improve the quality of expression data (including heritability), we propose AE-TWAS, which adds a transformation step before conducting standard TWAS. The transformation is composed of two steps by first splitting the whole transcriptome into co-expression networks (modules) and then using autoencoder (AE) to reconstruct the transcriptome data within each module. This transformation removes noise (including nonlinear ones) from the transcriptome data, paving the path for downstream TWAS. We applied AE-TWAS to the GTEx whole blood transcriptome data and GWAS data of five human diseases, showing two inspiring properties of AE-TWAS: (1) After transformation, the transcriptome data enjoy higher expression heritability at the low-heritability spectrum and possess higher connectivity within the modules. (2) The transferred transcriptome indeed enables better performance of TWAS; and moreover, the newly formed highly connected genes (i.e., hub genes) are more functionally relevant to diseases, evidenced by their functional annotations and overlap with TWAS hits. Taking together, we show that autoencoder transformation produces “better” transcriptome, which in turn enables improved expression-assisted genotype-phenotype association mapping. The impact of this work may be beyond the field of gene mapping: AE can be deemed as a nonlinear extension of principal component analysis (PCA) that is used for removing artifacts in expression data routinely. As such, this work may inspire more expression-based applications to be carried out after an appropriate AE-transformation, unlocking the use of AE-denoised transcriptome in many fields.

![My Image](Fig1.png)

## Installation
**Step1:** Preprocess transcriptome data by clustering genes into distinct modules using weighted gene co-expression network analysis (WGCNA).

We chose whole-blood gene expression data with 670 subjects as transcriptome and used gencode.v26.GRCh38.genes.gtf as gene model file downloaded from Genotype-Tissue Expression Project version 8 (GTEx v8) (https://gtexportal.org/home/datasets). Covariates, including genotyping principal components (PCs), were obtained from GTEx portal (https://gtexportal.org/home/datasets). For each gene, we adjusted the gene expression for the top five genotyping PCs, age, sex, sequencing platform, PCR protocol, and 60 confounding factors using a probabilistic estimation of expression residuals (PEER) analysis. There is a description of how to download and use the PEER tool: https://github.com/PMBio/peer/wiki/Tutorial. 

WGCNA is an R package consisting of a comprehensive collection of R functions for performing various aspects of weighted correlation network analysis. Users need to first install R and R studio, and then install the WGCNA package (https://cran.r-project.org/web/packages/WGCNA/WGCNA.pdf). We followed a step-by-step protocol for network construction and module detection. The power $\beta$ in adjacency matrix was set as 16, which was the lowest at which the scale-free topology fit index reached 0.8.

## Post analysis

## Contacts

## Copyright License (MIT Open Source)


