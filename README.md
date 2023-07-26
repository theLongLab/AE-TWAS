# Users’ Manual of AE-TWAS
## Overview
Transcriptome-wide association study (TWAS) is an emerging model leveraging gene expressions to direct genotype-phenotype association mapping. A key component in TWAS is the prediction of gene expressions; and many statistical approaches have been developed along this line. However, a problem is that many genes have low expression heritability, limiting the performance of any predictive model. In this work, hypothesizing that appropriate denoising may improve the quality of expression data (including heritability), we propose AE-TWAS, which adds a transformation step before conducting standard TWAS. The transformation is composed of two steps by first splitting the whole transcriptome into co-expression networks (modules) and then using autoencoder (AE) to reconstruct the transcriptome data within each module. This transformation removes noise (including nonlinear ones) from the transcriptome data, paving the path for downstream TWAS. We applied AE-TWAS to the GTEx whole blood transcriptome data and GWAS data of five human diseases, showing two inspiring properties of AE-TWAS: (1) After transformation, the transcriptome data enjoy higher expression heritability at the low-heritability spectrum and possess higher connectivity within the modules. (2) The transferred transcriptome indeed enables better performance of TWAS; and moreover, the newly formed highly connected genes (i.e., hub genes) are more functionally relevant to diseases, evidenced by their functional annotations and overlap with TWAS hits. Taking together, we show that autoencoder transformation produces “better” transcriptome, which in turn enables improved expression-assisted genotype-phenotype association mapping. The impact of this work may be beyond the field of gene mapping: AE can be deemed as a nonlinear extension of principal component analysis (PCA) that is used for removing artifacts in expression data routinely. As such, this work may inspire more expression-based applications to be carried out after an appropriate AE-transformation, unlocking the use of AE-denoised transcriptome in many fields.   

## Installation

## Post analysis

## Contacts

## Copyright License (MIT Open Source)


