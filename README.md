# Genetic control of non-coding RNAs in the human brain and their implications for complex traits

## Authors
- Li Chen<sup>1</sup> (chenli@westlake.edu.cn)
- Yazhou Guo<sup>1</sup>
- Junren Hou<sup>1</sup>
- Wen Yang<sup>1</sup>
- Ting Qi<sup>1,2</sup>
- Jian Yang<sup>1,2,*</sup> (jian.yang@westlake.edu.cn)

<sup>1</sup>School of Life Sciences, Westlake University, Hangzhou, Zhejiang 310024, China  
<sup>2</sup>Westlake Laboratory of Life Sciences and Biomedicine, Hangzhou, Zhejiang 310024, China  
<sup>*</sup>Correspondence: Jian Yang (jian.yang@westlake.edu.cn)

## Abstract
Non-coding RNAs (ncRNAs) play crucial roles in the regulation of gene expression, but their genetic underpinnings and roles in human traits and diseases remain largely elusive. Here, we identified 38,441 long non-coding RNAs (lncRNAs) and 23,548 circular RNAs (circRNAs) from RNA-seq data of 2,865 human cortex samples, of which 18,429 lncRNAs and all circRNAs were not reported in GENCODE. Expression quantitative trait locus (eQTL) analyses identified cis-eQTLs for 15,362 lncRNAs and 1,312 circRNAs. We showed that lncRNA- or circRNA-eQTLs were largely independent of, and had larger effects on average than, eQTLs of their adjacent or parental protein-coding genes (PCGs). The circRNA eQTLs were highly enriched in canonical splice sites, highlighting the crucial role of back-splicing in circRNA biogenesis. LncRNA eQTLs were enriched for heritability of brain-related complex traits and associated with 72 (11.2%) of the colocalized GWAS signals that showed no evidence of colocalization with PCG eQTLs or splicing QTLs identified in the same dataset. We showcased lncRNAs (e.g., those near VPS45, MAPT, and RGS6) and circRNAs (e.g., that for GRIN2A) that may be implicated in complex traits through genetic regulation of ncRNAs. Our study provides new insights into the genetic regulation of ncRNAs and their implications in brain-related complex traits.

## Repository Structure and Code Descriptions

### 1. Identification and quantification of lncRNAs and circRNAs
Contains scripts for identifying and quantifying lncRNAs and circRNAs in 10 primary brain RNA-seq datasets.

**Files:**
- `lncRNA_circRNA_identification.zip`

### 2. eQTL mapping analysis of lncRNAs and circRNAs
Contains command files for performing eQTL mapping for lncRNAs and circRNAs.

**Files:**
- `eQTL_mapping/HBCC.circRNA.eQTL.com.txt`
- `eQTL_mapping/HBCC_lncRNA_eQTL.command.txt`

### 3. Heritability enrichment analysis and mediated heritability estimation
Contains scripts for performing heritability enrichment analysis and mediated heritability estimation for lncRNAs and circRNAs.

**Files:**
- `Heritability_enrichment/LDSC_MESC.txt`

### 4. Integrative analysis of GWAS and eQTL data
Contains command files for integrative analysis of GWAS and eQTL data.

**Files:**
- `COLOC_SMR/com.txt`

### 5. Locus plot generation
Contains R scripts for generating locus plots for specific genomic regions of interest.

**Files:**
- `BrainMetalnc.22630-SCZ.R`
- `ENSG00000285184-SCZ.R`
- `GRIN2A-SCZ-circ.R`
- `noPCG-BrainMetalnc.28095-CP.R`
- `noPCG-BrainMetalnc.43391-SCZ.R`

### 6. Figure generation
Contains R code for generating figures presented in the manuscript.

**Files:**
- `all_command_figs.txt`

## Data Availability

### Online Resources
- **Interactive Query Tool**: An online tool for querying lncRNA- and circRNA-eQTLs is available at:  
  https://yanglab.westlake.edu.cn/data/brainmeta

- **Summary Statistics**: Full lncRNA- and circRNA-eQTL summary statistics are available at:  
  https://yanglab.westlake.edu.cn/pub_data.html

## Acknowledgements
This research was partially funded by the National Natural Science Foundation of China (U23A20165 and 32100493), "Pioneer & Leading Goose" R&D Program (2024SSYS0032 and 2022SDXHDX0001), the Leading Innovative and Entrepreneur Team Introduction Program (2021R01013), and the Westlake University Research Center for Industries of the Future (WU2022C002 and WU2023C010). 

We thank the Westlake University High-Performance Computing Center for assistance in computing. This study makes use of data from HRS (dbGaP accession: phs000428), GTEx (dbGaP accession: phs000424), PsychENCODE (Synapse accession: syn4921369), and AMP-AD (Synapse accession: syn5550382).
