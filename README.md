codes for generating the results in the manuscript  "Genetic control of non-coding RNAs in the human brain and their implications for complex traits".

Li Chen1, Yazhou Guo1, Junren Hou1, Wen Yang1, Ting Qi1,2, Jian Yang1,2,*

1School of Life Sciences, Westlake University, Hangzhou, Zhejiang 310024, China
2Westlake Laboratory of Life Sciences and Biomedicine, Hangzhou, Zhejiang 310024, China
*Correspondence: Jian Yang (jian.yang@westlake.edu.cn)

Abstract
Non-coding RNAs (ncRNAs) play crucial roles in the regulation of gene expression, but their genetic underpinnings and roles in human traits and diseases remain largely elusive. Here, we identified 38,441 long non-coding RNAs (lncRNAs) and 23,548 circular RNAs (circRNAs) from RNA-seq data of 2,865 human cortex samples, of which 18,429 lncRNAs and all circRNAs were not reported in GENCODE. Expression quantitative trait locus (eQTL) analyses identified cis-eQTLs for 15,362 lncRNAs and 1,312 circRNAs. We showed that lncRNA- or circRNA-eQTLs were largely independent of, and had larger effects on average than, eQTLs of their adjacent or parental protein-coding genes (PCGs). The circRNA eQTLs were highly enriched in canonical splice sites, highlighting the crucial role of back-splicing in circRNA biogenesis. LncRNA eQTLs were enriched for heritability of brain-related complex traits and associated with 72 (11.2%) of the colocalized GWAS signals that showed no evidence of colocalization with PCG eQTLs or splicing QTLs identified in the same dataset. We showcased lncRNAs (e.g., those near VPS45, MAPT, and RGS6) and circRNAs (e.g., that for GRIN2A) that may be implicated in complex traits through genetic regulation of ncRNAs. Our study provides new insights into the genetic regulation of ncRNAs and their implications in brain-related complex traits.


### Identification and quantification of lncRNAs and circRNAs in the 10 primary brain RNA-seq datasets
lncRNA_circRNA_identification.zip

### eQTL mapping analysis of lncRNAs and circRNAs
eQTL_mapping/HBCC.circRNA.eQTL.com.txt
eQTL_mapping/HBCC_lncRNA_eQTL.command.txt

### Heritability enrichment analysis and mediated heritability estimation for lncRNAs and circRNAs
Heritability_enrichment/LDSC_MESC.txt

### Integrative analysis of GWAS and eQTL data 
COLOC_SMR/com.txt

### Locus plot
BrainMetalnc.22630-SCZ.R
ENSG00000285184-SCZ.R
GRIN2A-SCZ-circ.R
noPCG-BrainMetalnc.28095-CP.R
noPCG-BrainMetalnc.43391-SCZ.R

### R codes for generating figures
all_command_figs.txt

### lncRNA or circRNA Summary statistics data availability
Online tool for querying lncRNA- and circRNA-eQTLs: https://yanglab.westlake.edu.cn/data/brainmeta. 
Full lncRNA- and circRNA-eQTL summary statistics: https://yanglab.westlake.edu.cn/pub_data.html.


Acknowledgements
This research was partially funded by the the National Natural Science Foundation of China (U23A20165 and 32100493), "Pioneer & Leading Goose" R&D Program (2024SSYS0032 and 2022SDXHDX0001), the Leading Innovative and Entrepreneur Team Introduction Program (2021R01013), and the Westlake University Research Center for Industries of the Future (WU2022C002 and WU2023C010). We thank the Westlake University High-Performance Computing Center for assistance in computing. This study makes use of data from HRS (dbGaP accession: phs000428), GTEx (dbGaP accession: phs000424), PsychENCODE (Synapse accession: syn4921369), and AMP-AD (Synapse accession: syn5550382). 
