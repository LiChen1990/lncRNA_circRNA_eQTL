library(dplyr)
library(reshape2)
library("data.table")
library("stringr")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("AnnotationHub")
library("org.Hs.eg.db")
library("plotgardener")
library("grid")


  

perl  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/get_input_for_locusplot.pl  /storage/yangjianLab/sharedata/GWAS_summary/01_Public/01_cojo/Cognitive_Performance_2018_NatGenet.txt  /storage/yangjianLab/guoyazhou/data/SNPs_database/dbSNP/dbSNP_b151/split_into_chrAuto/  /storage/yangjianLab/chenli/seq/PsychENCODE/CMC/merge/eQTL1/newlncRNA_exp.GeneLncRNA.hg19.bed  /storage/yangjianLab/guoyazhou/data/genecode_data/genecode_v43/v43lift37/gencode.v43lift37.annotation.gene.bed  BrainMetalnc.28095   1000000  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/BrainMeta.lncRNAeQTL.chr  > locus_input_data.BrainMetalnc.28095.CP.txt 


perl  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/get_input_for_locusplot.pl  /storage/yangjianLab/sharedata/GWAS_summary/01_Public/01_cojo/EA_NG_2018_excluding_23andMe.txt /storage/yangjianLab/guoyazhou/data/SNPs_database/dbSNP/dbSNP_b151/split_into_chrAuto/  /storage/yangjianLab/chenli/seq/PsychENCODE/CMC/merge/eQTL1/newlncRNA_exp.GeneLncRNA.hg19.bed  /storage/yangjianLab/guoyazhou/data/genecode_data/genecode_v43/v43lift37/gencode.v43lift37.annotation.gene.bed  BrainMetalnc.28095 1000000  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/BrainMeta.lncRNAeQTL.chr  > locus_input_data.BrainMetalnc.28095.EA.txt 


perl  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/get_input_for_locusplot.pl  /storage/yangjianLab/sharedata/GWAS_summary/01_Public/01_cojo/IQ_NG_2018.txt  /storage/yangjianLab/guoyazhou/data/SNPs_database/dbSNP/dbSNP_b151/split_into_chrAuto/  /storage/yangjianLab/chenli/seq/PsychENCODE/CMC/merge/eQTL1/newlncRNA_exp.GeneLncRNA.hg19.bed  /storage/yangjianLab/guoyazhou/data/genecode_data/genecode_v43/v43lift37/gencode.v43lift37.annotation.gene.bed  BrainMetalnc.28095  1000000  /storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/BrainMeta.lncRNAeQTL.chr  > locus_input_data.BrainMetalnc.28095.IQ.txt



 
####SCZ  
##Multiple genes in a single GWAS risk locus synergistically mediate aberrant synaptic development and function in human neurons
##lead snp 1886 chr17 44189067 1.185e-09 rs7225002 #2b5799
GWAS_trait_name="SCZ"
GWAS_file="/storage/yangjianLab/sharedata/GWAS_summary/01_Public/01_cojo/PGC3_SCZ_wave3_public_INFO80.txt"
probe=c("BrainMetalnc.28095")
probe_gene=c("ENSG00000136631")
highlight_probe=c("BrainMetalnc.28095","BrainMetalnc.31547")
window_size=500000
probe_list=probe
probe_gene_list=probe_gene
filename<-"/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/locus_input_data.ENSG00000264589.SCZ.txt"
#finemap_file<-"/storage/yangjianLab/chenli/seq/PsychENCODE/lncRNA_eQTL/enrichment/finemap/sczfinemap/scz.rs7225002.locus.CARMA.txt"


###For MAPT-AS1, some of the GWAS plots (e.g., neuroticism, EA and AD) do not seem to be very consistent with the lncRNA-eQTL signal. The GWAS signals for PD and SCZ do not seem to be genome-wide significant. If you really want to keep this example, I would suggest you move the CP GWAS result to the top of panel b, and remove GWAS plots for neuroticism, EA and AD.
## keep CP, PD, SCZ

dat<-read.table(file="locus_input_data.BrainMetalnc.28095.CP.txt",sep="\t",header=T);
gwasLocusData = dat[dat$type=="GWAS",]
start<-min(gwasLocusData$pos)-1000
end<-max(gwasLocusData$pos)+1000
chr<-paste0("chr",unique(gwasLocusData$chrom))

gwasLocusData$chrom=paste0("chr",gwasLocusData$chrom)
gwasLocusData<-gwasLocusData[,c("chrom","pos","P","SNP")]
gwasLocusData$color<-"#2b5799"
colnames(gwasLocusData)<-c("chrom","pos","p","SNP","color")
gwas_locus1=gwasLocusData


dat<-read.table(file="locus_input_data.BrainMetalnc.28095.EA.txt",sep="\t",header=T);
gwasLocusData = dat[dat$type=="GWAS",]
start<-min(gwasLocusData$pos)-1000
end<-max(gwasLocusData$pos)+1000
chr<-paste0("chr",unique(gwasLocusData$chrom))

gwasLocusData$chrom=paste0("chr",gwasLocusData$chrom)
gwasLocusData<-gwasLocusData[,c("chrom","pos","P","SNP")]
gwasLocusData$color<-"#2b5799"
colnames(gwasLocusData)<-c("chrom","pos","p","SNP","color")
gwas_locus2=gwasLocusData



filename<-"/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/locus_input_data.BrainMetalnc.28095.IQ.txt"
dat<-read.table(file=filename,sep="\t",header=T);
gwasLocusData = dat[dat$type=="GWAS",]
start<-min(gwasLocusData$pos)-1000
end<-max(gwasLocusData$pos)+1000
chr<-paste0("chr",unique(gwasLocusData$chrom))

gwasLocusData$chrom=paste0("chr",gwasLocusData$chrom)
gwasLocusData<-gwasLocusData[,c("chrom","pos","P","SNP")]
gwasLocusData$color<-"#2b5799"
colnames(gwasLocusData)<-c("chrom","pos","p","SNP","color")
gwas_locus3=gwasLocusData







lead_snp="rs2027349"
gwas_locus1[gwas_locus1$SNP==lead_snp,]
head(gwas_locus1)

leadSNP_p <- min(gwas_locus1$p)
leadSNP <- gwas_locus1[which(gwas_locus1$p == leadSNP_p), ]$SNP
gwas_locus1[gwas_locus1$SNP==leadSNP,]
#writeLines(gwas_locus1$SNP, paste0(leadSNP,".snplist.txt") )
chr
start
end
leadSNP
lead_snp<-leadSNP

#pip<-read.table(file=finemap_file,sep="\t",header=T);

str <- paste0("plink   --bfile /storage/yangjianLab/sharedata/LD_reference/UKB/genotype_10K/BED_ukbEUR_imp_v3_INFO0.8_maf0.01_mind0.05_geno0.05_hwe1e6_10K_chr",chr," --chr ",chr," --from-bp ", start, " --to-bp ", end, " --r2 --ld-window 10000000  --ld-window-kb 1000  --ld-window-r2 0 --ld-snp  ", leadSNP, " --out ",leadSNP) 
str


plink   --bfile /storage/yangjianLab/sharedata/LD_reference/UKB/genotype_10K/BED_ukbEUR_imp_v3_INFO0.8_maf0.01_mind0.05_geno0.05_hwe1e6_10K_chr16 --chr 16 --from-bp 52512821 --to-bp 54516395 --r2 --ld-window 10000000  --ld-window-kb 1000  --ld-window-r2 0 --ld-snp  rs8054299 --out rs8054299

plink   --bfile /storage/yangjianLab/sharedata/LD_reference/HRS/hrs_1kg_hwe1e-6_chr17  \
    --chr 17 \
    --from-bp 43375761  \
    --to-bp 44369215 \
    --r2 \
     --ld-window 10000000 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --ld-snp rs7225002 \
    --out rs7225002 
    
        

plink   --bfile /storage/yangjianLab/sharedata/LD_reference/UKB/genotype_10K/BED_ukbEUR_imp_v3_INFO0.8_maf0.01_mind0.05_geno0.05_hwe1e6_10K_chr17  \
    --chr 17 \
    --from-bp 43375761  \
    --to-bp 44369215 \
    --r2 \
    --ld-window 10000000 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --ld-snp rs7225002 \
    --out rs7225002
    
       
#ld<-read.table(file="rs12138231.ld",header=T);  
ld<-read.table(file=paste0(leadSNP,".ld"),header=T);  
colnames(ld)<-c("CHR_A","BP_A","SNP_A" ,"CHR_B", "BP_B", "SNP","R2")
dat_cov1<-merge(gwas_locus1,ld,by="SNP",all.x = TRUE)
gwas_locus1<-dat_cov1[,c("chrom","pos", "p" , "SNP","color","R2")]


index=which(gwas_locus1$R2 >=0.2)
gwas_locus1$color[index] = "#8BE9FF"
index=which(gwas_locus1$R2 >=0.4)
gwas_locus1$color[index] = "green"
index=which(gwas_locus1$R2 >=0.6)
gwas_locus1$color[index] = "orange"
index=which(gwas_locus1$R2 >=0.8)
gwas_locus1$color[index] = "red"


gwasLocusPlot=function(){
#  plot -------------------------------------------------------------------------
  pageCreate(width = 20, height = 30, default.units = "cm")
  
  #----plot position ----------------------------------------------------------------
  y0=1.5
  x0=2
  width0=14
  height0=2
  gap0=1
  #----gwas manhattan plot1 ------------------------------------------------------
  # y axis
  yMAX=round(max(-log10(gwas_locus1$p),na.rm = T))*1.25
  devbuf1 = yMAX/4
  y_axis_label=round(seq(0,yMAX,devbuf1),0)
  x_axis_label=round(seq(as.numeric(start),as.numeric(end), (as.numeric(end)-as.numeric(start))/4),0)
  
  ##https://phanstiellab.github.io/plotgardener/reference/plotManhattan.html
  manhattanPlot1 <- plotManhattan(
    data = gwas_locus1, 
    chrom = chr,
    chromstart = as.numeric(start),
    chromend = as.numeric(end),
    assembly = "hg19",
    fill = gwas_locus1$color,
    trans = "-log10",
    sigLine = TRUE, col = "darkgrey",
    baseline.color="black",
    cex=0.25,
    lty = 2, 
    range = c(0, yMAX),
    leadSNP = list(
     snp = leadSNP,
      pch = 18,
      cex = 0.5,
      fill = "red",
      fontsize = 18
    ),
    x = x0, y = y0, width = width0,
    height = height0,
    just = c("left", "top"),
    default.units = "cm"
  )
  
  ## Plot legend for R2
plotLegend(
    legend = c(
        "LD",
        paste("1.0", ">", "r^2",  "", ">=", "0.8"),
        paste("0.8", ">", "r^2",  "", ">=", "0.6"),
        paste("0.6", ">", "r^2",  "", ">=", "0.4"),
        paste("0.4", ">", "r^2",  "", ">=", "0.2"),
        paste("0.2", ">", "r^2",  "", ">=", "0")
    ),
    fill = c("grey","red","orange","green","#8BE9FF", "#2b5799" ), cex = 0.50,
    pch = c(19, 19,19,19,19,19), border = FALSE, x = x0+width0-3, y = y0-1.3,
    width = 1.5, height = 1.5, just = c("right", "top"),
    default.units = "cm"
)
  
  plotText(
    label = "CP", x = x0+0.2, y = y0, rot = 0,
    fontsize = 10, just = "left",
    default.units = "cm"
  )
  
  
  # annotation segments 
  annoSegments(
    x0 = unit(0, "npc"), y0 = 10,
    x1 = unit(1, "npc"), y1 = 10,
    plot = manhattanPlot1, default.units = "native",
    linecolor = "pink", lty = 2
  )
  
  ## Annotate genome label
  annoGenomeLabel(
    plot = manhattanPlot1, x = x0, y = y0+height0,
    fontsize = 8, scale = "Mb",
    just = c("left", "top"), default.units = "cm"
  )
  #----Annotate y-axis
  annoYaxis(
    plot = manhattanPlot1,
    at = y_axis_label,
    axisLine = TRUE, fontsize = 8
  )
  
#          ###annoHighlight
#  #https://phanstiellab.github.io/plotgardener/reference/annoHighlight.html
#  start_h<-gwas_locus1[gwas_locus1$SNP==lead_snp,]$pos
#  annoHighlight(
#    plot = manhattanPlot1,
#    chrom = chr,
#    chromstart = start_h-1000, chromend = start_h+1000,
#    y =  y0 , height = 25.5, just = c("left", "top"),
#    default.units = "cm",fill="red"
#) 
#
#plotText(
#    label = lead_snp, x = x0+6, y = y0 , rot = 0,
#    fontsize = 10, just = "left",
#    default.units = "cm"
#  )
  


###gwas_locus2
      probe_index=1
      gwas_locus2$pch <- 21
      gwas_locus2$CHR = chr
      gwas_locus2$color = "#2b5799"
      gwas_locus2 <- gwas_locus2[,c("CHR","pos","p","SNP","color")]
      colnames(gwas_locus2) <- c("chrom","pos","p","snp","color")
      
      #----plot lncRNA eQTL plot ----------------------------------------------------------------
      yMAX=round(max(-log10(gwas_locus2$p),na.rm = T))*1.25
      devbuf1 = yMAX/4
      y_axis_label=round(seq(0,yMAX,devbuf1),0)
      
      ##https://phanstiellab.github.io/plotgardener/reference/plotManhattan.html
      manhattanPlot2 <- plotManhattan(
        data = gwas_locus2, 
        chrom = chr,
        chromstart = start,
        chromend = end,
        assembly = "hg19",
        fill = gwas_locus2$color,
        trans = "-log10",
        sigLine = TRUE, col = "darkgrey",
        baseline.color="black",
        cex=0.25,
        lty = 2, 
        range = c(0, yMAX),
        # leadSNP = list(
        #   snp = lead_snp,
        #   pch = 18,
        #   cex = 0.6,
        #   fill = "#9632b8",
        #   fontsize = 8
        # ),
        x = x0, y =y0 + (height0 + gap0)*probe_index,
        width = width0,
        height = height0,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      plotText(
        label = "EA", x = x0+0.2, y = y0 + (height0 + gap0)*probe_index, rot = 0,
        fontsize = 10, just = "left",
        default.units = "cm"
      )
      
      annoSegments(
        x0 = unit(0, "npc"), y0 = 10,
        x1 = unit(1, "npc"), y1 = 10,
        plot = manhattanPlot2, default.units = "native",
        linecolor = "pink", lty = 2
      )
      
      
      ## Annotate genome label
      annoGenomeLabel(
        plot = manhattanPlot2, x = x0, y = y0 + (height0 + gap0)*probe_index + height0,
        fontsize = 8, scale = "Mb",
        just = c("left", "top"), default.units = "cm"
      )
      
      
      ## plot y-axis
      annoYaxis(
        plot = manhattanPlot2,
        at = y_axis_label,
        axisLine = TRUE, fontsize = 8
      )
      
###gwas_locus3
      probe_index=1+probe_index
      gwas_locus3$pch <- 21
      gwas_locus3$CHR = chr
      gwas_locus3$color = "#2b5799"
      gwas_locus3 <- gwas_locus3[,c("CHR","pos","p","SNP","color")]
      colnames(gwas_locus3) <- c("chrom","pos","p","snp","color")
      
      #----plot lncRNA eQTL plot ----------------------------------------------------------------
      yMAX=round(max(-log10(gwas_locus2$p),na.rm = T))*1.25
      devbuf1 = yMAX/4
      y_axis_label=round(seq(0,yMAX,devbuf1),0)
      
      ##https://phanstiellab.github.io/plotgardener/reference/plotManhattan.html
      manhattanPlot2 <- plotManhattan(
        data = gwas_locus3, 
        chrom = chr,
        chromstart = start,
        chromend = end,
        assembly = "hg19",
        fill = gwas_locus2$color,
        trans = "-log10",
        sigLine = TRUE, col = "darkgrey",
        baseline.color="black",
        cex=0.25,
        lty = 2, 
        range = c(0, yMAX),
        # leadSNP = list(
        #   snp = lead_snp,
        #   pch = 18,
        #   cex = 0.6,
        #   fill = "#9632b8",
        #   fontsize = 8
        # ),
        x = x0, y =y0 + (height0 + gap0)*probe_index,
        width = width0,
        height = height0,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      plotText(
        label = "IQ", x = x0+0.2, y = y0 + (height0 + gap0)*probe_index, rot = 0,
        fontsize = 10, just = "left",
        default.units = "cm"
      )
      
      annoSegments(
        x0 = unit(0, "npc"), y0 = 10,
        x1 = unit(1, "npc"), y1 = 10,
        plot = manhattanPlot2, default.units = "native",
        linecolor = "pink", lty = 2
      )
      
      
      ## Annotate genome label
      annoGenomeLabel(
        plot = manhattanPlot2, x = x0, y = y0 + (height0 + gap0)*probe_index + height0,
        fontsize = 8, scale = "Mb",
        just = c("left", "top"), default.units = "cm"
      )
      
      
      ## plot y-axis
      annoYaxis(
        plot = manhattanPlot2,
        at = y_axis_label,
        axisLine = TRUE, fontsize = 8
      )
      

                       
           
##
probe_list
probe_gene_list
probe_num=length(probe_list)
probe_gene_num=length(probe_gene_list)

  # ------------------------------------------------------------------------------
  # lncRNAeQTL -------------------------------------------------------------------
  
  if(probe_num>0  ){
    #probe_index=0
    for(probe in probe_list){
      probe_index=probe_index+1
      dat1<-dat[dat$type=="QTL",]
      gwas_locus2<-dat1[dat1$geneid==probe,]
      gwas_locus2$pch <- 21
      gwas_locus2$CHR = chr
      gwas_locus2$color = "#E8C32EFF"
      gwas_locus2 <- gwas_locus2[,c("CHR","pos","P","SNP","color")]
      colnames(gwas_locus2) <- c("chrom","pos","p","snp","color")
      
      #----plot lncRNA eQTL plot ----------------------------------------------------------------
      yMAX=round(max(-log10(gwas_locus2$p),na.rm = T))*1.25
      devbuf1 = yMAX/4
      y_axis_label=round(seq(0,yMAX,devbuf1),0)
      
      ##https://phanstiellab.github.io/plotgardener/reference/plotManhattan.html
      manhattanPlot2 <- plotManhattan(
        data = gwas_locus2, 
        chrom = chr,
        chromstart = start,
        chromend = end,
        assembly = "hg19",
        fill = gwas_locus2$color,
        trans = "-log10",
        sigLine = TRUE, col = "darkgrey",
        baseline.color="black",
        cex=0.25,
        lty = 2, 
        range = c(0, yMAX),
        # leadSNP = list(
        #   snp = lead_snp,
        #   pch = 18,
        #   cex = 0.6,
        #   fill = "#9632b8",
        #   fontsize = 8
        # ),
        x = x0, y =y0 + (height0 + gap0)*probe_index,
        width = width0,
        height = height0,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      plotText(
        label = probe, x = x0+0.2, y = y0 + (height0 + gap0)*probe_index, rot = 0,
        fontsize = 10, just = "left",
        default.units = "cm"
      )
      
      annoSegments(
        x0 = unit(0, "npc"), y0 = 10,
        x1 = unit(1, "npc"), y1 = 10,
        plot = manhattanPlot2, default.units = "native",
        linecolor = "pink", lty = 2
      )
      
      
      ## Annotate genome label
      annoGenomeLabel(
        plot = manhattanPlot2, x = x0, y = y0 + (height0 + gap0)*probe_index + height0,
        fontsize = 8, scale = "Mb",
        just = c("left", "top"), default.units = "cm"
      )
      
      
      ## plot y-axis
      annoYaxis(
        plot = manhattanPlot2,
        at = y_axis_label,
        axisLine = TRUE, fontsize = 8
      )
      

  


      
    }
  }



      
      
      

 # # ------------------------------------------------------------------------------
#  # eQTL -------------------------------------------------------------------
#   
#  if(probe_gene_num > 0 & probe_gene_list!="" ){
#    probe_gene_index=0
#    for(probe_gene in probe_gene_list){
#      probe_gene_index=probe_gene_index+1
#      dat1<-dat[dat$type=="QTL",]
#      gwas_locus3<-dat1[dat1$geneid==probe_gene,]
#      gwas_locus3$pch <- 21
#      gwas_locus3$CHR = chr
#      gwas_locus3$color = "#f06698"
#      gwas_locus3 <- gwas_locus3[,c("CHR","pos","P","SNP","color")]
#      colnames(gwas_locus3) <- c("chrom","pos","p","snp","color")
#      
#      
#      #----plot eQTL plot ----------------------------------------------------------------
#      yMAX=round(max(-log10(gwas_locus3$p),na.rm = T))*1.25
#      devbuf1 = yMAX/4
#      y_axis_label=round(seq(0,yMAX,devbuf1),0)
#      
#      manhattanPlot3 <- plotManhattan(
#        data = gwas_locus3, 
#        chrom = chr,
#        chromstart = start,
#        chromend = end,
#        assembly = "hg19",
#        fill = gwas_locus3$color,
#        trans = "-log10",
#        sigLine = TRUE, col = "darkgrey",
#        baseline.color="black",
#        cex=0.25,
#        lty = 2, 
#        range = c(0, yMAX),
#        # leadSNP = list(
#        #   snp = lead_snp,
#        #   pch = 18,
#        #   cex = 0.6,
#        #   fill = "#9632b8",
#        #   fontsize = 8
#        # ),
#        x = x0, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index,
#        width = width0,
#        height = height0,
#        just = c("left", "top"),
#        default.units = "cm"
#      )
#      
#      
#      plotText(
#        label = probe_gene, x = x0+0.2, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index, 
#        rot = 0,
#        fontsize = 10, just = "left",
#        default.units = "cm"
#      )
#      
#      annoSegments(
#        x0 = unit(0, "npc"), y0 = 10,
#        x1 = unit(1, "npc"), y1 = 10,
#        plot = manhattanPlot3, default.units = "native",
#        linecolor = "pink", lty = 2
#      )
#      
#      
#      ## Annotate genome label
#      annoGenomeLabel(
#        plot = manhattanPlot3, x = x0, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index + height0,
#        fontsize = 8, scale = "Mb",
#        just = c("left", "top"), default.units = "cm"
#      )
#      
#      
#      ## plot y-axis
#      annoYaxis(
#        plot = manhattanPlot3,
#        at = y_axis_label,
#        axisLine = TRUE, fontsize = 8)
#      
#      
#    }
#  }
#  
#  
#  ####GWAS finemap pip
#  
#     probe_index = probe_index+1;
#      pip$pch <- 21
#      pip$CHR = chr
#      pip$color = "#0000FF"
#      pip_locus <- pip[,c("CHR","start","PIP","SNP","color")]
#      colnames(pip_locus) <- c("chrom","pos","p","snp","color")
#      
#      
#      yMAX=round(max((pip_locus$p),na.rm = T))*1.25
#      devbuf1 = yMAX/4
#      y_axis_label=round(seq(0,yMAX,devbuf1),0)
#      
#      manhattanPlot4 <- plotManhattan(
#        data = pip_locus, 
#        chrom = chr,
#        chromstart = start,
#        chromend = end,
#        assembly = "hg19",
#        trans="",
#        fill = pip_locus$color,
#        sigLine = TRUE, col = "darkgrey",
#        baseline.color="black",
#        cex=0.25,
#        lty = 2, 
#        range = c(0, yMAX),
#        # leadSNP = list(
#        #   snp = lead_snp,
#        #   pch = 18,
#        #   cex = 0.6,
#        #   fill = "#9632b8",
#        #   fontsize = 8
#        # ),
#        x = x0, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index,
#        width = width0,
#        height = height0,
#        just = c("left", "top"),
#        default.units = "cm"
#      )
#      
#      
#      plotText(
#        label = "GWAS finemap", x = x0+0.2, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index, 
#        rot = 0,
#        fontsize = 10, just = "left",
#        default.units = "cm"
#      )
#      
#      annoSegments(
#        x0 = unit(0, "npc"), y0 = 10,
#        x1 = unit(1, "npc"), y1 = 10,
#        plot = manhattanPlot3, default.units = "native",
#        linecolor = "pink", lty = 2
#      )
#      
#      
#      ## Annotate genome label
#      annoGenomeLabel(
#        plot = manhattanPlot4, x = x0, y = y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index + height0,
#        fontsize = 8, scale = "Mb",
#        just = c("left", "top"), default.units = "cm"
#      )
#      
#      
#      ## plot y-axis
#      annoYaxis(
#        plot = manhattanPlot4,
#        at = y_axis_label,
#        axisLine = TRUE, fontsize = 8)
#      
#  
#  #----Plot y-axis label  --------------------------------------------------------
#  plotText(
#    label = "-log10(p-value)", 
#    x = x0-1.2, y = y0 + ( (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index + height0)/2, 
#    rot = 90,
#    fontsize = 10, just = "center",
#    default.units = "cm"
#  ) 
#  
#   plotText(
#    label = "PIP", 
#    x = x0-1.2, y = y0 + ( (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index + height0)-1, 
#    rot = 90,
#    fontsize = 10, just = "center",
#    default.units = "cm"
#  ) 
#  
#  
  
  
  
### gene model plot

########################## plot the gene labels #######################
  GeneRowNum = function(GENELIST) {
    
    BP_THRESH = 0.03; MAX_ROW = 3
    # get the start and end position
    GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
    START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
    
    STRLENGTH = nchar(as.character(GENELIST$GENE))
    MIDPOINT = (START1 + END1)/2
    
    
    
    START2 = MIDPOINT-bin0*STRLENGTH/10; END2 = MIDPOINT+bin0*STRLENGTH/10
    START = cbind(START1, START2); END = cbind(END1, END2);
    START = apply(START, 1, min); END = apply(END, 1, max)
    
    GENELIST = data.frame(GENELIST, START, END)
    GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
    START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
    
    # get the row index for each gene
    NBUF = dim(GENELIST)[1]
    ROWINDX = rep(1, NBUF)
    ROWEND = as.numeric(rep(0, MAX_ROW))
    MOVEFLAG = as.numeric(rep(0, NBUF))
    
    
    if(NBUF>1) {
      for( k in 2 : NBUF ) {
        ITERFLAG=FALSE
        if(START[k] < END[k-1]) {
          INDXBUF=ROWINDX[k-1]+1
        } else INDXBUF = 1
        
        if(INDXBUF>MAX_ROW) INDXBUF=1;
        
        REPTIME=0
        repeat{
          if( ROWEND[INDXBUF] > START[k] ) {
            ITERFLAG=FALSE
            INDXBUF=INDXBUF+1
            if(INDXBUF>MAX_ROW) INDXBUF = 1
          } else {
            ITERFLAG=TRUE
          }
          
          if(ITERFLAG) break;
          REPTIME = REPTIME+1
          if(REPTIME==MAX_ROW) break;
        }
        ROWINDX[k]=INDXBUF;
        
        if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
            | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
          MOVEFLAG[k] = 1
          SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
          MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
        }
        if(ROWEND[ROWINDX[k]]<END[k]) {
          ROWEND[ROWINDX[k]] = END[k]  }
      }
    }
    GENEROW = data.frame(as.character(GENELIST$GENE),
                         as.character(GENELIST$ORIENTATION),
                         as.numeric(GENELIST$GENESTART),
                         as.numeric(GENELIST$GENEEND),
                         ROWINDX, MOVEFLAG)
    colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
    return(GENEROW)
  }
  
  
  # -------------prepare the locus gene (lncRNA and gene code info) -----
  dat1<-dat[dat$type=="QTL",]
  dat1_lnc<-dat1[dat1$genetype=="lncRNA"  &  dat1$start > start & dat1$end < end,]
  locus_lncRNA_info=dat1_lnc[,c("chrom","start","end","geneid","genetype","strand","genename","genetype")]
  
  dat1_pcg<-dat1[ (dat1$genetype=="PCG" | dat1$genetype=="others")  &  dat1$start > start & dat1$end < end ,]
  locus_gene_info=dat1_pcg[,c("chrom","start","end","geneid","genetype","strand","genename","genetype")]
  
  
  
  
  
  # -------------- lncRNA labels ----
  gene_label_info=locus_lncRNA_info
  colnames(gene_label_info)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "CATEGORY","ORIENTATION","gene_name","gene_type");
  gene_label_info
  
  end
  start
  bin0=(end-start)/width0
  generow = GeneRowNum(gene_label_info)
  generow
  index=match(generow$GENE, gene_label_info$GENE, nomatch=0)
  generow$gene_name[which(index!=0)]=gene_label_info$gene_name[index]
  generow$gene_type[which(index!=0)]=gene_label_info$gene_type[index]
  num_row = max(as.numeric(generow$ROW));
  num_gene = dim(generow)[1]
  generow
  dist = 1
  yMax=y0 + (height0 + gap0)*probe_index + (height0 + gap0)*probe_gene_index + height0 + 1 + 0.8*num_row-3.5
  yMax
  
  for( k in 1 : num_row ) {
    generowbuf = generow[which(as.numeric(generow[,5])==k),]
    xstart = (as.numeric(generowbuf[,3]) - start)/bin0 + x0
    xend = (as.numeric(generowbuf[,4]) - start)/bin0 + x0
    snbuf = which(xend-xstart< 1e-3)
    if(length(snbuf)>0) {
      xstart[snbuf] = xstart[snbuf] - 0.0025
      xend[snbuf] = xend[snbuf] + 0.0025
    }
    
    xcenter = (xstart+xend)/2
    xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
    
    num_genebuf = dim(generowbuf)[1]
    for( l in 1 : num_genebuf ) {
      ofs=0.25
      if(l%%2==0) ofs=-0.25;
      m = num_row - k
      ypos = yMax - m*dist
      
      if(generowbuf[l,2]=="+"){
        gene_start=xstart[l]
        gene_end=xend[l]
      }else{
        gene_end=xstart[l]
        gene_start=xend[l]
      }
      
      if(as.character(generowbuf[l,"gene_name"]) %in% highlight_probe){
        plotSegments(
          x0 = gene_start, y0 = ypos, x1 =  gene_end, y1 = ypos,
          default.units = "cm",
          arrow = arrow(type = "open",length=unit(0.15,"cm")),
          lwd=1,lty=1,linecolor="orange",fill="red")
        
        genemove=0.01
        movebuf = as.numeric(generowbuf[l,6])*genemove
        
        plotText(
          label = substitute(genename, list(genename=as.character(generowbuf[l,"gene_name"]))), 
          x = xcenter[l]+movebuf, y = ypos + ofs, 
          rot = 0,
          fontsize = 6, just = "center",
          default.units = "cm",
          fontcolor="red",
          fontface="italic")
        
      }else{
      plotSegments(
        x0 = gene_start, y0 = ypos, x1 =  gene_end, y1 = ypos,
        default.units = "cm",
        arrow = arrow(type = "open",length=unit(0.15,"cm")),
        lwd=1,lty=1,linecolor="#E8C32EFF",fill="#E8C32EFF")
      
      genemove=0.01
      movebuf = as.numeric(generowbuf[l,6])*genemove
      
      plotText(
        label = substitute(genename, list(genename=as.character(generowbuf[l,"gene_name"]))), 
        x = xcenter[l]+movebuf, y = ypos + ofs, 
        rot = 0,
        fontsize = 6, just = "center",
        default.units = "cm",
        fontcolor="#E8C32EFF",
        fontface="italic")
      }
    }
  }
  
  
  # -------------- protein coding gene labels (non lncRNA labels, including others)----
  
  gene_label_info=locus_gene_info
  colnames(gene_label_info)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "CATEGORY","ORIENTATION","gene_name","gene_type");
  gene_label_info
 
 if(nrow(gene_label_info)!=0) {
  generow = GeneRowNum(gene_label_info)
  generow
  index=match(generow$GENE, gene_label_info$GENE, nomatch=0)
  generow$gene_name[which(index!=0)]=gene_label_info$gene_name[index]
  generow$gene_type[which(index!=0)]=gene_label_info$gene_type[index]
  num_row = max(as.numeric(generow$ROW));
  num_gene = dim(generow)[1]
  generow
  dist = 1
  yMax=yMax+ 1.0*num_row
  
  
  
  
  for( k in 1 : num_row ) {
    generowbuf = generow[which(as.numeric(generow[,5])==k),]
    xstart = (as.numeric(generowbuf[,3]) - start)/bin0 + x0
    xend = (as.numeric(generowbuf[,4]) - start)/bin0 + x0
    snbuf = which(xend-xstart< 1e-3)
    if(length(snbuf)>0) {
      xstart[snbuf] = xstart[snbuf] - 0.0025
      xend[snbuf] = xend[snbuf] + 0.0025
    }
    
    xcenter = (xstart+xend)/2
    xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
    
    num_genebuf = dim(generowbuf)[1]
    for( l in 1 : num_genebuf ) {
      ofs=0.25
      if(l%%2==0) ofs=-0.25;
      m = num_row - k
      ypos = yMax - m*dist
      
      if(generowbuf[l,2]=="+"){
        gene_start=xstart[l]
        gene_end=xend[l]
      }else{
        gene_end=xstart[l]
        gene_start=xend[l]
      }
      
      if(as.character(generowbuf[l,"gene_name"]) %in% highlight_probe){
        plotSegments(
          x0 = gene_start, y0 = ypos, x1 =  gene_end, y1 = ypos,
          default.units = "cm",
          arrow = arrow(type = "open",length=unit(0.15,"cm")),
          lwd=1,lty=1,linecolor="red",fill="red")
        
        genemove=0.01
        movebuf = as.numeric(generowbuf[l,6])*genemove
        
        plotText(
          label = substitute(genename, list(genename=as.character(generowbuf[l,"gene_name"]))), 
          x = xcenter[l]+movebuf, y = ypos + ofs, 
          rot = 0,
          fontsize = 6, just = "center",
          default.units = "cm",
          fontcolor="red",
          fontface="italic")
      }else{
      plotSegments(
        x0 = gene_start, y0 = ypos, x1 =  gene_end, y1 = ypos,
        default.units = "cm",
        arrow = arrow(type = "open",length=unit(0.15,"cm")),
        lwd=1,lty=1,linecolor="#f06698",fill="#f06698")
      
      genemove=0.01
      movebuf = as.numeric(generowbuf[l,6])*genemove
      
      plotText(
        label = substitute(genename, list(genename=as.character(generowbuf[l,"gene_name"]))), 
        x = xcenter[l]+movebuf, y = ypos + ofs, 
        rot = 0,
        fontsize = 6, just = "center",
        default.units = "cm",
        fontcolor="#f06698",
        fontface="italic")
      }}
    }
}
    
#
#######Hi-C bedpe
#		NDNAloops_pairs <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/Promoter-anchored_chromatin_loops.bed",sep="\t",header=T);	
#    NDNAloops_pairs$chr2=NDNAloops_pairs$chr1
#    NDNAloops_pairs<-NDNAloops_pairs[,c("chr1","start1","end1","chr2","start2","end2")]  
#  
#    
#    NDNAloops_pairs1 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/hNPC.RICseq.bedpe",sep="\t",header=T);
#    NDNAloops_pairs2 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GZ.loops.bedpe",sep="\t",header=T);
#    NDNAloops_pairs3 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/CP.loops.bedpe",sep="\t",header=T);
#    NDNAloops_pairs4 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/cortex_fetal.loops.bedpe",sep="\t",header=T);
#    NDNAloops_pairs5 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/cortex_adult.loops.bedpe",sep="\t",header=T);
#    NDNAloops_pairs6 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GSM3598048_motor.bedpe",sep="\t",header=T);
#    NDNAloops_pairs7 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GSM3106832_cortical.bedpe",sep="\t",header=T);
#    NDNAloops_pairs8 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GSM3598046_hippocampal.bedpe",sep="\t",header=T);
#    NDNAloops_pairs9 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GSM3598051_astrocyte.bedpe",sep="\t",header=T);
#    NDNAloops_pairs10 <-read.table(file="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/brain-cell-type-peak-files-master/PLACseq/NeuN.5k.2.peaks.bedpe",sep="\t",header=T);
#    
#    params <- pgParams(
#    chrom = chr, chromstart = start, chromend = end,
#    assembly = "hg19",
#    x = x0, just = c("left", "top"),
#    width = 2, length = 2, default.units = "cm"
#		)
#   
#  
#    plotPairsArches(
#    data = NDNAloops_pairs10,
#    params = params,width =width0, height = height0*0.5,
#    y = 19.5 ,
#    fill = "purple", linecolor = "purple", flip = TRUE
#		)   
#		
#		plotText(
#        label = "NeuN.5k", x = x0, y = 19.5 , 
#        rot = 0,
#        fontsize = 12, just = "left",
#        default.units = "cm"
#      )
#      
#     
#   
#   
#   ##params <- pgParams(
###    chrom = chr, chromstart = start, chromend = end,
###    assembly = "hg19",
###    x = x0, just = c("left", "top"),
###    width = 2, length = 2, default.units = "cm"
###		)
###   
###  
###    plotPairsArches(
###    data = NDNAloops_pairs9,
###    params = params,width =width0, height = height0*0.5,
###    y = 25 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "GSM3598051_astrocyte", x = x0, y = 25 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
###      
###      
### 
###    plotPairsArches(
###    data = NDNAloops_pairs8,
###    params = params,width =width0, height = height0*0.5,
###    y = 26 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "GSM3598046_hippocampal", x = x0, y = 26 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )  
###      
###      
###    
###     plotPairsArches(
###    data = NDNAloops_pairs7,
###    params = params,width =width0, height = height0*0.5,
###    y = 27 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "GSM3106832_cortical", x = x0, y = 27 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
###      
###
###
###     
###      plotPairsArches(
###    data = NDNAloops_pairs6,
###    params = params,width =width0, height = height0*0.5,
###    y = 28 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "GSM3598048_motor", x = x0, y = 28 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
###      
###      
###      
###       plotPairsArches(
###    data = NDNAloops_pairs5,
###    params = params,width =width0, height = height0*0.5,
###    y = 29 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "cortex_adult", x = x0, y = 29 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
###      
###      
###      
###       plotPairsArches(
###    data = NDNAloops_pairs4,
###    params = params,width =width0, height = height0*0.5,
###    y = 30 ,
###    fill = "purple", linecolor = "purple", flip = TRUE
###		)   
###		
###		plotText(
###        label = "cortex_fetal", x = x0, y = 16.5 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
#      
#      
#
#		
### bigwig
#
#		Astrocyte="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/Astrocyte.bw";		
#		Astrocyte_bw <- readBigwig(file = Astrocyte, chrom = chr,  chromstart =start,   chromend = end)		
#		GABA="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/ITL23_1.cov.bw";	
#			
#			
#		GABA="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/brain-cell-type-peak-files-master/human_NEUNnuclei_H3K4me3_epilepsy_hg19.ucsc.bigWig";
#		GABA="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/brain-cell-type-peak-files-master/human_NEUNnuclei_H3K27ac_epilepsy_pooled_hg19.ucsc.bigWig";
#		GABA="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/GABA.bw";	
#		GABA_bw <- readBigwig(file = GABA, chrom = chr,  chromstart =start,   chromend = end)
#		Glutamatergic="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/Glutamatergic.bw";		
#		Glutamatergic_bw <- readBigwig(file = Glutamatergic, chrom = chr,  chromstart =start,   chromend = end)
#		Oligodendrocyte="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/Oligodendrocyte.bw";		
#		Oligodendrocyte_bw <- readBigwig(file = Oligodendrocyte, chrom = chr,  chromstart =start,   chromend = end)
#		Microglia="/storage/yangjianLab/chenli/seq/PsychENCODE/metaeQTL/SMR/locos_plot/Microglia.bw";		
#		Microglia_bw <- readBigwig(file = Microglia, chrom = chr,  chromstart =start,   chromend = end)
#		
#		library("RColorBrewer")
#		signalList <- list(Astrocyte_bw, GABA_bw,Glutamatergic_bw, Oligodendrocyte_bw,Microglia_bw)
#    multisignal <- plotMultiSignal(signalList, chrom = chr,
#                               chromstart = start, chromend = end,
#                               linecolor = c(brewer.pal(n = 9,"YlGnBu")[3],
#                               							 brewer.pal(n = 9,"YlGnBu")[4],
#                                             brewer.pal(n = 9,"YlGnBu")[5],
#                                             brewer.pal(n = 9,"YlGnBu")[6],
#                                             brewer.pal(n = 9,"YlGnBu")[7]),
#                               label = c("Astrocyte", "GABA",
#                                         "Glutamatergic", "Oligodendrocyte","Microglia"),
#                               assembly = "hg19",
#                               x = x0, y = 21,
#                               width =width0, height = height0*3,
#                               default.units = "cm",
#                               gapdistance = 0.1)
#		
#		##region <- pgParams(
###    chrom = chr,
###    chromstart = start, chromend = end,
###    assembly = "hg19",
###    range = c(0, 45)
###		)
###
###		
###		signal1 <- plotSignal(
###		    data = Astrocyte_bw, params = region,
###		    x = x0, y = 19  , width =width0, height = height0-1,
###		    just = c("left", "top"),  linecolor = "#aec7e8", default.units = "cm"
###		)
###   
###    plotText(
###        label = "Astrocyte", x = x0, y = 19 , 
###        rot = 0,
###        fontsize = 10, just = "left",
###        default.units = "cm"
###      )
#      
#		
	 			
    
			    ###plotGenomeLabel
			    plotGenomeLabel(
			      chrom = chr, chromstart = start, chromend = end,
			      assembly = "hg19",
			      x = x0, y = 19, length = width0, scale = "Mb",
			      just = c("left", "top"), default.units = "cm",
			      fontsize = 8
			    )
			    
			    
##			    ## Set genomic coordinates
##			    # https://phanstiellab.github.io/plotgardener/reference/plotGenes.html
##			paramsbig <- pgParams(
##			    chrom = chr,
##			    chromstart = start, chromend = end,
##			    assembly = "hg19", width = 7
##			)
##			## Set colors
##			cols <- c("#41B6C4", "#225EA8")
##
##			## Plot genes big
##			genesPlot <- plotGenes(
##			    params = paramsbig, fill = cols,
##			    fontcolor = cols,
##			    x = x0, y = 25, width =width0, height = height0,
##			    just = c("left", "top"),
##			    default.units = "cm"
##			)
##    

  pageGuideHide()	
  
}


pdf("noPCG-BrainMetalnc.28095-CP.R.pdf",  width=10,height=13)  
gwasLocusPlot()
dev.off()


	