rm(list=ls())

library('coloc')
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]
QTL = as.numeric(args[3])
GWAS_case = as.numeric(args[4])
GWAS_control = as.numeric(args[5])

# Rscript /storage/yangjianLab/chenli/scripts/GWASscripts/ColocAnalysisCC.R  SCZ-CMC.lncRNAeQTL.chr4_COLOC_INPUT.txt SCZ-CMC.lncRNAeQTL.chr4.coloc  2865  76755  243649
# input = "SCZ-CMC.lncRNAeQTL.chr4_COLOC_INPUT.txt"
# output = "SCZ-CMC.lncRNAeQTL.chr4.coloc"
# QTL = 2865
# GWAS_case = 76755
# GWAS_control = 243649



Dat <- read.table(input,header=T,stringsAsFactors = F)
Dat <- as.data.frame(Dat)
names(Dat) <- c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_Chr","Probe_bp","Gene","Orientation","b","se","p","freq2","beta2","se2","p2","N2","minorAF")
#Dat$varbeta1 <- (Dat$se)*(Dat$se)
#Dat$varbeta2 <- (Dat$se2)*(Dat$se2)

calcu_std_b_se<-function(z,p,n)
{
	std_b_hat=z/sqrt(2*p*(1-p)*(n+z^2))
	std_se=1/sqrt(2*p*(1-p)*(n+z^2))
	res<-data.frame(std_b_hat,std_se);
	return(res)
}

##https://www.rdocumentation.org/packages/coloc/versions/5.1.0/topics/check_dataset  
Probes <- unique(Dat$Probe)
RES <- c("Probe","Probe_Chr","Probe_bp","eGene_topSNP","eGene_topSNP_p","GWAS_topSNP","GWAS_topSNP_p","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
for (i in (1:length(Probes)))
{
		TEST <- Dat[Dat$Probe==Probes[i],]
		TEST$minorAF=ifelse(TEST$Freq>0.5,1-TEST$Freq,TEST$Freq)
		
		
		##QTLSdata
		#std_eff1<-calcu_std_b_se(TEST$b/TEST$se,TEST$Freq,QTL)
		#TEST$b<-std_eff1[,1]
		#TEST$se<-std_eff1[,2]
		TEST$varbeta1 <- (TEST$se)*(TEST$se)
		D1=TEST[c("b","varbeta1","SNP","BP","minorAF")]
		colnames(D1) <- c("beta","varbeta","snp","position","MAF")
		D1$type <- "quant"
		D1$N <- QTL
		
		##GWASdata
		TEST$freq2=TEST$Freq;
		#std_eff1<-calcu_std_b_se(TEST$beta2/TEST$se2,TEST$freq2,GWAS_case+GWAS_control)
		#TEST$beta2<-std_eff1[,1]
		#TEST$se2<-std_eff1[,2]
		TEST$varbeta2 <- (TEST$se2)*(TEST$se2)
		
		D2=TEST[,c("minorAF","SNP","BP","p2")]
		colnames(D2) <- c("MAF","snp","position","pvalues")
		D2$type <- "cc"
		D2$s <- GWAS_case/GWAS_control
		D2$N <- (GWAS_case+GWAS_control)
		res1 <- c(TEST$Gene[1],TEST$Probe_Chr[1],TEST$Probe_bp[1],TEST$SNP[which.min(TEST$p)],min(TEST$p),TEST$SNP[which.min(TEST$p2)],min(TEST$p2))
		#my.res <- coloc.abf(dataset1=list(beta=TEST$b, MAF=TEST$minorAF, varbeta=TEST$varbeta1, N=QTL,type="quant"),
		#                    dataset2=list(pvalues=TEST$p2,MAF=TEST$minorAF, N=GWAS_case+GWAS_control,s=GWAS_case/(GWAS_case+GWAS_control),type="cc"))
		my.res <- coloc.abf(dataset1=list(beta=TEST$b, MAF=TEST$minorAF, varbeta=TEST$varbeta1, N=QTL,type="quant"),
		                    dataset2=list(beta=TEST$beta2, MAF=TEST$minorAF, varbeta=TEST$varbeta2, N=GWAS_case+GWAS_control,type="cc"))
		res2 <- as.data.frame(my.res$summary)
		res3 <- c(res1,res2$`my.res$summary`)
		RES <- rbind(RES,res3)
}

RES2 <- as.data.frame(RES[c(2:nrow(RES)),])
names(RES2) <- RES[1,]
##RESSIG <- RES2[as.numeric(RES2$PP.H4.abf)>0.8,]

write.table(file=paste0(output,".txt"),RES2,col.names = TRUE, row.names = FALSE, quote = FALSE,  sep = "\t")

