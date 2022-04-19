===============
QTL_Rice_Cold_Tolerance
===============

:Author: Michael Hall
:Date:   4/13/2022

QTL-Rice-Cold-Tolerance
==========

QTLseqr is an R package for QTL mapping using NGS Bulk Segregant
Analysis.

QTLseqr is still under development and is offered with out any
guarantee.

**For more detailed instructions please read the vignette**\ `here <https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf>`__
---------------------------------------------------------------------------------------------------------------------------------------------

For updates read the `NEWS.md <https://github.com/bmansfeld/QTLseqr/blob/master/NEWS.md>`__
-------------------------------------------------------------------------------------------

Installation
============

.. raw:: html

   <!-- You can install and update QTLseqr by using our [drat](http://dirk.eddelbuettel.com/code/drat.html) repository hosted on our github page: -->

.. raw:: html

   <!-- ```{r drat-install, eval = FALSE} -->

.. raw:: html

   <!-- install.packages("QTLseqr", repos = "http://bmansfeld.github.io/drat") -->

.. raw:: html

   <!-- ``` -->

.. raw:: html

   <!-- OR You can install QTLseqr from github with: -->

You can install QTLseqr from github with:

.. code:: r

   # install devtools first to download packages from github
   install.packages("devtools")

   # use devtools to install QTLseqr
   devtools::install_github("PBGLMichaelHall/QTLseqr")

**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).

**If you use QTLseqr in published research, please cite:**

   Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
   analysis with next-generation sequencing *The Plant Genome*
   `doi:10.3835/plantgenome2018.01.0006 <https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006>`__

We also recommend citing the paper for the corresponding method you work
with.

QTL-seq method:

   Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
   C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
   L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
   quantitative trait loci in rice by whole genome resequencing of DNA
   from two bulked populations. *Plant J*, 74: 174–183.
   `doi:10.1111/tpj.12105 <https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105>`__

G prime method:

   Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
   Segregant Analysis Using Next Generation Sequencing. *PLOS
   Computational Biology* 7(11): e1002255.
   `doi.org/10.1371/journal.pcbi.1002255 <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255>`__

Abstract
--------

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is
efficient in detecting quantitative trait loci (QTL). Despite the
popularity of NGS-BSA and the R statistical platform, no R packages are
currently available for NGS-BSA. We present QTLseqr, an R package for
NGS-BSA that identifies QTL using two statistical approaches: QTL-seq
and G’. These approaches use a simulation method and a tricube smoothed
G statistic, respectively, to identify and assess statistical
significance of QTL. QTLseqr, can import and filter SNP data, calculate
SNP distributions, relative allele frequencies, G’ values, and
log10(p-values), enabling identification and plotting of QTL.

Examples:
=========

Load/install libraries
======================

.. code:: r 
   
   install.packages(“tinytex”) 
   install.packages(“vcfR”) 
   install.packages(“tidyr”) 
   install.packages(“ggplot2”)
   devtools::install_github(“PBGLMichaelHall/QTLseqr”,force = TRUE)   
   library(QTLseqr) 
   library(tinytex) 
   library(vcfR) 
   library(tidyr)
   library(ggplot2)

::

   # Set the Working Directory to where VCF file is stored in file system

.. code:: r 

setwd("/home/michael/Desktop/RiceCold2")

Vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTLseqR functions will only take SNPS, ie, length of REF and ALT== 1
==================================================================================================================================================================================

.. code:: r

vcf <- read.vcfR(file = "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf.gz"

.. figure:: ../images/1.png
   :alt: 

.. code:: r


   #Convert to tidy data frame
   VCF_TIDY <- vcfR2tidy(vcf)


Call the Parser
===============

.. code:: r

   QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "D2_F2_tt",LowBulk = "D2_F2_TT", filename = Hall)

.. figure:: ../images/3.png
   :alt: 


Invoke unique command to extract Sample names reverse comapatible to the VCF
============================================================================

.. code:: r

   unique(VCF_TIDY$gt$Indiv)


.. code:: r

   #Set High bulk and Low bulk sample names and parser generated file name
   #The file name is generated from the QTLParser_1_MH function in line 119

   HighBulk <- "ET-pool-385"
   LowBulk <- "ES-pool-430"
   file <- "Hall.csv"

   #Choose which chromosomes/contigs will be included in the analysis,

   Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")


   df <-
     importFromTable(
       file = file,
       highBulk = HighBulk,
       lowBulk = LowBulk,
       chromList = Chroms
     ) 



.. code:: r

#plot histograms associated with filtering arguments such as mamximum and minumum Total Depths and reference Allele Frequency to determine cut off     values 
   ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)
   ggsave(filename = "Depth_Histogram.png",plot=last_plot())

.. figure:: ../images/8.png
   :alt: 

.. code:: r

   ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))
   ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

.. figure:: ../images/9.png
   :alt: 



.. code:: r

   #Filter SNPs based on some criteria 
   df_filt <- filterSNPs( SNPset = df,
   refAlleleFreq = 0.20, minTotalDepth = 100, maxTotalDepth = 400,
   minSampleDepth = 40, 
   # minGQ = 0 )

.. figure:: ../images/10.png
   :alt: 

.. code:: r

   #Run G' analysis
   df_filt<-runGprimeAnalysis_MH(
     SNPset = df_filt,
     windowSize = 1e6,
     outlierFilter = "deltaSNP",
     filterThreshold = 0.1)

.. figure:: ../images/11.png
   :alt: 

 

G’ Distribution Plot
====================

.. code:: r

   #The plot reveals a skewed G Prime statistic with a really small variance. Perhaps it is due to the small number of variants called.
   #In addition, Hampels outlier filter in the second argument, can also be changed to "deltaSNP"
   
   plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")

.. figure:: ../images/12.png
   :alt: 


.. code:: r

   #We can see raw data before and after our filtering step
   
   plotGprimeDist_MH(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)

.. figure:: ../images/13.png
   :alt: 

.. code:: r
   

   #Run QTLseq analysis
   df_filt2 <- runQTLseqAnalysis_MH(
     SNPset = df_filt,
     windowSize = 1e6,
     popStruc = "F2",
     bulkSize = c(430, 385),
     replications = 10000,
     intervals = c(95, 99)
   )

.. figure:: ../images/14.png
   :alt: 



Plot G Statistic Distribution
=============================

.. code:: r

   hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")

.. figure:: ../images/15.png
   :alt:



.. code:: r

   #Plot Snps as a function of chromosome and position values
   
   plotQTLStats(SNPset = df_filt2, var = "nSNPs")
   ggsave(filename = "nSNPs.png",plot = last_plot())

.. figure:: ../images/16.png
   :alt: 

 

.. code:: r

   #Using QTLStats funciton plot Gprime Statistic with False Discovery Rate Threhshold as a third argument boolean operator as TRUE. The q value is used as FDR threshold null value is 0.05%.
   
   plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
   ggsave(filename = "GPrime.png",plot = last_plot())

.. figure:: ../images/17.png
   :alt: 

  

.. code:: r

   #Again using plotQTLStats change second argument varaible to deltaSNP and plot.
   
   plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
   ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

.. figure:: ../images/18.png
   :alt: 

 

.. code:: r

   #Finally with plotQTLStats plot negLog10Pval
   
   plotQTLStats(SNPset = df_filt, var = "negLog10Pval",plotThreshold = TRUE,q=0.15)
   ggsave(filename = "negLog10Pval.png",plot = last_plot())

.. figure:: ../images/19.png
   :alt: 

   

.. code:: r

   #Add subset argument to focus on particular chromosomes one, three, four, and six.
   #The reason is due to signficant QTL regions
   plotQTLStats(SNPset = df_filt, var = "Gprime",plotThreshold = TRUE,q=0.05,subset = c("NC_029256.1","NC_029258.1","NC_029259.1","NC_029261.1"))

.. figure:: ../images/20.png
   :alt:



Use RMVP package to view SNPs on chromosomes/contigs
====================================================

.. code:: r

   #install.packages("rMVP")
   library(rMVP)
   sample<-"Semi_Dwarfism_in_Sorghum"
   pathtosample <- "/home/michael/Desktop/QTLseqr/extdata/subset_freebayes_D2.filtered.vcf.gz"
   out<- paste0("mvp.",sample,".vcf")
   memo<-paste0(sample)
   dffile<-paste0("mvp.",sample,".vcf.geno.map")

   message("Making MVP data S1")
   MVP.Data(fileVCF=pathtosample,
         #filePhe="Phenotype.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out=out)
         
   message("Reading MVP Data S1")
   df <- read.table(file = dffile, header=TRUE)
   message("Making SNP Density Plots")
   MVP.Report.Density(df[,c(1:3)], bin.size = 1000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


.. figure:: ../images/21.png
   :alt: 

 

Export summary CSV
==================

.. code:: r

   QTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

Preview the Summary QTL
=======================

.. figure:: ../images/22.png
   :alt: 

 

.. code:: r

   #Use the function to plot allele frequencies per chromosome
   #Second argument size specifes size of scalar factor on nSNPs and if you have a relatively small SNP set .001 is a good startin point otherwise set to 1
   Obs_Allele_Freq(SNPSet = df_filt)

.. figure:: ../images/23.png
   :alt:
   
.. figure:: ../images/233.png
   :alt:

 

.. code:: r

   ##Use the function to investigate chromosomal region of interest
   Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = "Chr04", threshold = .90)

.. figure:: ../images/24.png
   :alt: 

.. code:: r

   setwd("/home/michael/Desktop/QTLseqr/extdata")
   # Theory and Analytical Framework of Sampling from BSA
   par(mfrow=c(1,1))
   # Define Ranges of Success
   # Sample Size from High Bulk sn = 385
   success <- 0:770
   # The Difference between realized and Expected Frequencies 
   # ns : Sample Size taken from Low Bulk
   # 2(ns)p1_star ~ Binomial(2(ns),p1)
   # p1 Expected Frequencies
   # Expected Frequencies:
   # E(n1) = E(n2) = E(n3) = E(n4) = C/2 = 110
   # We prefer for accuracy to have ns >> C >> 1
   plot(success, dbinom(success, size = 770, prob = .50), type = "h",main="Binomial Sampling from Diploid Orgainism from High Bulk",xlab="2(ns)(p1_STAR)",ylab="Density")

.. figure:: ../images/25.png
   :alt: 


.. code:: r


   # ns : Sample Size from High Bulk
   # 2(ns)p2_star ~ Binomial(2(ns),p2)
   # p2 Expected Frequencies
   success <- 0:860
   plot(success, dbinom(success, size = 860, prob = 0.5), type = "h",main="Binomial Sampling from Diploid Organism from Low Bulk",xlab="2(n2)(p2_STAR)",ylab="Density")

.. figure:: ../images/26.png
   :alt: 

 

.. code:: r



   par(mfrow=c(1,1))
   #Define Ranges of Success (Allele Frequencies High and Low)
   success <- 0:100
   #n1|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p1_STAR)), type = 'h',main="n1|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n1|(n3/n1+n3)",ylab="Prob")

.. figure:: ../images/27.png
   :alt: 

 

.. code:: r

# Filter outliers
TT <- TT %>% filter(AD_REF. <= 500)

hist(TT$AD_REF., probability = FALSE,main="Histogram of Actually Realized n1 Values",xlab="n1",breaks = "Sturges")



.. figure:: ../images/28.png
   :alt: 

  

.. code:: r

   #n2|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*(1-p2_STAR)), type='h', main="n2|p2_STAR ~ Poisson(C[[1-p2_STAR])",xlab="n2|(n4/n2+n4)",ylab="Prob")

.. figure:: ../images/29.png
   :alt: 



.. code:: r

tt <- tt %>% filter(AD_REF. <= 500)
hist(tt$AD_REF., probability = TRUE, main = "Histogram of Actually Realized n2 Values",xlab="n2")

.. figure:: ../images/30.png
   :alt: 

 

.. code:: r

   #n3|p1_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p1_STAR),type='h',main="n3|p1_STAR ~ Poisson(C[1-p1_STAR])",xlab="n3|(n3/n1+n3)",ylab="Prob")

.. figure:: ../images/31.png
   :alt: 


.. code:: r


TT <- TT %>% filter(AD_ALT. <= 300)
hist(TT$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n3 Values",xlab="n3")


.. figure:: ../images/32.png
   :alt:

.. code:: r

   #n4|p2_star ~ Poisson(lambda)
   plot(success, dpois(success, lambda = C*p2_STAR), type = 'h',main="n4|p2_STAR ~ Poisson(C[1-p2_STAR])",xlab="n4|n4/(n2+n4)",ylab="Prob")

.. figure:: ../images/33.png
   :alt: 

.. code:: r

   hist(tt$AD_ALT., probability = TRUE, main="Histogram of Acutally Realized n4 Values",xlab="n4")

.. figure:: ../images/34.png
   :alt: 

.. code:: r

   #Assuming average sequencing coverage (C) expected values for n1,n2,n3,n4
   E(n1) = E(n2) = E(n3) = E(n4) = C/2 = 35




# Read in the csv file from High bulk tt
tt<-read.table(file = "ET-pool-385.csv",header = TRUE,sep = ",")
# Calculate average Coverage per SNP site
mean(tt$DP)
# Find REalized frequencies
p1_STAR <- sum(tt$AD_ALT.) / sum(tt$DP)

# Read in the csv file from Low Bulk TT
TT<-read.table(file ="ES-pool-430.csv",header = TRUE,sep=",")
# Calculate average Coverage per SNP sit
mean(TT$DP)
# Find Realized frequencies
p2_STAR <- sum(TT$AD_ALT.) / sum(TT$DP)
# Take the average of the Averages
C <-(mean(tt$DP)+mean(TT$DP))/2
C<-round(C,0)
#Average Coverage
70
C/2 = 35



p2 >> p1 QTL is present
=======================

However, ns >> C >> 1 is TRUE 
=================================






