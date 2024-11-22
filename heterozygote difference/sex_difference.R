library(sqldf)
# chrX
data <- read.table("heterozygote_difference_chrX.txt",header = F)
data <- data[,c('V1','V2','V3','V4','V5','V6','V10','V11','V12')]
names(data) <- c('chr','position','allele_freq','male_geno0','male_heter','male_geno1','female_geno0',
                 'female_heter','female_geno1')

plotXdata=function(winsize,daf1,daf2,samplesize)
{
  tmp = sqldf(sprintf("select floor(position/%d) as window,sum(male_geno0*male_geno1) as malediff,
                      sum((male_geno0+male_geno1)*(male_geno0+male_geno1-1))/2 as maletimes,
                      sum(female_heter) as femalediff,sum(female_geno0+female_heter+female_geno1) as femaletimes from data 
                      where allele_freq>=%d and allele_freq<=%d and male_heter=0 and
                      male_geno0+male_geno1+female_geno0+female_heter+female_geno1>=%d group by window",winsize,daf1,daf2,samplesize),drv="SQLite")
  tmp$malediff/tmp$maletimes
  tmp$maleratio = tmp$malediff/tmp$maletimes
  tmp$femaleratio = tmp$femalediff/tmp$femaletimes
  tmp$ratio <- (tmp$femaleratio-tmp$maleratio)*2/(tmp$femaleratio+tmp$maleratio)
  tmp
}
#SNPs with MAF more than 0.3
plotXdata(100000,6,9,8)
#SNPs with MAF less than 0.3
plotXdata(100000,0,6,8)

# Autosomes
data <- read.table("heterozygote_difference_Autosomes.txt",header = F)
data <- data[,c('V1','V2','V3','V4','V5','V6','V10','V11','V12')]
names(data) <- c('chr','position','allele_freq','male_geno0','male_heter','male_geno1','female_geno0',
                 'female_heter','female_geno1')

plotAdata=function(winsize,daf1,daf2,samplesize)
{
  tmp = sqldf(sprintf("select floor(position/%d) as window,
  sum(male_geno0*male_heter*2 + male_heter*male_geno1*2 + male_heter*(male_heter-1) + male_geno0*male_geno1*4) as malediff,
  sum((male_geno0+male_heter+male_geno1)*(male_geno0+male_heter+male_geno1-1))*2 as maletimes,
  sum(female_heter) as femalediff,
  sum(female_geno0+female_heter+female_geno1) as femaletimes 
  from data where allele_freq>=%d and allele_freq<=%d and 
  male_geno0+male_heter+male_geno1+female_geno0+female_heter+female_geno1>=%d group by window",winsize,daf1,daf2,samplesize),drv="SQLite")
  tmp$maleratio = tmp$malediff/tmp$maletimes
  tmp$femaleratio = tmp$femalediff/tmp$femaletimes
  tmp$ratio = (tmp$femaleratio-tmp$maleratio)*2/(tmp$femaleratio+tmp$maleratio)
  tmp
}

#SNPs with MAF more than 0.3
plotAdata(100000,6,9,8)
#SNPs with MAF less than 0.3
plotAdata(100000,0,6,8)


