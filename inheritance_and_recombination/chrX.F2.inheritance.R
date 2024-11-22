# This R script is written to determine the genomic transmission of X chromosome from F0 to F1 and F2.
# It takes shapeit2's output as an input.
# It can be used to determine the paternal and maternal alleles in each of F2.
# Slight modifications can be made and then it can apply to autosomal studies.
# Author: Hai-Bing Xie
# Email:  xiehb@mail.kiz.ac.cn
#

require(sqldf)

setwd("E:\\xiehb_sync\\F2\\data\\sscrofa11.data")

#load data
haplotype = read.table("chrX.phased.duohmm.haps")
sample = read.table("chrX.phased.duohmm.sample",skip=1)

f2id = as.data.frame(as.matrix(sample[-c(1:65),2]))
f2male = f2id[f2id%%2==1]
f2female = f2id[f2id%%2==0]

# pedigree saves the family information
# the values are the sample id
#-------------------------------------------------------------------------  
# table format (7 columns):
#-------------------------------------------------------------------------  
#       F2   F1male   F0male   F0female   F1female  F0male  F0female
#     ^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^
#     f2id          family1                       family2
#-------------------------------------------------------------------------  

f1father = sqldf("select b.V4 from f2id a, sample b where a.V1=b.V2",drv="SQLite")
f1mother = sqldf("select b.V5 from f2id a, sample b where a.V1=b.V2",drv="SQLite")
family1 = sqldf("select b.V2,b.V4,b.V5 from f1father a, sample b where a.V4=b.V2",drv="SQLite")
family2 = sqldf("select b.V2,b.V4,b.V5 from f1mother a, sample b where a.V5=b.V2",drv="SQLite")
pedigree=cbind(f2id,family1,family2)

pedigree_tmp = as.data.frame(matrix(as.vector(as.matrix(pedigree)),byrow=T,ncol=1))

allsample = as.vector(sample[[2]])
allsample = as.data.frame(cbind(1:length(allsample),allsample))

allsamplepos = sqldf("select a.V1 from pedigree_tmp b left join allsample a on b.V1 = a.allsample",drv="SQLite")
allsamplepos = as.vector(as.matrix(allsamplepos))

#pedigree_info saves the column # for each individual
pedigree_info = matrix(allsamplepos,byrow=F,nrow=nrow(pedigree))

# pedigree_haps saves the column # for chromosomes in the haplotype (data frame)
pedigree_haps1 = (pedigree_info-1)*2+4
pedigree_haps2 = (pedigree_info-1)*2+5
pedigree_haps = cbind(pedigree_haps1,pedigree_haps2)[,c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]

# add SNP coordinates column # to pedigree_haps 
pedigree_haps = cbind(3,pedigree_haps)
pedigree_haps = cbind(1,pedigree_haps)
pedigree_haps[is.na(pedigree_haps)]=1

#determine the difference between two haplotypes
# m and n are numbers indicating the column # in the family
GetScore = function(m,n)
{
  return(abs(pedigree_haplotype[[m]]-pedigree_haplotype[[n]]))
}

# To figure out which part of the haplotype in a chromosome (x) is inherited from what part of the haplotypes of parental chromosomes y and z.
# The primary role of this function is to identify putative recombination events
# The return value is the haplotype of x coded by y and z
# x, y, z is a column # specifying different haplotypes in the data frame for a specific family
StepHap = function(x,y,z,n)
{
  mismatchY = GetScore(x,y)
  mismatchZ = GetScore(x,z)

  cbind(mismatchY, mismatchZ)[1:100,] 
  finalseq = 0
  localsites = n
  
  finalseq=rep(0,localsites)
  lastparent = 0
  lastblockend = 0
  
  while(length(finalseq[finalseq==0])>0)
  {
    vectorY = which(mismatchY[(lastblockend+1):length(mismatchY)]==1)
    vectorZ = which(mismatchZ[(lastblockend+1):length(mismatchZ)]==1)
    
    firstY = min(vectorY) + lastblockend
    firstZ = min(vectorZ) + lastblockend
    firstY
    firstZ
    
    cbind(mismatchY, mismatchZ)
    
    if(is.infinite(firstY) && is.infinite(firstZ)==F)
    {
      finalseq[(lastblockend+1):length(mismatchY)] = y
      break
    }
    if(is.infinite(firstZ) && is.infinite(firstY)==F)
    {
      finalseq[(lastblockend+1):length(mismatchZ)] = z
      break
    }
    if(is.infinite(firstY) && is.infinite(firstZ))
    {
      finalseq[(lastblockend+1):length(mismatchZ)] = lastparent
      break
    }
    
    if(firstY==firstZ)
    {
      if(lastparent>0)
      {
        finalseq[(lastblockend+1):firstY] = lastparent
      }
      else
      {
        finalseq[(lastblockend+1):firstY] = y
      }
      lastblockend = firstY
    } else
    {
      if(firstZ > firstY)
      {
        finalseq[(lastblockend+1):(firstZ-1)] = z
        lastblockend = firstZ - 1
        lastparent = z
      } else
      {
        finalseq[(lastblockend+1):(firstY-1)] = y
        lastblockend = firstY - 1
        lastparent = y
      }
    }
  }
  return(finalseq)
}

#--------------------------------------------------------------------
# This function is to correct genotyping errors & switching errors
# It takes two steps
correcting_haplotype = function(p_haplotype,col1,col2,onehap,columnflag,minfraglen,mindiffs)
{
  alleles = unique(onehap)
  #step1
  #remove genotyping errors
  if(length(onehap[onehap==alleles[1]])<=3)
  {
    onehap[onehap==alleles[1]] = alleles[2]
    return(onehap)
  }
  if(length(onehap[onehap==alleles[2]])<=3)
  {
    onehap[onehap==alleles[2]] = alleles[1]
    return(onehap)
  }
  
  #step2
  #remove switching errors with short genomic spans
  finderr=1
  i=1:(length(onehap)-1)
  j=i+1
  while(finderr==1)
  {
    finderr=0
    seps = c(0,which(onehap[i]!=onehap[j]),length(onehap))
    if(length(seps)>=3)
    {
      frags = cbind(seps[1:(length(seps)-1)]+1,seps[2:length(seps)])
      fraglen = frags[,2]-frags[,1]+1
      shortspanid = which(fraglen<=minfraglen)
      if(length(shortspanid)>0)
      {
        for(k in 1:length(shortspanid))
        {
          #if(sum(abs(p_haplotype[frags[k,1]:frags[k,2],col1] - p_haplotype[frags[k,1]:frags[k,2],col2]))<mindiffs)
          #{
            finderr=1
            onehap[frags[shortspanid[k],1]:frags[shortspanid[k],2]] = columnflag - onehap[frags[shortspanid[k],1]:frags[shortspanid[k],2]]
          #}
        }
      }
    }
  }

  #remove switching-cross with a few informative markers
  seps = c(0,which(onehap[i]!=onehap[j]),length(onehap))
  if(length(seps)>=3)
  {
    frags = cbind(seps[1:(length(seps)-1)]+1,seps[2:length(seps)])
    fraglen = frags[,2]-frags[,1]+1
    hapdiff = pedigree_haplotype[[(columnflag-1)/2]]-pedigree_haplotype[[(columnflag+1)/2]]
    for(k in 1:nrow(frags))
    {
      fragdiff = sum(abs(hapdiff[frags[k,1]:frags[k,2]]))
      if(fragdiff<=mindiffs&(frags[k,2]-frags[k,2])<=minfraglen)
      {
        onehap[frags[k,1]:frags[k,2]] = columnflag - onehap[frags[k,1]:frags[k,2]]
      }
    }
  }
  return(onehap)
}

#--------------------------------------------------------
#   Constructing Haplotype
#--------------------------------------------------------
# Recode F2 haplotype using F1 haplotypes
F2RecodedByF1 = haplotype
for(pedigreeid in 1:nrow(f2id))
{
  cat(sprintf("pedigreeid=%d\n",pedigreeid))
  pedigree_haplotype = haplotype[,pedigree_haps[pedigreeid,]]
  if(f2id[pedigreeid,1]%%2==1)
  {
    x=4
    y=11
    z=12
    chrhap1 = StepHap(x,y,z,nrow(haplotype))
    F2RecodedByF1[,pedigree_haps[pedigreeid,][3]] = correcting_haplotype(pedigree_haplotype,11,12,chrhap1,23,20,4)
    F2RecodedByF1[,pedigree_haps[pedigreeid,][4]] = correcting_haplotype(pedigree_haplotype,11,12,chrhap1,23,20,4)
  } else
  {
    score1 = sum(GetScore(3,5)+GetScore(4,11))
    score2 = sum(GetScore(3,5)+GetScore(4,12))
    score3 = sum(GetScore(3,6)+GetScore(4,11))
    score4 = sum(GetScore(3,6)+GetScore(4,12))
    score = c(score1,score2,score3,score4)
    optpos = which(score==min(score))
    optpos
    chrhap1=0
    chrhap2=0
    if(optpos[1]<=4)
    {
      x=3
      y=5
      z=6
      chrhap1 = StepHap(x,y,z,nrow(haplotype))
      x=4
      y=11
      z=12
      chrhap2 = StepHap(x,y,z,nrow(haplotype))
    } else
    {
      x=4
      y=5
      z=6
      chrhap1 = StepHap(x,y,z,nrow(haplotype))
      x=3
      y=11
      z=12
      chrhap2 = StepHap(x,y,z,nrow(haplotype))
    }
    F2RecodedByF1[,pedigree_haps[pedigreeid,][3]] = correcting_haplotype(pedigree_haplotype,5,6,chrhap1,11,20,4)
    F2RecodedByF1[,pedigree_haps[pedigreeid,][4]] = correcting_haplotype(pedigree_haplotype,11,12,chrhap2,23,20,4)  
  }
}

for(pedigreeid in 1:nrow(f2id))
{
  x = cbind(f2id[pedigreeid,1],F2RecodedByF1[,c(1,3)],F2RecodedByF1[,pedigree_haps[pedigreeid,][3]], F2RecodedByF1[,pedigree_haps[pedigreeid,][4]])
  write.table(x, "f2inheritance.chrx.txt",append=T,quote=F,row.names=F,col.names=F,sep="\t")
}
