#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <sstream>

using namespace std;

using std::ios;
using std::cout;
using std::endl;

long sampleocc1[10][1024][3];
long sampleocc2[10][1024][3];
vector<string> headervec;
int GenotypeInfoSize=3;

vector<string> split(const string& s, const string& delim, const bool keep_empty = true)
{
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin();
    string::const_iterator subend;
    while (true)
    {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if(keep_empty || !temp.empty())
        {
            result.push_back(temp);
        }
        if (subend == s.end())
        {
            break;
        }
        substart = subend + delim.size();
    }
    return result;
}

void SplitString(const string& str, const string& delimiters, vector<string> &elems, bool skip_empty=false)
{
    string::size_type pos, prev = 0;
    while((pos=str.find_first_of(delimiters,prev))!=string::npos)
    {
        if(pos>prev)
        {
            if(skip_empty && 1== pos - prev)
                break;
            elems.emplace_back(str,prev,pos-prev);
        }
        prev = pos + 1;
    }
    if(prev < str.size())
        elems.emplace_back(str, prev, str.size() - prev);
}

vector<string> LoadFileLinesIntoVector(string filename)
{
    ifstream infile;
    string data;
    vector<string> retdata;

    infile.open(filename.c_str());
    while (getline(infile,data))
    {
        if(data.length()==0)
        {
            continue;
        }
        retdata.push_back(data);
    }
    infile.close();

    return retdata;
}

vector<int> FindSamplesInVCFHeader(string header, vector<string> sampleid)
{
    int i,j;
    vector<int> retpos;
    vector<string> cols;

    SplitString(header,"\t",cols,false);
    //cols = split(header,"\t",true);
    for(i=0;i<sampleid.size();i++)
    {
        for(j=0;j<cols.size();j++)
        {
            if(cols[j] == sampleid[i])
            {
                retpos.push_back(j);
                break;
            }
        }
        if(j==cols.size())
            retpos.push_back(-1);
    }
    return(retpos);
}

vector<int> GetSampleGenotypeCodeWithDepth(vector<string> genotype,vector<int> samplepos,int mindepth,int maxdepth)
{
    int i;
    vector<int> retval;
    vector<string> fieldinfo;
    char base1,base2;
    int depth;
    int samplesize;
    samplesize = samplepos.size();

    for(i=0;i<samplesize;i++)
    {
        fieldinfo.clear();
        SplitString(genotype[samplepos[i]],":",fieldinfo,false);
        depth = atoi(fieldinfo[2].c_str());
        if( depth<mindepth || depth>maxdepth )
            continue;
        base1 = genotype[samplepos[i]].c_str()[0];
        base2 = genotype[samplepos[i]].c_str()[2];
        if( base1 =='0'&& base2 =='0')
        {
          retval.push_back(0);
          retval.push_back(samplepos[i]);
          retval.push_back(depth);
          continue;
        }
        if( base1=='1' && base2=='1')
        {
          retval.push_back(2);
          retval.push_back(samplepos[i]);
          retval.push_back(depth);
          continue;
        }
        if((base1=='0'&& base2=='1')||(base1=='1'&& base2=='0'))
        {
          retval.push_back(1);
          retval.push_back(samplepos[i]);
          retval.push_back(depth);
          continue;
        }
        //retval.push_back(-1);
        //retval.push_back(samplepos[i]);
    }
    return(retval);
}

vector<int> CountSampleGenotype(vector<int> genotype)
{
    int i;
    vector<int> retval(3);
    int samplesize;
    samplesize = genotype.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype[i]>=0)
        retval[genotype[i]]++;
    }
    return(retval);
}

void CountSampleGenotypeStat(vector<int> genotype, int popid, int DAFid)
{
    int i;
    int samplesize;
    samplesize = genotype.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype[i]>=0)
      {
        if(popid == 1)
            sampleocc1[DAFid][genotype[i+1]][genotype[i]]++;
        else
            sampleocc2[DAFid][genotype[i+1]][genotype[i]]++;
      }
    }
}

vector<float> GetSampleGenotypeMeanDepth(vector<int> genotype)
{
    int i;
    int samplesize;
    vector<float> sumdepth(3);
    samplesize = genotype.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype[i]>=0)
      {
        sumdepth[genotype[i]] = sumdepth[genotype[i]] + genotype[i+2];
      }
    }
    return(sumdepth);
}

float CountMinorAlleleFreq(vector<int> genotype1, vector<int> genotype2)
{
    float allsize;
    float allelecount;
    float allelefreq;
    allsize = 0.0;
    allelecount = 0.0;
    int i;
    int samplesize;
    
    samplesize = genotype1.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype1[i]>=0)
      {
        allsize = allsize+2;
        allelecount += genotype1[i];
      }
    }
    samplesize = genotype2.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype2[i]>=0)
      {
        allsize = allsize+2;
        allelecount += genotype2[i];
      }
    }
    allelefreq = allelecount/allsize;
    if(allelefreq>0.5)
        allelefreq = 1.0 - allelefreq;
    return(allelefreq);
}

float CountSexChrMinorAlleleFreq(vector<int> genotype1, vector<int> genotype2)
{
    float allsize;
    float allelecount;
    float allelefreq;
    allsize = 0.0;
    allelecount = 0.0;
    int i;
    int samplesize;
    
    samplesize = genotype1.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype1[i]==0||genotype1[i]==2)
      {
        allsize = allsize+1;
        allelecount += genotype1[i]/2;
      }
    }
    samplesize = genotype2.size();
    for(i=0;i<samplesize;i=i+GenotypeInfoSize)
    {
      if(genotype2[i]>=0)
      {
        allsize = allsize+2;
        allelecount += genotype2[i];
      }
    }
    allelefreq = allelecount/allsize;
    if(allelefreq>0.5)
        allelefreq = 1.0 - allelefreq;
    return(allelefreq);
}


vector<float> calcHWEdev(vector<int> genotypecount1,vector<int> genotypecount2)
{
    int samplecount,pcount,qcount;
    float expecthomo1,expecthomo2,expecthet;
    float p,q;
    float homo1diff,homo2diff,hetdiff;
    vector<float> hwe_expect;
    
    samplecount = genotypecount1[0] + genotypecount1[1] + genotypecount1[2];
    pcount = genotypecount1[0]*2 + genotypecount1[1];
    qcount = samplecount*2 - pcount;
    p = ((float)pcount)/(samplecount*2);
    q = 1 - p;
    expecthomo1 = samplecount * p * p;
    expecthomo2 = samplecount * q * q;
    expecthet = samplecount * p * q * 2;
    homo1diff = genotypecount1[0] - expecthomo1;
    homo2diff = genotypecount1[2] - expecthomo2;
    hetdiff = genotypecount1[1] - expecthet;
    
    hwe_expect.push_back(expecthomo1);
    hwe_expect.push_back(expecthet);
    hwe_expect.push_back(expecthomo2);
    hwe_expect.push_back(homo1diff);
    hwe_expect.push_back(hetdiff);
    hwe_expect.push_back(homo2diff);

    samplecount = genotypecount2[0] + genotypecount2[1] + genotypecount2[2];
    pcount = genotypecount2[0]*2 + genotypecount2[1];
    qcount = samplecount*2 - pcount;
    p = ((float)pcount)/(samplecount*2);
    q = 1 - p;
    expecthomo1 = samplecount * p * p;
    expecthomo2 = samplecount * q * q;
    expecthet = samplecount * p * q * 2;
    homo1diff = genotypecount2[0] - expecthomo1;
    homo2diff = genotypecount2[2] - expecthomo2;
    hetdiff = genotypecount2[1] - expecthet;
    
    hwe_expect.push_back(expecthomo1);
    hwe_expect.push_back(expecthet);
    hwe_expect.push_back(expecthomo2);
    hwe_expect.push_back(homo1diff);
    hwe_expect.push_back(hetdiff);
    hwe_expect.push_back(homo2diff);
    
    return(hwe_expect);
}

int main(int argc,char *argv[])
{
  if(argc!=8)
  {
    cout << "Usage: "<<argv[0]<<" sexchrflag mindepth1 maxdepth1 mindepth2 maxdepth2 popfile1 popfile2\n";
    return 0;
  }
  string popfile1, popfile2;
  vector<string> pop1_individuals,pop2_individuals;
  vector<int> pop1_individuals_columns,pop2_individuals_columns;
  vector<float> hwe_diff;
  int pop1size, pop2size;
  int flag;
  int mindepth1,maxdepth1,mindepth2,maxdepth2;

  string lastchr, thischr;
  long windowsize;
  long lastwindow, thiswindow;
  int homo1,het1,homo2,het2;

  int i, j, k;

  string linedata;
  vector<string> data_columns;
  vector<int> genotype1,genotype2;
  vector<int> genotypecount1,genotypecount2;

  int DAFid;
  float allelefreq;
  vector<float> meandepth1, meandepth2;

  flag = atoi(argv[1]);
  mindepth1 = atoi(argv[2]);
  maxdepth1 = atoi(argv[3]);
  mindepth2 = atoi(argv[4]);
  maxdepth2 = atoi(argv[5]);
  popfile1 = argv[6];
  popfile2 = argv[7];

  pop1_individuals = LoadFileLinesIntoVector(popfile1);
  pop2_individuals = LoadFileLinesIntoVector(popfile2);
  
  while( getline(cin,linedata) )
  {
      if(linedata.find("#CHROM",0)==0)
          break;
  }
  
  SplitString(linedata,"\t",headervec,false);
  pop1_individuals_columns = FindSamplesInVCFHeader(linedata, pop1_individuals);
  pop2_individuals_columns = FindSamplesInVCFHeader(linedata, pop2_individuals);

  pop1size = pop1_individuals_columns.size();
  pop2size = pop2_individuals_columns.size();

  cout << fixed << setprecision(2);
  
  while( getline(cin,linedata) )
  {
     data_columns.clear();
     SplitString(linedata,"\t",data_columns,false);

     genotype1 = GetSampleGenotypeCodeWithDepth(data_columns,pop1_individuals_columns,mindepth1,maxdepth1);
     genotype2 = GetSampleGenotypeCodeWithDepth(data_columns,pop2_individuals_columns,mindepth2,maxdepth2);

     genotypecount1 = CountSampleGenotype(genotype1);
     genotypecount2 = CountSampleGenotype(genotype2);
     if ( genotypecount1[0]==0 && genotypecount1[1]==0 && genotypecount2[0]==0 && genotypecount2[1]==0 )
        continue;
     if ( genotypecount1[2]==0 && genotypecount1[1]==0 && genotypecount2[2]==0 && genotypecount2[1]==0 )
        continue;
     if( (genotypecount1[0]==0 && genotypecount1[1]==0 && genotypecount1[2]==0 ) || (genotypecount2[0]==0 && genotypecount2[1]==0 && genotypecount2[2]==0 ) )
        continue;
 
     meandepth1 = GetSampleGenotypeMeanDepth(genotype1);
     meandepth2 = GetSampleGenotypeMeanDepth(genotype2);

     if(flag==1)
        allelefreq = CountSexChrMinorAlleleFreq(genotype1, genotype2);
     else
        allelefreq = CountMinorAlleleFreq(genotype1, genotype2);
     
     DAFid = static_cast<int>(allelefreq/0.05);
     if(DAFid==10) DAFid = 9;

     cout << data_columns[0] << "\t" << data_columns[1] << "\t" << DAFid 
     << "\t" << genotypecount1[0] << "\t" << genotypecount1[1] << "\t" << genotypecount1[2] << "\t"
     << (genotypecount1[0]>0?static_cast<float>(meandepth1[0])/static_cast<float>(genotypecount1[0]):0) 
     << "\t" << (genotypecount1[1]>0?static_cast<float>(meandepth1[1])/static_cast<float>(genotypecount1[1]):0) 
     << "\t" << (genotypecount1[2]>0?static_cast<float>(meandepth1[2])/static_cast<float>(genotypecount1[2]):0) 
     << "\t" << genotypecount2[0] << "\t" << genotypecount2[1] << "\t" << genotypecount2[2] 
     << "\t" << (genotypecount2[0]>0?static_cast<float>(meandepth2[0])/static_cast<float>(genotypecount2[0]):0) 
     << "\t" << (genotypecount2[1]>0?static_cast<float>(meandepth2[1])/static_cast<float>(genotypecount2[1]):0) 
     << "\t" << (genotypecount2[2]>0?static_cast<float>(meandepth2[2])/static_cast<float>(genotypecount2[2]):0);

     //hwe_diff = calcHWEdev(genotypecount1,genotypecount2);
     //for(j=0;j<12;j++)
     //{
     //   cout << "\t" << hwe_diff[j];
     //}

     cout << endl;
  }
  
  return 1;
}



