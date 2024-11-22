#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <sstream>

using namespace std;

using std::ios;
using std::cout;
using std::endl;

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

vector<int> GetSampleGenotypeCode(vector<string> genotype,vector<int> samplepos)
{
    int i;
    vector<int> retval;

    for(i=0;i<samplepos.size();i++)
    {
        if(genotype[samplepos[i]].c_str()[0]=='0'&&genotype[samplepos[i]].c_str()[2]=='0')
        {
          retval.push_back(0);
          continue;
        }
        if(genotype[samplepos[i]].c_str()[0]=='1'&&genotype[samplepos[i]].c_str()[2]=='1')
        {
          retval.push_back(2);
          continue;
        }
        if((genotype[samplepos[i]].c_str()[0]=='0'&&genotype[samplepos[i]].c_str()[2]=='1')||(genotype[samplepos[i]].c_str()[0]=='1'&&genotype[samplepos[i]].c_str()[2]=='0'))
        {
          retval.push_back(1);
          continue;
        }
        retval.push_back(-1);
    }
    return(retval);
}

vector<int> CountSampleGenotype(vector<int> genotype)
{
    int i;
    /*
    vector<int> retval;
    for(i = 0; i<3; i++)
		retval.push_back(0);
    */
    vector<int> retval(3);
    
    for(i=0;i<genotype.size();i++)
    {
	   if(genotype[i]>=0)
		retval[genotype[i]]++;
    }
    return(retval);
}


int main(int argc,char *argv[])
{
  if(argc!=4)
  {
    cout << "Usage: "<<argv[0]<<" windowsize popfile1 popfile2\n";
    return 0;
  }
  string popfile1, popfile2;
  vector<string> pop1_individuals,pop2_individuals;
  vector<int> pop1_individuals_columns,pop2_individuals_columns;
  int pop1size, pop2size;
  vector<int> pop1_genotype,pop2_genotype;

  string lastchr, thischr;
  long windowsize;
  long lastwindow, thiswindow;
  long comparisons;
  long diffs,sitediff,snps,hethomodiffs,hethetdiffs,homohomodiffs;

  unordered_map<string,int> diffs_map,hethomodiffs_map,hethetdiffs_map,homohomodiffs_map,comparisons_map;
  stringstream hash_key;
  string hash_str;

  int snpflag;
  int i, j;

  string linedata;
  vector<string> data_columns;
  vector<int> genotype1,genotype2;
  vector<int> genotypecount1,genotypecount2;
  int genotypecount[3];
  
  windowsize = atol(argv[1]);
  popfile1 = argv[2];
  popfile2 = argv[3];

  pop1_individuals = LoadFileLinesIntoVector(popfile1);
  pop2_individuals = LoadFileLinesIntoVector(popfile2);

  while( getline(cin,linedata) )
  {
      if(linedata.find("#CHROM",0)==0)
          break;
  }
  pop1_individuals_columns = FindSamplesInVCFHeader(linedata, pop1_individuals);
  pop2_individuals_columns = FindSamplesInVCFHeader(linedata, pop2_individuals);

  pop1size = pop1_individuals_columns.size();
  pop2size = pop2_individuals_columns.size();
  
  lastwindow=0;
  thiswindow=0;
  lastchr="";

  cout << "Chr" << "\t" << "Position" << "\tPop1size\tPop2size\t" << "SNPs" << "\t" << "Comparisons" << "\t" << "Difference" << "\t" << "HetHet_Difference" << "\t" << "Het_Homo_Difference" << "\t" << "Homo_Homo_Difference" << "\t" << "Dxy" << "\t" << "SNP-wise dxy"  <<  endl;
  
  while( getline(cin,linedata) )
  {
    //snpflag = 0;
    //data_columns = split(linedata,"\t",true);
    data_columns.clear();
    SplitString(linedata,"\t",data_columns,false);
    thischr = data_columns[0];
    thiswindow = atol(data_columns[1].c_str())/windowsize;
    if(lastchr=="")
    {
      lastchr = thischr;
      lastwindow = thiswindow;
      comparisons = 0;
      hethetdiffs = 0;
      hethomodiffs = 0;
      homohomodiffs = 0;
      diffs = 0;
      snps = 0;
     }
     if( (thischr!=lastchr) || (lastwindow!=thiswindow) )
     {
         if(comparisons>0&&snps>0)
            cout << lastchr << "\t" << lastwindow*windowsize << "\t" << pop1size << "\t" << pop2size << "\t" << snps << "\t" << comparisons << "\t" << diffs << "\t" << hethetdiffs << "\t" << hethomodiffs << "\t" << homohomodiffs << "\t" << double(diffs)/double(windowsize)/(double(comparisons)/double(snps)) << "\t"  << double(diffs)/double(comparisons) << endl;
         else
            cout << lastchr << "\t" << lastwindow*windowsize << "\t" << pop1size << "\t" << pop2size << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;

         lastchr = thischr;
         lastwindow = thiswindow;
         comparisons = 0;
		 hethetdiffs = 0;
		 hethomodiffs = 0;
		 homohomodiffs = 0;
         diffs = 0;
         snps = 0;
     }
      
     genotype1 = GetSampleGenotypeCode(data_columns,pop1_individuals_columns);
     genotype2 = GetSampleGenotypeCode(data_columns,pop2_individuals_columns);
      
     genotypecount1 = CountSampleGenotype(genotype1);
     genotypecount2 = CountSampleGenotype(genotype2);
     if( (genotypecount1[0]==0 && genotypecount1[1]==0 && genotypecount1[2]==0 ) || (genotypecount2[0]==0 && genotypecount2[1]==0 && genotypecount2[2]==0 ))
		continue;
	 genotypecount[0] = genotypecount1[0] + genotypecount2[0];
	 genotypecount[1] = genotypecount1[1] + genotypecount2[1];
	 genotypecount[2] = genotypecount1[2] + genotypecount2[2];
      
     if(genotypecount[1]==0)
     {
		 if(genotypecount[0]==0 || genotypecount[2]==0)
			continue;
     }
     snps ++;

	 hash_key.str("");
     hash_key <<  genotypecount1[0] << "_" << genotypecount1[1] << "_" << genotypecount1[2] << "_" << genotypecount2[0] << "_" << genotypecount2[1] << "_" << genotypecount2[2];
     hash_str = hash_key.str();
     
     if(diffs_map.find(hash_str)!=diffs_map.end())
     {
		diffs += diffs_map[hash_str];
		hethomodiffs += hethomodiffs_map[hash_str];
		hethetdiffs += hethetdiffs_map[hash_str];
		homohomodiffs += homohomodiffs_map[hash_str];
		comparisons += comparisons_map[hash_str];
		continue;
     }

	 diffs_map[hash_str] = 0;
	 hethomodiffs_map[hash_str] = 0;
	 hethetdiffs_map[hash_str] = 0;
	 homohomodiffs_map[hash_str] = 0;
	 comparisons_map[hash_str] = 0;

     comparisons+=genotypecount1[0]*genotypecount2[0]*4;
     comparisons+=genotypecount1[2]*genotypecount2[2]*4;
     comparisons_map[hash_str] += genotypecount1[0]*genotypecount2[0]*4;
     comparisons_map[hash_str] += genotypecount1[2]*genotypecount2[2]*4;
     
     sitediff = genotypecount1[0]*genotypecount2[1]*2;
     hethomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff*2;
	 diffs_map[hash_str] += sitediff;
	 hethomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff*2;
     
     sitediff = genotypecount1[0]*genotypecount2[2]*4;
     homohomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff;
	 diffs_map[hash_str] += sitediff;
	 homohomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff;

     sitediff = genotypecount1[1]*genotypecount2[0]*2;
     hethomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff*2;
	 diffs_map[hash_str] += sitediff;
	 hethomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff*2;

     sitediff = genotypecount1[1]*genotypecount2[1]*2;
     hethetdiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff*2;
	 diffs_map[hash_str] += sitediff;
	 hethetdiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff*2;
	 
     sitediff = genotypecount1[1]*genotypecount2[2]*2;
     hethomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff*2;
	 diffs_map[hash_str] += sitediff;
	 hethomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff*2;

     sitediff = genotypecount1[2]*genotypecount2[0]*4;
     homohomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff;
	 diffs_map[hash_str] += sitediff;
	 homohomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff;

     sitediff = genotypecount1[2]*genotypecount2[1]*2;
     hethomodiffs += sitediff;
     diffs += sitediff;
     comparisons+=sitediff*2;
	 diffs_map[hash_str] += sitediff;
	 hethomodiffs_map[hash_str] += sitediff;
	 comparisons_map[hash_str] += sitediff*2;
     
      /*
      for(i=0; i< pop1size; i++)
      {
          if( genotype1[i]<0)
              continue;
         
          for(j=0; j < pop2size; j++)
          {
             if( genotype2[j]<0)
                  continue;

             snpflag = 1;
             if(genotype1[i]==1 && genotype2[j]==1)
                sitediff = 2;
             else
                sitediff = abs(genotype1[i] - genotype2[j])*2;
             
             comparisons+=4;
             diffs+=sitediff;
          }
      }
      if(snpflag==1)
        snps++;
      */
  }
  if(comparisons>0 && snps>0)
  {
     cout << lastchr << "\t" << lastwindow*windowsize << "\t" << pop1size << "\t" << pop2size << "\t" << snps << "\t" << comparisons << "\t" << diffs << "\t" << hethetdiffs << "\t" << hethomodiffs << "\t" << homohomodiffs << "\t" << double(diffs)/double(windowsize)/(double(comparisons)/double(snps)) << "\t" << double(diffs)/double(comparisons) << endl;
  }
  
  return 1;
}
