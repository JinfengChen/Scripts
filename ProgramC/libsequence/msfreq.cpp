#include <Sequence/SimData.hpp>
#include <Sequence/SimParams.hpp>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
//to compile:
//g++ -o msfreq msfreq.cpp -lsequence -03
//ms 10 10 -t 10 -r 10 10 | ./msfreq

using namespace std;
using namespace Sequence;

int main (int argc, char *argv[])
{
   SimParams p;
   p.fromfile(stdin);
   
   SimData data;
   std::ios_base::sync_with_stdio(true);
   int rv;
   
   while((rv = data.fromfile(stdin)) != EOF)
   {
      cout << "sample" << endl;
      vector<unsigned> freqSpec(data.size()-1,0u);
      for (unsigned site =0; site < data.numsites(); ++site)
      {
         unsigned nc =0;
         for (unsigned seq =0; seq < data.size(); ++seq)
         {
             nc += (data[seq][site] == '1') ? 1 : 0; //for every SNP site, calculate number of derived mutation
         }
         freqSpec[nc-1]++; //count for the frequency of allele frequency, how many times this allele frequency is present in this region
      }
      for (unsigned freq = 0;freq < freqSpec.size();++freq)
      {
         fprintf(stdout,"%d %d\n",freq+1,freqSpec[freq]);
      }
   }
   return 0;
}

