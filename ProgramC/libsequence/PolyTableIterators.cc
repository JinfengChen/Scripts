/*! \include PolyTableIterators.cc */

/*
  Example of how to use the iterators
  of class Sequence::PolyTable. In
  this example, Sequence::PolySites is used,
  which is derived from Sequence::PolyTable
*/
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;


int main(int argc, char **argv)
{
  const char *infilename = argv[1];

  vector<Sequence::Fasta> data;
  
  Sequence::Alignment::GetData(data,infilename);

  if ( Sequence::Alignment::IsAlignment(data) &&
       Sequence::Alignment::validForPolyAnalysis(data.begin(),data.end()) )
    {
      Sequence::PolySites SNPs(data);

      //1. use PolyTable::pos_iterator to access site positions
      
      //Print positions of segregating sites
      cout << "Original position order:\n";
      copy(SNPs.pbegin(),SNPs.pend(),
	   ostream_iterator<double>(cout," "));
      cout << endl;

      //2. PolyTable::pos_iterator can be accessed in const-
      // and non-const contexts, allowing us to do things like
      // permute site positions.  Note that this only permutes
      // the positions, not the states associated with them.
      // This allows, amongst other things, to calculate the
      // significance of linkage-disequilibrium measures by
      // a permutation test.
      random_shuffle(SNPs.pbegin(),SNPs.pend());

      cout << "Permuted positions:\n";
      copy(SNPs.pbegin(),SNPs.pend(),
	   ostream_iterator<double>(cout," "));
      cout <<'\n'<< endl;

      //3. Access the individuals using iterators.
      //This iterator type is PolyTable::data_iterator,
      //which can be accessed in both const- and non-const
      //contexts.
      cout << "Original data:\n";
      copy(SNPs.begin(),SNPs.end(),ostream_iterator<string>(cout,"\n"));
      cout <<'\n'<< endl;

      //4. We can use PolyTable::data_iterator in
      //non-const contexts to do things like permute
      //the order of the haplotypes.  One application
      //of this would be to assess the significance of 
      //population structure statistics (Fst and the like)
      //by permutation tests
      random_shuffle(SNPs.begin(),SNPs.end());
      cout << "Permuted data:\n";
      copy(SNPs.begin(),SNPs.end(),ostream_iterator<string>(cout,"\n"));
      cout << endl;

      //5. There is a special iterator, PolyTable::const_site_iterator,
      //which gives access to a single SNP in a const-context.  So far,
      //only const access is supported.  The value type of this iterator
      //is equivalent to std::pair<double,std::string>.  We will use this
      //iterator to calculate nucleotide diversity for the sample.  The 
      //defnitiion of nucleotide diversity does not depend on either the order
      //of the SNPs nor the arrangement of haplotypes, so we're able to do the 
      //calculation on the data that we've permuted.
      //We handle missing data by adjusting the sample size for each site.
      //This implementation is contrived, and we'd really use
      //Sequence::stateCounter or Sequence::makeCountList to handle the 
      //counting of the nucleotides for us.
      double pi = 0.;
      Sequence::PolySites::const_site_iterator itr = SNPs.sbegin();
      while( itr < SNPs.send() )
	{
	  //count up numbers of A,G,C,T, and N using std::count
	  unsigned A = count(itr->second.begin(),itr->second.end(),'A');
	  unsigned G = count(itr->second.begin(),itr->second.end(),'G');
	  unsigned C = count(itr->second.begin(),itr->second.end(),'C');
	  unsigned T = count(itr->second.begin(),itr->second.end(),'T');
	  unsigned N = count(itr->second.begin(),itr->second.end(),'N');
	  double SH = 1.; //SH = site heterozygosity

	  //sample size at this site, w/o missing data
	  unsigned n = itr->second.length() - N; 

	  //For each character state, the homozygosity is the probability
	  //of sampling the observed count for that state (k) out of n alleles times the 
	  //probability of sampling it again out of n-1 alleles, given that you've sampled
	  //it once (i.e. (k-1)/(n-1).
	  //Pi is the sum of 1-homozygosity accross all segregating sites
	  SH -= (A>0) ? (double(A)/double(n))*(double(A-1)/double(n-1)) : 0.;
	  SH -= (G>0) ? (double(G)/double(n))*(double(G-1)/double(n-1)) : 0.;
	  SH -= (C>0) ? (double(C)/double(n))*(double(C-1)/double(n-1)) : 0.;
	  SH -= (T>0) ? (double(T)/double(n))*(double(T-1)/double(n-1)) : 0.;
	  pi += SH;
	  ++itr;
	}
      cout << "Pi (for the entire region) = "<< pi << endl;
    }
}
