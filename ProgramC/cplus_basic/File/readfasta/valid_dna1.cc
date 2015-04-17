/*! \include valid_dna.cc */
#include <Sequence/Seq.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/SeqRegexes.hpp>
#include <Sequence/SeqProperties.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

int main(int argc, char **argv)
{
  if (argc == 1)
  {
     cout<< "Usage: ./valid pop.fa"<<endl;
     exit(1);
  }
  std::ifstream in(argv[1]);
  Sequence::Fasta seq;
  while (!in.eof())
    {
      in >> seq;
      Sequence::Fasta temp;
      temp=seq;
      temp.Revcom();
      std::cout << Sequence::validSeq(seq.begin(),seq.end())
		<< '\t'
		<< Sequence::validSeq(seq.begin(),seq.end(),Sequence::full_dna_alphabet)
		<< '\t'
		<< (std::find_if(seq.begin(),seq.end(),Sequence::ambiguousNucleotide())
		    != seq.end())
                << '\t'
                << seq.GetName()<<'\t'<<seq.length()
                << '\n'
                << "Sequence"<<'\n'<<seq.GetSeq()
                << '\n'
                << "Revcom seq"<<'\n'<<temp.GetSeq()<<'\n'
                << "Subseq,start=1,len=100"<<'\n'<<seq.substr(1,100)<<'\n'
		<<'\n';
    }
}
