/*! \include valid_dna.cc */
#include <Sequence/Fasta.hpp>
#include <Sequence/SeqRegexes.hpp>
#include <Sequence/SeqProperties.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>

int main(int argc, char **argv)
{
  std::ifstream in(argv[1]);
  Sequence::Fasta seq;
  while (!in.eof())
    {
      in >> seq;
      std::cout << Sequence::validSeq(seq.begin(),seq.end())
		<< '\t'
		<< Sequence::validSeq(seq.begin(),seq.end(),Sequence::full_dna_alphabet)
		<< '\t'
		<< (std::find_if(seq.begin(),seq.end(),Sequence::ambiguousNucleotide())
		    != seq.end())
                << '\t'
                << seq.GetName()
		<<'\n';
    }
}
