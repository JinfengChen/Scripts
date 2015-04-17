#include <Sequence/Translate.hpp>
#include <iostream>
#include <string>

const char alphabet[4] = {'A','G','C','T'};

int main(int argc, char **argv)
{
  std::string codon;
  codon.resize(3);
  for (unsigned first = 0 ; first < 4 ; ++first)
    {
        for (unsigned second = 0 ; second < 4 ; ++second)
	  {
	    for (unsigned third = 0 ; third < 4 ; ++third)
	      {
		codon[0] = alphabet[first];
		codon[1] = alphabet[second];
		codon[2] = alphabet[third];
		std::cout << codon
			  << '\t' 
			  << Sequence::Translate(codon.begin(),codon.end())
			  << std::endl;
	      }
	  }
    }
}
