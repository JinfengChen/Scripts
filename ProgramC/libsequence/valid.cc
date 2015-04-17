#include <Sequence/Fasta.hpp>
#include <Sequence/SeqUtilities.hpp>
#include <fstream> 
#include <iostream>

int main (int argc, char **argv)
{
    std::ifstream in(argv[1]);
    Sequence::Fasta fseq;
    while (!in.eof()){
       in >> fseq;
       std::cout<< "5"// "Revcom" the sequence
       << "\t"
       << fseq.GetName()
       <<"\n";
    }
}

