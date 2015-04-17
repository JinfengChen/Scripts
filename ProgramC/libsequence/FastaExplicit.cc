#include <Sequence/Fasta.hpp>
#include <Sequence/Clustalw.hpp>
#include <Sequence/Alignment.hpp>

/*! \include FastaExplicit.cc
  an example of explicit template instantiation
*/
using namespace std;
namespace Sequence
{
  namespace Alignment
  {
     template void GetData( std::vector<Sequence::Fasta> &, const char *);
     template istream & GetData( std::vector<Sequence::Fasta> &, istream &);
     template istream & ReadNObjects( std::vector<Sequence::Fasta> &, unsigned, istream &);
     template bool Gapped(const std::vector<Sequence::Fasta> &);
     template bool IsAlignment(const std::vector<Sequence::Fasta> &);
     template unsigned UnGappedLength(const std::vector<Sequence::Fasta> &);
     template void RemoveGaps(std::vector<Sequence::Fasta> &);
     template void RemoveTerminalGaps(std::vector<Sequence::Fasta> &);
     template std::vector<Sequence::Fasta> 
     Trim(const std::vector<Sequence::Fasta> &, const std::vector<int> &);
     template std::vector<Sequence::Fasta> 
     TrimComplement(const std::vector<Sequence::Fasta> &, const std::vector<int> & );
  }
  template class ClustalW<Fasta>;
}

using namespace Sequence;
using namespace Alignment;

int main(int argc, char **argv)
{
  return 0;
}
