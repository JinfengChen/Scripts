#include <vector>
#include <iostream>
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>

//using declarations make ti more explicit where we get our routines from
using std::vector;
using std::cerr;
using std::endl;
using Sequence::Fasta;
using Sequence::Alignment::GetData;
using Sequence::Alignment::IsAlignment;

int main (int argc, char **argv)
{
   vector<Fasta> data;
   char *infile =argv[1];
   GetData(data, infile);

   if (! IsAlignment(data))
   {
      cerr << "file not align!" << endl;
   }
   else
   {
      cerr << "file aligned" << endl;
   }
return 0;
}
