#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>

#if defined (__GNUG__) && __GNUC__ > 3
#include <Sequence/FastaExplicit.hpp>
#else
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#endif
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/SeqExceptions.hpp>

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;

struct params
{
  char * infile;
  bool purgeMissing;
  bool purgeMultiHits;
  bool haveOutgroup;
  unsigned outgroup;
};

int parseargs (int argc, char **argv, params *args)
{
  args->infile = NULL;
  args->purgeMissing = false;
  args->purgeMultiHits = false;
  args->haveOutgroup = false;
  args->outgroup = 0;

  extern int optind;
  int c;
  while ((c = getopt (argc, argv, "i:o:nb")) != -1)
  {
      cout << c << endl;
      switch (c)
      {
         case 'i':
             args->infile = optarg;
             break;
         case 'o':
             args->haveOutgroup = true;
             args->outgroup = atoi(optarg);
             break;
         case 'n':
             args->purgeMissing = true;
             break;
         case 'b':
             args->purgeMultiHits = true;
             break;
         default:
             cerr << "Huh?" << endl;
             exit(1);
             break;
      }
  }
  return optind;
}

int main (int argc, char **argv)
{
   params args;
   parseargs (argc,argv,&args); 
   if (args.infile == NULL)
   { 
      cerr << "error: no infile specificed!\n";
      exit(1);
   }
   vector<Fasta> alignment;
   try
   {
      GetData(alignment, args.infile);
   }
   catch(SeqException &e)
   {
      cerr << e << endl;
      exit(1);
   }
   if (args.outgroup >= alignment.size())
   {
      cerr << "error: outgroup sequence out of range!\n";
      exit(1);
   }
   if (! IsAlignment(alignment))
   {
      cerr << "error: data not aligned!\n";
      exit(1);
   }
   PolySites SNPtable(alignment, args.purgeMultiHits, true, args.purgeMissing);
   PolySNP calculator(&SNPtable, args.haveOutgroup, args.outgroup);
   cout << "Tajima's D for file " << args.infile
        << " is: " << calculator.TajimasD() << endl;
   return 0;
}





