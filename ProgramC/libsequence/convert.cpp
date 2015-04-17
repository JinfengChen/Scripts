#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <Sequence/Fasta.hpp>
#include <Sequence/Clustalw.hpp>
#include <Sequence/SeqExceptions.hpp>

using namespace std;
using std::ostream_iterator;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::copy;
using std::exception;
using Sequence::Fasta;
using Sequence::ClustalW;
using Sequence::SeqException;

int main (int argc, char **argv)
{
    char *infile = argv[1];
    char *outfile= argv[2];
    
    ClustalW<Fasta> alignment;
    ifstream in(infile);
    try
    {
       in >> alignment;
    }
    catch (SeqException &e)
    {
       cerr << e << endl;
       exit(10);
    }
    catch (exception &e)
    {
       cerr << e.what() << endl;
       exit(10);
    }
    
    for (unsigned seq =0 ; seq < alignment.size(); ++seq) 
    { 
        cout << alignment[seq] << '\n';
    }

    ofstream out(outfile);
    copy(alignment.begin(),alignment.end(), ostream_iterator<Fasta>(out,"\n"));
    return 0;
}

