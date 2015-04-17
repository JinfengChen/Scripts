#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>


using namespace std;

#include "hmm.h"

int main(int argc, char* argv[])
{
  Hmm hmm;

  if (argc<2) {
    cerr << "USAGE: genseq NAME N" << endl
	 << "generates N observation sequences using the HMM with the given NAME" << endl;
    exit(1);
  }
  hmm.loadProbs(argv[1]);
  int seqs = atoi(argv[2]);
  hmm.genSeqs(cout, seqs);
}
  
