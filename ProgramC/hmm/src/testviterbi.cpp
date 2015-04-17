#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

#include "hmm.h"

int main(int argc, char* argv[])
{
  Hmm hmm;

  hmm.loadProbs(argv[1]);
  vector<vector<unsigned long>*> seqs;
  hmm.readSeqs(cin, seqs);

  for (unsigned int i = 0; i<seqs.size(); i++) {
    vector<unsigned long>& seq = *seqs[i];
    for (unsigned int j =0; j<seq.size(); j++) {
      hmm.addObservation(seq[j]);
    }
    vector<Transition*> path;
    double jointProb = hmm.viterbi(path);
    cout << "P(path)=" << exp(jointProb-hmm.obsProb()) << endl
	 << "path: " << endl;
    for (unsigned int i = 0; i<path.size(); i++) {
      Transition* trans = path[i];
      if (trans==0) continue;
      cout << hmm.getStr(trans->_obs) << '\t' 
	   << hmm.getStr(trans->_to->state()) << endl;
    }
    hmm.reset();
  }    
}
  
