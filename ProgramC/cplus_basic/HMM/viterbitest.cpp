#include <iostream>
#include <cstdlib>
#include <vector>
#include <deque>

#include "hmm.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
       cout<<"Usage: viterbitest test.hmm test.seq"<<"\n"<<endl;
       exit(0);
    }
    HMM hmm;
    hmm.readHMM(argv[1]); //read hmm
    //hmm.printdata();
    vector< vector<string> > seqs;
    hmm.readsequence(argv[2], seqs); //read sequence
    for(int i=0;i<seqs.size();i++){
       cout<<"Observation index:"<<i<<endl;
       hmm.printdata(seqs[i]);
       vector<string> state(seqs[i].size()+1);
       double prob=hmm.viterbi(seqs[i], state); //predict hidden state
       cout<<"Hiden sequence: "<<"P="<<prob<<endl;

       for(vector<string>::iterator it= state.begin()+1;it!=state.end();it++)
       {
          cout<<*it;
       }
       cout<<endl;
    }
}
