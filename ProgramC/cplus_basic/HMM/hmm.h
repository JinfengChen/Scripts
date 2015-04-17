#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <deque>

using namespace std;

/*
            typedef struct {
               vector<string> M; //number of observation
               vector<string> N; //number of state
               vector< vector<double> > A; //transition probability
               vector< vector<double> > B; //emission probability
               vector<double> pi; // initial probability of state
            } HMMparameter;
*/

class HMM { // define class 
      public:
            typedef struct {
               vector<string> M; //number of observation
               vector<string> N; //number of state
               vector< vector<double> > A; //transition probability
               vector< vector<double> > B; //emission probability
               vector<double> pi; // initial probability of state
               map<string,int> str2id;// store the observation's index a->1,b->2, z->n
               map<int,string> id2str;//store the hiden's index 1->A,2->B,n->Z
            } HMMparameter;

            //typedef vector< vector<string> > OB; //observation sequence, two dimesion matrix of multi sequence

            //OB sequence;
            HMMparameter hmm; // hmm define here, so it can be acessed in all member functions of class HMM
            void readHMM(string name); /*read hmm file, which include observation/hidden state 
                             and transition/emmision/initial probability*/
            double viterbi(vector<string> obseq, vector<string> &state); /* viterbi algorithm to predict
                                           the most possible hidden state of observed sequence*/
            void readsequence(string name, vector< vector<string> > &sequence); // read sequence
            void printdata();
            void printdata(vector<string> seq);
            HMM (); // constructor with parameter
            ~HMM () {}; // destructor, no need to implemented
      private:
};

