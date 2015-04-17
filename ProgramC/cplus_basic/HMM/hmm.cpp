#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <deque>

#include "hmm.h"

using namespace std;

//output the hmm parameter to screen
//test version, just output observations
void HMM::printdata()
{
   cout<<"Printing Data:"<<endl;
   for(vector<string>::iterator it=hmm.M.begin();it!=hmm.M.end();it++)
   {
      cout<<"Observation:"<<endl;
      cout<<*it<<endl;
   }

}

//print data in vector, this function is overload of printdata
void HMM::printdata(vector<string> seq)
{
   //cout<<"Printing Data:"<<endl;
   for(vector<string>::iterator it=seq.begin();it!=seq.end();it++)
   {
      cout<<*it;
   }
   cout<<endl;
}



//read sequence of observation and store the sequence in two dimension matrix
//one should use the sequence in "for cycles"
void HMM::readsequence(string name, vector< vector<string> > &sequence)
{
   ifstream sequencefile(name.c_str());
   if (!sequencefile)
   {
      cout<<"Error in read sequence file: "<<name<<endl;
   }
   else
   {
      cout<<"Reading observed sequence......"<<endl;
      string line;
      while(getline(sequencefile,line))
      {
         if (line.size() > 0) // no blanc line
         {
             stringstream ss(line); // make string "line" to sstream
             string temp;
             vector<string> t;
             while( ss>>temp ){ // reading ss into temp one by one
               t.push_back(temp); 
             }
             ss.clear();
             sequence.push_back(t);
         }
      }
      cout<<"Number of sequence:"<<sequence.size()<<endl;
   }
}


//viterbi algorithm
double HMM::viterbi(vector <string> obseq, vector<string> &state) 
{
   cout<<"Doing prediction using viterbi algorithm......"<<endl;
   int i,j; //state indicator
   int t; // time index
   int m=obseq.size(),n=hmm.N.size();
   cout<<"observation sequence M:"<<m<<"\n"<<"state N:"<<n<<endl;
   double delta[m][n]; // pai[i]*B[i][j], probability of state i at t time
   int psi[m][n]; // best t-1 state of t time i state
   int maxvalindex; // index for the best state of t-1 time
   double maxval,val; //
   double prob=0.0;//probability of best hidden sequence
   vector<int> stateindex(m+1); // store the state index for final path
 
   //we start state and sequence from 1, so the last element of vector is index by m or n.
   //initialization,t=1
   for (i=1;i<=n;i++) //for every state in time 1
   {
       delta[1][i]=hmm.pi[i-1]*hmm.B[i-1][hmm.str2id[obseq[0]]];//t=1, 0=1-1
       psi[1][i]=0;
   }

   //recursion
   for (t=2;t<=m;t++) //for sequence from 2 to m
   {
       for(j=1;j<=n;j++) //for each state at t time
       {
          maxval=0.0;
          maxvalindex=1;
          for (i=1;i<=n;i++) // for each state at t-1 time
          {
              val=delta[t-1][i]*hmm.A[i-1][j-1];
              if (val > maxval)
              {
                 maxval=val;
                 maxvalindex=i;
              }
          }
          delta[t][j]=maxval*hmm.B[j-1][hmm.str2id[obseq[t-1]]];
          psi[t][j]=maxvalindex;
       }
   }
   
   //termination
   for (i=1;i<=n;i++) //for each state t=m
   { 
       if (delta[m][i] > prob)//find the best state with largest delta
       {
          prob=delta[m][i];
          stateindex[m]=i;
          state[m]=hmm.id2str[stateindex[m]-1];
       }
   }

   //path backtracking
   for(t=m-1;t>=1;t--)
   {
      stateindex[t]=psi[t+1][stateindex[t+1]];
      state[t]=hmm.id2str[stateindex[t]-1];
   }
   return prob;
}

//modify viterbi to allow variation in emmision probability
double HMM::viterbi(vector <string> obseq, vector<string> &state)
{
   cout<<"Doing prediction using viterbi algorithm......"<<endl;
   int i,j; //state indicator
   int t; // time index
   int m=obseq.size(),n=hmm.N.size();
   cout<<"observation sequence M:"<<m<<"\n"<<"state N:"<<n<<endl;
   double delta[m][n]; // pai[i]*B[i][j], probability of state i at t time
   int psi[m][n]; // best t-1 state of t time i state
   int maxvalindex; // index for the best state of t-1 time
   double maxval,val; //
   double prob=0.0;//probability of best hidden sequence
   vector<int> stateindex(m+1); // store the state index for final path

   //we start state and sequence from 1, so the last element of vector is index by m or n.
   //   //initialization,t=1
   for (i=1;i<=n;i++) //for every state in time 1
   {
       delta[1][i]=hmm.pi[i-1]*hmm.B[i-1][hmm.str2id[obseq[0]]];//t=1, 0=1-1
       psi[1][i]=0;
   }
   //recursion
   for (t=2;t<=m;t++) //for sequence from 2 to m
   {
       for(j=1;j<=n;j++) //for each state at t time
       {
          maxval=0.0;
          maxvalindex=1;
          for (i=1;i<=n;i++) // for each state at t-1 time
          {
              val=delta[t-1][i]*hmm.A[i-1][j-1];
              if (val > maxval)
              {
                 maxval=val;
                 maxvalindex=i;
              }
          }
          delta[t][j]=maxval*hmm.B[j-1][hmm.str2id[obseq[t-1]]];
          psi[t][j]=maxvalindex;
       }
   }
   //termination
   for (i=1;i<=n;i++) //for each state t=m
   {
       if (delta[m][i] > prob)//find the best state with largest delta
       {
          prob=delta[m][i];
          stateindex[m]=i;
          state[m]=hmm.id2str[stateindex[m]-1];
       }
   }
   //path backtracking
   for(t=m-1;t>=1;t--)
   {
      stateindex[t]=psi[t+1][stateindex[t+1]];
      state[t]=hmm.id2str[stateindex[t]-1];
   }
   return prob;
}



//read hmm file and store the data in HMMparameter struct
/*
M: observation
1 2
N: hidden stage
1 2 3
A: transition matrix
0.333 0.333 0.333
0.333 0.333 0.333
0.333 0.333 0.333
B: emmision matrix
0.5   0.5  
0.75  0.25
0.25  0.75
pi: initial probability
0.333 0.333 0.333
*/
void HMM::readHMM(string name)
{
   cout<<"Reading HMM model......"<<endl;
   ifstream hmmfile(name.c_str()); 
   if (!hmmfile)
   {
      cout<<"Error in read hmm file: "<<name<<endl;
   }
   else
   {
      string line;
      int i=0,j=0,m=0,n=0;
      while(getline(hmmfile,line))
      {
         if (line.size() > 0) // no blanc line
         {
           if (line.find("M:") != line.npos) // find observation line
           {
             getline(hmmfile,line); // read observation in next line
             stringstream ss(line); // make string "line" to sstream
             string temp;
             while( ss >> temp ){ // reading ss into temp one by one
               hmm.M.push_back(temp);
               //cout<<hmm.M[i]<<endl;
               hmm.str2id[temp]=i;
               i++;
             }
             ss.clear();
           }

           if (line.find("N:") != line.npos)
           {
             getline(hmmfile,line); // read state in next line
             stringstream ss(line); // make string "line" to sstream
             string temp;
             while( ss >> temp ){ // reading ss into temp one by one
               hmm.N.push_back(temp);
               hmm.id2str[j]=temp;
               //cout<<hmm.N[j]<<endl;
               j++;
             }
             ss.clear();
           }
           if (line.find("pi:") != line.npos)
           {
             getline(hmmfile,line); // read initial probability in next line
             stringstream ss(line); // make string "line" to sstream
             string temp;
             while( ss >> temp ){ // reading ss into temp one by one
               hmm.pi.push_back(atof(temp.c_str()));
               //cout<<hmm.pi[m]<<endl;
               m++;
             }
             ss.clear();
           }

           if (line.find("A:") != line.npos)   
           {
             for (int x=0;x<j;x++)
             {
                getline(hmmfile,line);
                stringstream ss(line);
                string temp;
                vector<double> a1;
                for (int y=0;y<j;y++)
                {
                   ss >> temp;
                   //cout<<temp<<endl;
                   a1.push_back(atof(temp.c_str()));
                }
                hmm.A.push_back(a1);
             }
           }        

           if (line.find("B:") != line.npos)
           {
             for (int x=0;x<j;x++)
             {
                getline(hmmfile,line);
                stringstream ss(line);
                string temp;
                vector<double> b1;
                for (int y=0;y<i;y++)
                {
                   ss >> temp;
                   //cout<<temp<<endl;
                   b1.push_back(atof(temp.c_str()));
                }
                hmm.B.push_back(b1);
             }  
           }

         }
      }
   }
}

HMM::HMM()
{}

