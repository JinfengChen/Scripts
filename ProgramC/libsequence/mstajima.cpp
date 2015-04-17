#include <Sequence/SimData.hpp>
#include <Sequence/SimParams.hpp>
#include <Sequence/PolySIM.hpp>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
//to compile:
//g++ -o msfreq msfreq.cpp -lsequence -03
//ms 10 10 -t 10 -r 10 10 | ./msfreq

using namespace std;
using namespace Sequence;

int main (int argc, char *argv[])
{
   SimParams p;
   p.fromfile(stdin);
   
   SimData data;
   std::ios_base::sync_with_stdio(true);
   int rv;
   
   while((rv = data.fromfile(stdin)) != EOF)
   {
      cout << "sample" << endl;
      PolySIM calc(&data);
      cout << calc.TajimasD() << endl;
   }
   return 0;
}

