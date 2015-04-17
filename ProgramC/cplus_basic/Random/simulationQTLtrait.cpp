#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string.h>

using namespace std;

double carmen_gaussian_random(double mean, double std)
{
  //srand((unsigned int)time(NULL)); // use time to generate rand number
  const double norm = 1.0 / (RAND_MAX + 1.0);
  double u = 1.0 - rand() * norm;                  /**//* can't let u == 0 */
  double v = rand() * norm;
  double z = sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
  return mean + std * z;
} 

int main(int argc, char *argv[])
{
   if(argc < 2)
   {
           cout<<"Simulation QTL trait for input samples"<<endl;
           cout<<"Usage:"<<endl;
           cout<<"./simulation sampleRIL.list\n\n"<<endl;
           cout<<"Example file:"<<endl;
           cout<<"Sample1\nSample2\nSample3\n"<<endl;
           exit(1); 
   }
   char outfile[] = ".trait";
   ofstream fout(strcat(argv[1],outfile));
   fout<<"Sample"<<"\t"<<"Trait1\t"<<"Trait2\t"<<"Trait3"<<endl;
   ifstream fin("sampleRIL.list");
   string sample;
   while(getline(fin,sample)){
       double U1  =carmen_gaussian_random (100, 10);
       double U2  =carmen_gaussian_random (30, 10);
       double U3  =carmen_gaussian_random (20, 4);
       fout<<sample<<"\t"<<U1<<"\t"<<U2<<"\t"<<U3<<endl;
   }
   fin.close();
   fout.close();
}
