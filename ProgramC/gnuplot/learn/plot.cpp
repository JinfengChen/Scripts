#include<iostream>
#include<fstream>

using namespace std;

int main ()
{
    ofstream outfile("Preplot.plt");
    if (!outfile) 
    {
       cout<<"outfile failed"<<endl;
    }
    else
    {
       outfile<<"set xlabel \"Month\""<<'\n'
              <<"set ylabel \"Mass (mm)"<<'\n'
              <<"set term pdfcairo lw 2 font \"Times_New_Roman,8\""<<'\n'
              <<"set output \"precipitation.pdf\""<<'\n'
              <<"plot \"precipitation.dat\" u 1:2 w lp pt 5 title \"Beijing\","
              <<"\"precipitation.dat\" u 1:3 w lp pt 7 title \"Shanghai\""<<endl;
       outfile.close();
    }
}
