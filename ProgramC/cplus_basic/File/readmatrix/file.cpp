#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>

using namespace std;

int main ()
{
    ifstream infile("table.txt");
    ofstream outfile("output.txt");
    if (!infile)
    {
       cout<<"Error opening"<<"table.txt"<<endl;
    }
    else
    {
       string line;
       vector< vector<string> > data;// two dimension vector
       while(getline(infile, line)){
          if (line.size() > 0) // not a blank line
          {
             char *temp=new char[line.size()+1]; //must initialize the char pointer
             strcpy(temp,line.c_str());//copy "const char *" to "char *"
             outfile<<temp<<endl;
             char *word;
             vector<string> tempv;
             while(word=strsep(&temp,"\t")){ // seperate temp with "\t"
                tempv.push_back(word);
                outfile<<word<<endl;
             }
             data.push_back(tempv);
             outfile<<"Read line from file:"<<line<<endl;
          }
       }
       for (int i=0; i<data.size();i++)
       {
           cout<<"col2: "<<data[i][1]<<endl;
           int x=atoi(data[i][1].c_str())+40;//convert string to int before adding
           cout<<"col2+40: "<<x<<endl;
       }
       for (vector< vector<string> >::iterator it=data.begin(); it!=data.end();++it)
       {
           cout<<"iterator read: "<<(*it)[1]<<endl;
       }
    }
    infile.close();
    outfile.close();
}
