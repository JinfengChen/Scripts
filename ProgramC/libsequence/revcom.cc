#include <Sequence/Fasta.hpp> //declaration of Sequence::Fasta
#include <iostream>
using namespace std;

int main (int argc, char **argv)
{
    Sequence::Fasta fseq; //declare a variable of type Sequence::Fasta
    while(!cin.eof()) //read through the file
    {
       cin >> fseq;
       fseq.Revcom(); //"Revcom" the sequence
       cout << fseq << endl; //write results to stdout
    }
    return 0;
}

