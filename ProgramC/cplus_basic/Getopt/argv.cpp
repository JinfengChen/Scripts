#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

using namespace std;

int main (int argc, char *argv[])
{
        char *chr; // chromosome
        int max_sd = 1000000000;
        int fisher = 0;
        
        if(argc < 3){
                fprintf(stderr, "Usage:\n");
                fprintf(stderr, "./argv chr01 1 200\n\n");
                return 1;
        }
                
        chr = argv[1];
        sscanf(argv[2],"%d",&fisher); // convert char to int: sscanf(source str, "type", address)
        sscanf(argv[3],"%d",&max_sd);
       
        cout<<"o "<<chr<<"\n"<<endl;
        cout<<"m "<<max_sd<<"\n"<<endl;
        cout<<"f "<<fisher<<"\n"<<endl;
 
}
