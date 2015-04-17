#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

using namespace std;

int main(int argc, char *argv[])
{
        int c;
	char *chr; // chromosome
	int max_sd = 1000000000;
        int fisher = 1;
	while((c = getopt(argc, argv, "o:m:f:")) >= 0){
		switch(c) {
			case 'o': chr = strdup(optarg); break; //strdup copy the string address of "o" to chr
			case 'm': max_sd = atoi(optarg); break; //atoi convert the string of "m" to int
			case 'f': fisher = 1; break; // gave fisher a value 1
                        default: fprintf(stderr, "Unrecognized option '-%c'.\n", c); // 
		}
	}
	if(argc < 7){
		fprintf(stderr, "\n");
		fprintf(stderr, "getopt -o chr01 -f 1 -m 200\n\n");
		fprintf(stderr, "Options: \n");
		fprintf(stderr, "       -o STRING       operate on a single chromosome [all chromosome]\n");
		fprintf(stderr, "       -f INT          number [%d]\n", fisher);		 
		fprintf(stderr, "       -m INT          maximum SV size [%d]\n", max_sd);		 
		return 1;
	}
        cout<<"o "<<chr<<"\n"<<endl;
        cout<<"m "<<max_sd<<"\n"<<endl; 
        cout<<"f "<<fisher<<"\n"<<endl;
}
