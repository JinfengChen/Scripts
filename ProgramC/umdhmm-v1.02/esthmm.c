/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   22 February 1988 
**      File:   esthmm.c
**      Purpose: estimate HMM parameters from observation. 
**      Organization: University of Maryland
**
**      $Id: esthmm.c,v 1.1 1998/02/23 07:49:45 kanungo Exp kanungo $
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>


static char rcsid[] = "$Id: esthmm.c,v 1.1 1998/02/23 07:49:45 kanungo Exp kanungo $";

void Usage(char *name);

int main (int argc, char **argv)
{
	int 	T;
	HMM  	hmm;
	int	N;
	int	M;
	double 	**alpha; 
	double	**beta;
	double	**gamma;
	int	*O;
	int	iflg=0, sflg=0, nflg=0, mflg=0, errflg =0, vflg=0;
	int	c;
	int	seed; /* seed for random number generator */
	char	*hmminitfile;
	int	niter;
	double	logprobinit, logprobfinal;
	FILE	*fp;
        extern char *optarg;
        extern int optind, opterr, optopt;
 

        while ((c= getopt(argc, argv, "vhI:S:N:M:")) != EOF)
                switch (c) {
		case 'v': 
			vflg++;	
			break;
		case 'h': 
			Usage(argv[0]);
			exit(1);
			break;
                case 'S':
                        /* set random number generator seed */
                        if (sflg)
                                errflg++;
                        else {
                                sflg++;
                                sscanf(optarg, "%d", &seed);
                        }
                        break;
                case 'N':  
                        /* set random number generator seed */
                        if (nflg) 
                                errflg++; 
                        else { 
                                nflg++;  
                                sscanf(optarg, "%d", &N);
                        } 
                        break;   
                case 'M':  
                        /* set random number generator seed */
                        if (mflg) 
                                errflg++; 
                        else { 
                                mflg++;  
                                sscanf(optarg, "%d", &M);
                        } 
                        break;   
                case 'I':  
                        /* set random number generator seed */
                        if (iflg) 
                                errflg++; 
                        else { 
                                iflg++;  
				hmminitfile = optarg;
                        } 
                        break;   
                case '?':
                        errflg++;
                }

	/* you can initialize the hmm model three ways:
           i) with a model stored in a file, which also sets 
	      the number of states N and number of symbols M.
           ii) with a random model by just specifyin N and M
              on the command line.
           iii) with a specific random model by specifying N, M
              and seed on the command line. 
        */

	if (iflg) {
		/* model being read from a file */
		if (((sflg || nflg) || mflg)) errflg++;
	}
	else if ((!nflg) || (!mflg)) { 
		/* Model not being intialied from file */ 
		/* both N and M should be specified */
		errflg++; 
	}

	
        if ((argc - optind) != 1) errflg++; /* number or arguments not okay */



	if (errflg) {
		Usage(argv[0]);
		exit (1);
	}
		
	/* read the observed sequence */
	fp = fopen(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found \n", argv[optind]);
		exit (1);
	}
	ReadSequence(fp, &T, &O); 
	fclose(fp);


	/* initialize the hmm model */
	if (iflg) { 
		fp = fopen(hmminitfile, "r");
		if (fp == NULL) {
			fprintf(stderr, "Error: File %s not found \n", 
				hmminitfile);
                	exit (1);
		}
		ReadHMM(fp, &hmm);
		fclose(fp);
	}
	else if (sflg)  
		InitHMM(&hmm, N, M, seed);
	else {
		seed = hmmgetseed();
		InitHMM(&hmm, N, M, seed);
	}

	/* allocate memory */
	alpha = dmatrix(1, T, 1, hmm.N);
	beta = dmatrix(1, T, 1, hmm.N);
	gamma = dmatrix(1, T, 1, hmm.N);

	/* call Baum Welch */
	BaumWelch(&hmm, T, O, alpha, beta, gamma, &niter, 
		&logprobinit, &logprobfinal);

	if (vflg) {
		if (sflg) fprintf(stderr, "RandomSeed: %d\n", seed);
		fprintf(stderr, "Number of iterations: %d\n", niter);
		fprintf(stderr, "Log Prob(observation | init model): %E\n",
			logprobinit);	
		fprintf(stderr, "Log Prob(observation | estimated model): %E\n",
			logprobfinal);	
	}


	/* print the answer */
	PrintHMM(stdout, &hmm);

	/* free memory */
	free_ivector(O, 1, T);
	free_dmatrix(alpha, 1, T, 1, hmm.N);
	free_dmatrix(beta, 1, T, 1, hmm.N);
	free_dmatrix(gamma, 1, T, 1, hmm.N);
	FreeHMM(&hmm);
}

void Usage(char *name)
{
	printf("Usage error. \n");
        printf("Usage1: %s [-v] -N <num_states> -M <num_symbols> <file.seq>\n", 
		name);
        printf("Usage2: %s [-v] -S <seed> -N <num_states> -M <num_symbols> <file.seq>\n", 
		name);
        printf("Usage3: %s [-v] -I <mod.hmm> <file.seq>\n", 
		name);
        printf("  N - number of states\n");
        printf("  M - number of symbols\n");
        printf("  S - seed for random number genrator\n");
        printf("  I - mod.hmm is a file with the initial model parameters\n");
        printf("  file.seq - file containing the obs. seqence\n");
	printf("  v - prints out number of iterations and log prob\n");
}
