/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   4 May 1999 
**      File:   testfor.c
**      Purpose: driver for testing the Forward, ForwardWithScale code.
**      Organization: University of Maryland
**
**	$Id$
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

int main (int argc, char **argv)
{
	int 	t, T; 
	HMM  	hmm;
	int	*O;	/* observation sequence O[1..T] */
	double **alpha;
	double *scale;
	double 	proba, logproba; 
	FILE	*fp;

	if (argc != 3) {
		printf("Usage error \n");
		printf("Usage: testfor <model.hmm> <obs.seq> \n");
		exit (1);
	}
	
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}
	ReadHMM(fp, &hmm);
	fclose(fp);

	fp = fopen(argv[2], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[2]);
		exit (1);
	}
	ReadSequence(fp, &T, &O);
	fclose(fp);


	alpha = dmatrix(1, T, 1, hmm.N);
	scale = dvector(1, T);

	printf("------------------------------------\n");
	printf("Forward without scaling \n");
	Forward(&hmm, T, O, alpha, &proba); 
	fprintf(stdout, "log prob(O| model) = %E\n", log(proba));

	printf("------------------------------------\n");
	printf("Forward with scaling \n");

	ForwardWithScale(&hmm, T, O, alpha, scale, &logproba); 

	fprintf(stdout, "log prob(O| model) = %E\n", logproba);
	printf("------------------------------------\n");
	printf("The two log probabilites should identical \n");
	printf("(within numerical precision). When observation\n");
	printf("sequence is very large, use scaling. \n");
	
	free_ivector(O, 1, T);
	free_dmatrix(alpha, 1, T, 1, hmm.N);
	free_dvector(scale, 1, T);
	FreeHMM(&hmm);
}

