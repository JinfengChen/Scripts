/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   backward.c
**      Purpose: Backward algorithm for computing the probabilty
**              of observing a sequence given a HMM model parameter.
**      Organization: University of Maryland
**
**      $Id: backward.c,v 1.3 1998/02/23 07:56:05 kanungo Exp kanungo $
*/

#include <stdio.h>
#include "hmm.h"
static char rcsid[] = "$Id: backward.c,v 1.3 1998/02/23 07:56:05 kanungo Exp kanungo $";

void Backward(HMM *phmm, int T, int *O, double **beta, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
        double sum;
 
 
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
                beta[T][i] = 1.0;
 
        /* 2. Induction */
 
        for (t = T - 1; t >= 1; t--) {
                for (i = 1; i <= phmm->N; i++) {
                        sum = 0.0;
                        for (j = 1; j <= phmm->N; j++)
                                sum += phmm->A[i][j] *
                                        (phmm->B[j][O[t+1]])*beta[t+1][j];
                        beta[t][i] = sum;
 
                }
        }
 
        /* 3. Termination */
        *pprob = 0.0;
        for (i = 1; i <= phmm->N; i++)
                *pprob += beta[1][i];
 
}

void BackwardWithScale(HMM *phmm, int T, int *O, double **beta, 
	double *scale, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
	double sum;
 
 
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
                beta[T][i] = 1.0/scale[T]; 
 
        /* 2. Induction */
 
        for (t = T - 1; t >= 1; t--) {
                for (i = 1; i <= phmm->N; i++) {
			sum = 0.0;
                        for (j = 1; j <= phmm->N; j++)
                        	sum += phmm->A[i][j] * 
					(phmm->B[j][O[t+1]])*beta[t+1][j];
                        beta[t][i] = sum/scale[t];
 
                }
        }
 
}
