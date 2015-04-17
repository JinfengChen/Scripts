/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   baumwelch.c
**      Purpose: Baum-Welch algorithm for estimating the parameters
**              of a HMM model, given an observation sequence. 
**      Organization: University of Maryland
**
**	Update: 
**	Author: Tapas Kanungo
**	Date:	19 April 1999
**	Purpose: Changed the convergence criterion from ratio
**		to absolute value. 
**
**      $Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $
*/

#include <stdio.h> 
#include "nrutil.h"
#include "hmm.h"
#include <math.h>

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001 
void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,
	double **gamma, int *pniter, 
	double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb,  threshold;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;

	deltaprev = 10e-70;

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
	*plogprobinit = logprobf; /* log P(O |intial model) */
	BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
	ComputeGamma(phmm, T, alpha, beta, gamma);
	ComputeXi(phmm, T, O, alpha, beta, xi);
	logprobprev = logprobf;

	do  {	

		/* reestimate frequency of state i in time t=1 */
		for (i = 1; i <= phmm->N; i++) 
			phmm->pi[i] = .001 + .999*gamma[1][i];

		/* reestimate transition matrix  and symbol prob in
		   each state */
		for (i = 1; i <= phmm->N; i++) { 
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++) 
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++) 
					numeratorA += xi[t][i][j];
				phmm->A[i][j] = .001 + 
						.999*numeratorA/denominatorA;
			}

			denominatorB = denominatorA + gamma[T][i]; 
			for (k = 1; k <= phmm->M; k++) {
				numeratorB = 0.0;
				for (t = 1; t <= T; t++) {
					if (O[t] == k) 
						numeratorB += gamma[t][i];
				}

				phmm->B[i][k] = .001 +
						.999*numeratorB/denominatorB;
			}
		}

		ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
		BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
		ComputeGamma(phmm, T, alpha, beta, gamma);
		ComputeXi(phmm, T, O, alpha, beta, xi);

		/* compute difference between log probability of 
		   two iterations */
		delta = logprobf - logprobprev; 
		logprobprev = logprobf;
		l++;
		
	}
	while (delta > DELTA); /* if log probability does not 
                                  change much, exit */ 
 
	*pniter = l;
	*plogprobfinal = logprobf; /* log P(O|estimated model) */
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta, 
	double **gamma)
{

	int 	i, j;
	int	t;
	double	denominator;

	for (t = 1; t <= T; t++) {
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}

		for (i = 1; i <= phmm->N; i++) 
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

void ComputeXi(HMM* phmm, int T, int *O, double **alpha, double **beta, 
	double ***xi)
{
	int i, j;
	int t;
	double sum;

	for (t = 1; t <= T - 1; t++) {
		sum = 0.0;	
		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha[t][i]*beta[t+1][j]
					*(phmm->A[i][j])
					*(phmm->B[j][O[t+1]]);
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) 
				xi[t][i][j]  /= sum;
	}
}

double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N, 1, N);
	return xi;
}

void FreeXi(double *** xi, int T, int N)
{
	int t;



	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N, 1, N);

	xi ++;
	free(xi);

}
