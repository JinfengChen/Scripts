/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   22 February 1998 
**      File:   sequence.c
**      Purpose: Routines for generating, reading and writing sequence of
**		observation symbols.
**      Organization: University of Maryland
**	
**	Update:
**	Author: Tapas Kanungo
**	Purpose: To make calls to generic random number generators
**		and to change the seed everytime the software is executed.
**
**      $Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"

static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q)
{
        int     t = 1;
        int     q_t, o_t;

	hmmsetseed(seed); 
 
        q[1] = GenInitalState(phmm);
        O[1] = GenSymbol(phmm, q[1]);
 
        for (t = 2; t <= T; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                O[t] = GenSymbol(phmm, q[t]);
        }
}

int GenInitalState(HMM *phmm)
{
        double val, accum;
        int i, q_t;
 
        val = hmmgetrand();
        accum = 0.0;
        q_t = phmm->N;
        for (i = 1; i <= phmm->N; i++) {
                if (val < phmm->pi[i] + accum) {
                        q_t = i;
                        break;
                }
                else {
                                accum += phmm->pi[i];
                }
        }
 
        return q_t;
}

int GenNextState(HMM *phmm, int q_t)
{
        double val, accum;
        int j, q_next;
 
        val = hmmgetrand();
        accum = 0.0;
        q_next = phmm->N;
        for (j = 1; j <= phmm->N; j++) {
                if ( val < phmm->A[q_t][j] + accum ) {
                        q_next = j;
                        break;
                }
                else
                        accum += phmm->A[q_t][j];
        }
 
        return q_next;
}
int GenSymbol(HMM *phmm, int q_t)
{
        double val, accum;
        int j, o_t;
 
        val = hmmgetrand();
        accum = 0.0;
        o_t = phmm->M;
        for (j = 1; j <= phmm->M; j++) {
                if ( val < phmm->B[q_t][j] + accum ) {
                       o_t = j;
                       break;
                }
                else
                        accum += phmm->B[q_t][j];
        }
 
        return o_t;
}
 
void ReadSequence(FILE *fp, int *pT, int **pO)
{
        int *O;
        int i;
 
        fscanf(fp, "T= %d\n", pT);
        O = ivector(1,*pT);
        for (i=1; i <= *pT; i++)
                fscanf(fp,"%d", &O[i]);
        *pO = O;
}
 
void PrintSequence(FILE *fp, int T, int *O)
{
        int i;
 
        fprintf(fp, "T= %d\n", T);
        for (i=1; i <= T; i++) 
                fprintf(fp,"%d ", O[i]);
	printf("\n");
 
}

