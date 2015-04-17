/////////////////////////////////////////////////////////////////////////
//Copyright (C) 2003 Dekang Lin, lindek@cs.ualberta.ca
//
//Permission to use, copy, modify, and distribute this software for any
//purpose is hereby granted without fee, provided that the above
//copyright notice appear in all copies and that both that copyright
//notice and this permission notice appear in supporting documentation.
//No representations about the suitability of this software for any
//purpose is made. It is provided "as is" without express or implied
//warranty.
//
/////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "logprobs.h"

/** The input array contains a set of log probabilities lp1, lp2, lp3
    ... The return value should be the log of the sum of the
    probabilities: log(e^lp1 + e^lp2 + e^lp3 + ...) */
double sumLogProb(vector<double>& logprobs)
{
  double max = 0;
  unsigned int i;
  for (i = 0; i<logprobs.size(); i++) {
    if (i==0 || logprobs[i]>max)
      max = logprobs[i];
  }
  if (isinf(max)) // the largest probability is 0 (log prob= -inf)
    return max;   // return log 0
  double p = 0;
  for (i = 0; i<logprobs.size(); i++) {
    p += exp(logprobs[i]-max);
  }
  return max + log(p);
}

/** returns log (e^logprob1 + e^logprob2). */
double sumLogProb(double logprob1, double logprob2)
{
  if (isinf(logprob1) && isinf(logprob2)) 
    return logprob1; // both prob1 and prob2 are 0, return log 0.
  if (logprob1>logprob2)
    return logprob1+log(1+exp(logprob2-logprob1));
  else
    return logprob2+log(1+exp(logprob1-logprob2));
}

