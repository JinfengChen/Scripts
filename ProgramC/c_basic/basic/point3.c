#include "stdio.h"
void main()
{
  int m[5] = {5,4,3,2,1};
  int *p = m;
  printf ("%d,%d,%d",m[2],*(m+2),*(p+2));
}
