#include "stdio.h"
void main ()
{ 
   int m,n,* pm,* pn, * temp;
   pm=&m;
   pn=&n;
   scanf("%d%d",&m,&n);
   printf ("%d,%d\n",m,n);
   temp = pm;
   pm = pn;
   pn = temp;
   printf ("%d,%d\n",m,n);
   printf ("%d,%d\n",* pm,* pn);
}
