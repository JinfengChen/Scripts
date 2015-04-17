#include "stdio.h"
void main()
{
   float data[5], *pf=data;
   int k;
   printf("input datas:\n");
   for (k=0;k<5;k++)
       scanf("%f",pf+k);
   for (k=0;k<5;k++)
   {
       printf("%f.1f,",*pf);
       pf++;
   }
   pf=data;
   printf("\n");
   for (k=0;k<5;k++)
       printf("%f.1f,",* (pf+k));
   printf("\n");

}
