#include "stdio.h"
void main()
{
  static int a[3][4] = {{1,2,3,4},{5,6,7,8},{9,10,11,12}};
  int *p,i,j,n;
  n=sizeof(a)/sizeof(int);
  printf("%d,\n",sizeof(a));
  printf("%d,\n",sizeof(int));
  for(p=a[0];p<a[0]+n;p++)
     printf("%d,",*p);
}
