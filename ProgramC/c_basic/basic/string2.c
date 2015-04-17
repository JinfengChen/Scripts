#include <stdio.h>
void main()
{
  char a[]="Tsinghua University", b[30],*p1,*p2;
  p1=a;
  p2=b;
  for(;*p1 != '\0'; p1++,p2++)
     *p2 = *p1;
  *p2 = '\0';
  printf("String a is:%s\n",a);
  printf("String b is:%s\n",b);
}

