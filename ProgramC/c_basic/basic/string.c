#include <stdio.h>
void main()
{
  char a[]="Tsinghua University", b[30];
  int i,j;
  for (i=0;*(a+i) != '\0';i++)
     *(b+i)=*(a+i);
  *(b+i) = '\0';
  printf("String a is:%s\n",a);
  printf("String b is:");
  printf("%c",b[1]);
  printf("\n");
  for (i=0;b[i] != '\0';i++)
     printf("%c",b[i]);
  printf("\n");
}
