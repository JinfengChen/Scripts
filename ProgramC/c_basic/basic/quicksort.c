#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void PrintArray (char *a[], int from, int to)
{
   int i;
   for(i=from; i<= to;i++)
      printf ("%4s",a[i]);
   printf("\n");
}

void QuickSort (char *a[], int left, int right)
{
   int min, max;
   char *norm, *tmp;
   min = left;
   max = right;
   norm = a[(left+right)/2];
   do
   {
      while(strcmp(a[min],norm)<0 && min<right)
        min++;
      while(strcmp(a[max],norm)>0 && max>left)
        max--;
      if(min<=max)
      {
         tmp=a[min];a[min]=a[max];a[max]=tmp;
         min++;
         max--;
      }
      PrintArray(a,left,right);
   }while(min<=max);
   if (left<max)
      QuickSort(a,left,max);
   if (right>min)
      QuickSort(a,min,right);
}



void main ()
{
   char * a[]={"abc","ccc","bca","acb","bbb"};
   int n=5;
   PrintArray(a,0,n-1);
   QuickSort(a,0,n-1);
   PrintArray(a,0,n-1);
}
