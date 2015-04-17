#include <stdio.h>
#define MAXN 7
#define NUM 29
int num[NUM];
int lottery[MAXN];
void combine (int n, int m)
{
     int i,j;
     for (i=n;i>=m;i--)
     {
         lottery[m-1]=num[i-1];
         printf("i:%d\tm:%d\t",i,m);
         if (m>1)
            combine(i-1,m-1);
         else
         {
            for (j=MAXN-1;j>=0;j--)
                printf("%3d",lottery[j]);
            printf ("\n");
         }
     }
}

int main()
{
    int i,j;
    for(i=0;i<NUM;i++)
       num[i]=i+1;
    for(i=i;i<MAXN;i++)
       lottery[i]=0;
    combine(NUM,MAXN);
    return 0;
}
