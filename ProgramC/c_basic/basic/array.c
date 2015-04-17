#include <stdio.h>
void main()
{
    int f[11], k;
    f[0]=0;f[1]=1;f[2]=2;
    for (k=3;k<11;k++){
        f[k]=f[k-1]+2*f[k-2]*f[k-3];
    }
    for (k=0;k<11;k++){
        printf ("%d,\n",f[k]);
    }
}

