#include <stdio.h>
void hanoi (int n, char A, char B, char C)
{
   if (n == 1)
   {
      printf("Move disk %d from %c to %c\n",n,A,C);
   }else
   {
      hanoi(n-1,A,C,B);
      printf("Move disk %d from %c to %c\n",n,A,C);
      hanoi(n-1,B,A,C);
   }
}

void main()
{
   int n;
   scanf("%d",&n);
   hanoi(n,'A','B','C');
}
