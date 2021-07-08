#include <stdio.h>
#include <stdlib.h>

// N is for series number, X is the input value
double Legendre(int N, double X)
{
  int n;
  double P0 = 1.0;
  double P1 = X;
  double P_b = -1;   // result, Pn+1
  double P_f = -1;   // Pn-1
  double P_ff = -1;  // Pn-2

  // check N
  if (N == 0)
  {
    return X;
  }
  else if (N < 0)
  {
    printf("ERROR: N must be positive...\n");
    return -1;
  }
  // calculate the Legendre number.
  for (n = 1; n < N; n++)
  {
    if (n == 1)
    {
      P_b = ((2 * n + 1)*X*P1 - n*P0) / (n + 1);
      P_f = P_b;
      P_ff = P1;
      continue;
    }
    P_b = ((2 * n + 1)*X*P_f - n*P_ff) / (n + 1);
    P_ff = P_f;
    P_f = P_b;
  }
  return P_b;
}

int main(void)
{
  // make a test
  printf("%.8f\n", Legendre(6, 1.98));
  return 0;
}
