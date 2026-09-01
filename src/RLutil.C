/* Program file RLutil.C
* Different general routines
* Last change: KL 20.03.15
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94
* Version 4.0, September 2017
*
============================================
*/
// #include "/usr/include/sys/types.h"
#include <sys/types.h>

#include <string.h>
#include <math.h>
#include <complex>

#include <RLexact.h>
#include "Functions.h"

/* Global variables declared elsewhere */

extern long long Nunique;

/* Functions defined elsewhere */

#ifndef NEVER
double rot_mat[3][3];
#endif /* M_SYM */

/* reverse string */
void reverse(char s[])
{
  long long c, i, j;

  for (i = 0, j = strlen(s) - 1; i < j; i++, j--)
  {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

/* itoa convert n to characters in s */
void itoa(long long n, char s[])
{
  long long i, sign;

  if ((sign = n) < 0)
    n = -n;
  i = 0;
  do
  {
    s[i++] = n % 10 + '0';
  } while ((n /= 10) > 0);
  if (sign < 0)
    s[i++] = '-';
  s[i] = '\0';
  reverse(s);
}

#ifndef NEVER
void FillRotationMatrix(double *vector, struct FLAGS *input_flags)
{
  double costh, sinth, cosphi, sinphi;

  costh = vector[Z];
  sinth = sqrt(1 - SQR(costh));
  if (sinth == 0)
  {
    cosphi = 1;
    sinphi = 0;
  }
  else
  {
    cosphi = vector[X] / sinth;
    sinphi = vector[Y] / sinth;
  }
if (input_flags->TEST_ROTATION)
{
  LogMessageCharInt("Rotation angles: cos(th): ", costh);
  LogMessageCharInt("Rotation angles: sin(th): ", sinth);
  LogMessageCharInt("Rotation angles: cos(phi): ", cosphi);
  LogMessageCharInt("Rotation angles: sin(phi): ", sinphi);
  LogMessageChar("\n");
// SQR(cosphi)+SQR(sinphi)
}
  rot_mat[X][X] = cosphi * costh;
  rot_mat[X][Y] = sinphi * costh;
  rot_mat[X][Z] = -sinth;
  rot_mat[Y][X] = -sinphi;
  rot_mat[Y][Y] = cosphi;
  rot_mat[Y][Z] = 0;
  rot_mat[Z][X] = sinth * cosphi;
  rot_mat[Z][Y] = sinth * sinphi;
  rot_mat[Z][Z] = costh;
if (input_flags->TEST_ROTATION)
{
  LogMessageChar3Vector("The rotation matrix is: \n", rot_mat[X][X], rot_mat[X][Y], rot_mat[X][Z]);
  LogMessageChar3Vector("\n", rot_mat[Y][X], rot_mat[Y][Y], rot_mat[Y][Z]);
  LogMessageChar3Vector("\n", rot_mat[Z][X], rot_mat[Z][Y], rot_mat[Z][Z]);
  LogMessageChar("\n");
  double test_vec[3];
  test_vec[X] = vector[X];
  test_vec[Y] = vector[Y];
  test_vec[Z] = vector[Z];
  RotateVector(test_vec);
  LogMessageChar3Vector("Rotation of vector", vector[X], vector[Y], vector[Z]);
  LogMessageChar3Vector(" gives ", test_vec[X], test_vec[Y], test_vec[Z]);
  LogMessageChar("\n");
}
}
#endif /* M_SYM */

#ifndef NEVER
void NormalizeVector(double *vector)
{
  double r = sqrt(SQR(vector[X]) + SQR(vector[Y]) + SQR(vector[Z]));

  if (r == 0)
    return;

  vector[X] = vector[X] / r;
  vector[Y] = vector[Y] / r; /* is understood as r_ij where i<j */
  vector[Z] = vector[Z] / r;

  return;
}

void RotateVector(double *vector)
{
  double vx = vector[X], vy = vector[Y], vz = vector[Z];
  long long i;

  for (i = X; i <= Z; i++)
  {
    vector[i] = rot_mat[i][X] * vx + rot_mat[i][Y] * vy + rot_mat[i][Z] * vz;
  }
  return;
}
#endif /* M_SYM */

void Bubblesort(double *key, double *value, long long length)
{
  // quick and dirty sort algorithm. Should not be used for large datasets

  double temp;

  for (long long i = 0; i < length; i++)
  {
    for (long long j = i; j < length; j++)
    {
      if (key[j] < key[i])
      {
        temp = key[i];
        key[i] = key[j];
        key[j] = temp;
        if (value != NULL)
        {
          temp = value[i];
          value[i] = value[j];
          value[j] = temp;
        }
      }
    }
  }
}
