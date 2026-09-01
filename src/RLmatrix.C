/* Program file RLmatrix.C -
* Performing the matrix diagonalization
* Last change: KL 26.03.14
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include "Functions.h"
#include <RLexact.h>
#include <nr.h>
#include <cnr.h>

// long long CrossMatrix();   /* NOT IN USE !!! */

/* Global variables defined in RLexact.c */
extern unsigned long long *unique;
extern long long *uniq_k;
extern long long Nsym;
extern long long Nunique, Nuniq_k;
extern komplex *gs, *smpgs, *szqgs, *tmp;
// extern komplex **hamilton;
extern double *energies;
extern long long *mag;
extern double *magnetisation;
#ifdef CROSS
extern double *cross, *sz;
#endif /* CROSS */
extern double cosine[], sine[], invsqrt[];
extern long long Nocc[], Ncoup, Nspins;
extern long long m;

unsigned long long i_min;

// TODO: Find out what uniqk has been changed to ...

/* Matrix_gs(k) diagonalizes the Hamilton matrix for
k[] and given h/m and finds the ground state.
The routine returns the energy of the ground state,
and evec is the ground state vector */
double Matrix_gs(komplex **hamil, long long *uniqk, long long k[NSYM],
                 komplex *evec, struct FLAGS *input_flags)
{
  time_t time_single;
  double emin;
  long long i;
  long long j;

  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageChar(" Matrix_gs called: \n");
  }
  if (input_flags->VERBOSE_TIME_LV1)
  {
    time_stamp(&time_single, START, "\n timer of intiallizing and filling matrix");
  }
  BuildCycle(k, input_flags);
  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageCharInt(" Nuniq_k: ", Nuniq_k);
    hamil[1][1] = zero;
    LogMessageChar(" , hamil accessed, Fill zeros H=\n");
  }
  if (Nuniq_k == 0)
    return LARGE_NUMBER;
  for (i = 1; i <= Nuniq_k; i++)
  { // Nuniq_k
    for (j = 1; j <= Nuniq_k; j++)
    { // Nuniq_k
      hamil[i][j] = zero;
      if (input_flags->MATRIX_MESSAGES)
      {
        LogMessageCharDouble("  (", real(hamil[i][j]));
        LogMessageCharDouble("+ i", imag(hamil[i][j]));
        LogMessageChar(")");
      }
    }
    if (input_flags->MATRIX_MESSAGES)
    {
      LogMessageChar("\n");
    }
  }
  // FillHamilton(k, uniqk, hamil); Old Version
  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageChar("\n Next step, fill hamiltonian and diagonalize it \n H= \n");
  }
  FillHamilSparse(hamil, k, input_flags);
  if (input_flags->MATRIX_MESSAGES)
  {
    for (i = 1; i <= Nuniq_k; i++)
    {
      for (j = 1; j <= Nuniq_k; j++)
      {
        LogMessageCharDouble("  (", real(hamil[i][j]));
        LogMessageCharDouble("+ i", imag(hamil[i][j]));
        LogMessageChar(")");
      }
      LogMessageChar("\n");
    }
    LogMessageChar("\n");
  }
  if (input_flags->VERBOSE_TIME_LV1)
  {
    time_stamp(&time_single, STOP, " ");
    time_stamp(&time_single, START, "timer of total Diagonalization time");
  }

  Diagonalize(hamil, Nuniq_k, energies, input_flags);

  if (input_flags->VERBOSE_TIME_LV1)
  {
    time_stamp(&time_single, STOP, " Total diagonalization time for this symmetry=S");
  }

  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageChar("\n Hamiltonian has been diagonalized \nH = \n");
    for (i = 1; i <= Nuniq_k; i++)
    {
      for (j = 1; j <= Nuniq_k; j++)
      {
        LogMessageCharDouble("  (", real(hamil[i][j]));
        LogMessageCharDouble("+ i", imag(hamil[i][j]));
        LogMessageChar(")");
      }
      LogMessageChar("\n");
    }
    LogMessageChar("\n");
  }

  if (input_flags->find_mag)
  {
    CalculateMatrixM(hamil, magnetisation, input_flags);

    if (input_flags->TEST_MATRIXMAG)
    {
      LogMessageChar("Testing: \n");
      for (int i = 0; i < Nuniq_k; i++)
      {
        LogMessageCharDouble("magn = ", magnetisation[i]);
        LogMessageCharDouble(", E =", energies[i]);
        LogMessageChar("\n");
      }
    }
  }

  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageChar(" Matrix_gs, Hamilton diagonalized \n");
  }

  /* Find index of smallest eigenvalue */
  emin = LARGE_NUMBER;
  for (i = 0; i < Nuniq_k; i++)
  {
    if (input_flags->TEST_ENERGIES)
    {
      WriteEnergy(energies[i]);
    }
    if (energies[i] < emin)
    {
      emin = energies[i];
      i_min = i;
    }
  }

  if (input_flags->MATRIX_MESSAGES)
  {
    LogMessageCharDouble(" Smallest eigenvalue:", emin);
    LogMessageCharInt(" number:", i_min);
  }

  /* Now: Save ground state in evec  */
  for (i = 0; i < Nuniq_k; i++)
  {
    evec[i] = hamil[i + 1][i_min + 1];

    if (input_flags->TEST_EVEC)
    {
      LogMessageCharDouble(" ", real(evec[i]));
      LogMessageCharDouble(" + i ", imag(evec[i]));
      LogMessageChar("\n");
    }
  }

  if (input_flags->TEST_EIG)
  {
    Eigenvector_test(k, evec, tmp, input_flags);
  }
  return emin;
}

void CalculateMatrixM(komplex **matrix, double *mvector, struct FLAGS* input_flags)
{
  long long v;
  double sum;
  double norm;

  for (v = 0; v < Nuniq_k; v++)
  {
    // mvector[v]=CalculateM(matrix[v+1]); //this accessing the v+1'th row - we want the column!

    sum = 0;
    norm = 0;
    for (int j = 0; j < Nunique; j++)
    {
      sum += mag[uniq_k[j]] * sqrabs(matrix[j + 1][v + 1]);
      // assuming bloody vectors arent normalized
      norm += sqrabs(matrix[j + 1][v + 1]);
      if (input_flags->TEST_CALC_M)
      {
        LogMessageCharInt("j=", j);
        LogMessageCharDouble(", m_uniq=", mag[uniq_k[j]]);
        LogMessageCharDouble(", uniq_k=", uniq_k[j]);
        LogMessageCharDouble(", norm=", sqrabs(matrix[j + 1][v + 1]));
        LogMessageCharDouble(", matrix=", real(matrix[j + 1][v + 1]));
        LogMessageCharDouble("+ i", imag(matrix[v + 1][j + 1]));
        LogMessageChar("\n");
      }
    }

    if (input_flags->TEST_CALC_M)
    {
      LogMessageCharDouble("Sum = ", sum);
      LogMessageCharDouble("Norm = ", norm);
      LogMessageChar("\n");
    }
    mvector[v] = sum; // sqrt(norm);
  }

  return;
}

double CalculateM(komplex *state, struct FLAGS* input_flags) // currently obsolete SJ 100616
{

  long long j;
  double sum = 0;

  for (j = 0; j < Nunique; j++)
  {
    sum += mag[uniq_k[j]] * sqrabs(state[j + 1]);
    if (input_flags->TEST_CALC_M)
    {
      LogMessageCharInt("j=", j);
      LogMessageCharDouble(", m_uniq=", mag[uniq_k[j]]);
      LogMessageCharDouble(", uniq_k=", uniq_k[j]);
      LogMessageCharDouble(", norm=", sqrabs(state[j + 1]));
      LogMessageCharDouble(", state=", real(state[j + 1]));
      LogMessageCharDouble("+ i", imag(state[j + 1]));
      LogMessageChar("\n");
    }
  }

  if (input_flags->TEST_CALC_M)
  {
    LogMessageCharDouble("Sum = ", sum);
    LogMessageChar("\n");
  }

  return sum;
}
