/* Program file RLhamil.C -
* Applying the Hamilton operator to a given state
* Last change: KL 25.09.12
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// #include <nr.h>
#include <cnr.h>
#include <RLexact.h>
#include <complex>
#include "Functions.h"
#include "Max_problem.h"

/* Global variables defined in RLexact.c */
extern unsigned long long *unique;
extern long long Nsym, Nspins;
extern long long Nunique, Nuniq_k;
extern long long *uniq_k;
// extern long long *Nocc;
extern long long Nsymvalue[];
extern long long hamil_coup[NCOUP][2];
extern double Jzz[], Jxy[], Janis[];
extern double Jdip[], geom_13[], r_vector[NCOUP][3];
extern long long m;
extern double h, field[3];
extern double sine[], cosine[], sqroot[];
extern long long Ncoup;
extern double Jr[NRING];
extern long long ring_coup[NRING][4];

/* Regional variables in this file */
unsigned long long bitmap, new_state;
long long n_2, u_occ;
unsigned long long index1, index2;
komplex this_;

void Hamil_Zeeman(unsigned long long bitmap, unsigned long long *new_state, long long i, int *nelem, long long *totcount, int *T, komplex *J, double *diag, FILE *indexfile, FILE *Tfile, FILE *Jfile, struct FLAGS *input_flags)
{ // only supports fields along one coordinate axis
  if (input_flags->TEST_HAMZEE)
  {
    LogMessageCharDouble("\nIn Hamil_Zeeman. h =", h);
    LogMessageCharInt(" bitmap =", bitmap);
  }
  if (h != 0) // if there is no field, this doesnt matter
  {

    if (input_flags->TEST_HAMZEE)
    {
      LogMessageCharDouble("\nfield[0] =", field[0]);
      LogMessageCharDouble(", field[1] =", field[1]);
      LogMessageCharDouble(", field[2] =", field[2]);
    }
    unsigned long long mask0, s0;
    if (field[2] > 0) /* Mag field along z*/
    {
      double sz;
      sz = Count(bitmap, input_flags) - Nspins / 2;
      if (input_flags->TEST_HAMZEE)
      {
        LogMessageCharInt("\nHamzee Sz: Count(bitmap) =", Count(bitmap, input_flags));
        LogMessageCharInt(" , and sz =", sz);
        LogMessageChar("\n");
      }
      *diag -= h * sz * field[2];
    }

    if (field[0] > 0 || field[1] > 0) // Magnetic field "transverse", i.e. along x or y
    {                                 // Here J is simply the coupling strength
      for (int k = 0; k < Nspins; k++)
      {
        mask0 = ((unsigned long long)1) << hamil_coup[k][0];
        s0 = (bitmap & mask0) != 0;
        if (input_flags->TEST_HAMZEE)
        {
          LogMessageCharInt("k =", k);
          LogMessageCharInt(" and s0 =", s0);
          LogMessageChar("\n");
        }

        if (s0 == 0) // S+, spin can be raised
        {
          *new_state = (bitmap | mask0);
          *J = -h / 2 * field[0] - (h / 2 * I) * field[1];
          if (input_flags->TEST_HAMZEE)
          {
            LogMessageCharInt("Hamzee: S+: From state ", bitmap);
            LogMessageCharInt("to state ", *new_state);
            LogMessageCharDouble("with J=", real(*J));
            LogMessageCharDouble("+ i", imag(*J));
            LogMessageChar("\n");
          }
        }
        else // then s0==1, S-, spin can be lowered
        {
          *new_state = (bitmap & ~(mask0));
          *J = -h / 2 * field[0] + (h / 2 * I) * field[1];

          if (input_flags->TEST_HAMZEE)
          {
            LogMessageCharInt("Hamzee: S-: From state ", bitmap);
            LogMessageCharInt("to state ", *new_state);
            LogMessageCharDouble("with J=", real(*J));
            LogMessageCharDouble("+ i", imag(*J));
            LogMessageChar("\n");
          }
        }
        WriteCouplingFiles(bitmap, *new_state, T, *J, i, i, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
      }
    }
  }
}

void Hamil2_sparse(unsigned long long bitmap, unsigned long long *new_state,
                   long long i, long long j, int *nelem, long long *totcount,
                   int *T, komplex *J, double *diag, FILE *indexfile,
                   FILE *Tfile, FILE *Jfile, struct FLAGS *input_flags)
{
  if (input_flags->TEST_HAM2)
    LogMessageChar("Now entering Hamil2_sparse function \n");
  unsigned long long mask0, mask1, s0, s1;

  mask0 = ((unsigned long long)1) << hamil_coup[j][0];
  mask1 = ((unsigned long long)1) << hamil_coup[j][1];
  s0 = (bitmap & mask0) != 0;
  s1 = (bitmap & mask1) != 0;
  if (input_flags->TEST_HAM2)
  {
    LogMessageCharInt(" Coupling between ", hamil_coup[j][0]);
    LogMessageCharInt(" and ", hamil_coup[j][1]);
    LogMessageCharInt(", (s0,s1)= (", s0);
    LogMessageCharInt(", ", s1);
    LogMessageChar(")\n");
    LogMessageCharInt("Also, mask0 = ", mask0);
    LogMessageCharInt(" and mask1 = ", mask1);
  }
  if ((s0 + s1) == 1) /* Spins are of opposite sign */
    *diag -= Jzz[j] * 0.25;
  else
    *diag += Jzz[j] * 0.25;

  if (s0 == 0) /* case s0 down; S+.. terms (and Sz.. for DIPOLE) */
  {

    if (s1 == 0) /* down down: S+S+ term */
    {
      if (!input_flags->m_sym)
      {
        *new_state = (bitmap | mask0 | mask1);
        *J = Janis[j] / 2;
      }
    }
    else /* down up: S+S- terms */
    {
      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jxy[j] / 2.0;
    } /* if s1==0. */
  }
  else /* case s0 up; S-.. terms (and Sz.. for DIPOLE) */
  {

    if (s1 == 0) /* up down: S-S+ term (and SzS+) */
    {
      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jxy[j] / 2.0;
    }
    else /* up up: S-S- terms (and SzS-) */
    {
      if (!input_flags->m_sym)
      {
        *new_state = (bitmap & ~mask0) & ~mask1;
        *J = Janis[j] / 2;
      }
    } /* if s1==0.. */
  } /* if s0==0.. */

  if (input_flags->TEST_HAM2)
  {
    LogMessageCharDouble("with J=", real(*J));
    LogMessageCharDouble("+ i", imag(*J));
    LogMessageChar("\n");
  }

  WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
}

void Hamil4_sparse(unsigned long long bitmap, unsigned long long *new_state,
                   long long i, long long j, int *nelem, long long *totcount,
                   int *T, komplex *J, double *diag,
                   FILE *indexfile, FILE *Tfile, FILE *Jfile, struct FLAGS *input_flags)
{
  unsigned long long mask0, mask1, mask2, mask3, s0, s1, s2, s3;

  mask0 = ((unsigned long long)1) << ring_coup[j][0];
  mask1 = ((unsigned long long)1) << ring_coup[j][1];
  mask2 = ((unsigned long long)1) << ring_coup[j][2];
  mask3 = ((unsigned long long)1) << ring_coup[j][3];

  s0 = (bitmap & mask0) != 0;
  s1 = (bitmap & mask1) != 0;
  s2 = (bitmap & mask2) != 0;
  s3 = (bitmap & mask3) != 0;

  if (input_flags->TEST_HAM4)
  {
    LogMessageChar(" Ring coupling: ");
    LogMessageInt(s0);
    LogMessageInt(s1);
    LogMessageInt(s2);
    LogMessageInt(s3);
    LogMessageChar("\n");
  }

  if (s0 + s1 + s2 + s3 == 0 || s0 + s1 + s2 + s3 == 4 || s0 + s1 + s2 + s3 == 2) /* Hzzzz */
  {
    *diag += Jr[j] / 16.0;
    if (input_flags->TEST_HAM4)
    {
      LogMessageCharDouble("  Hzzz term, J = ", Jr[j] / 16.0);
      LogMessageChar("\n");
    }
  }
  else
  {
    *diag -= Jr[j] / 16.0;
    if (input_flags->TEST_HAM4)
    {
      LogMessageCharDouble("  Hzzz term, J = ", -Jr[j] / 16.0);
      LogMessageChar("\n");
    }
  }

  if (s0 + s1 + s2 + s3 == 2) /* m=0 rings*/
  {
    if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 1) /*Hpmpm, only alternating rings*/
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hpmpm term");
        LogMessageChar("\n");
      }
      *new_state = ((((bitmap | mask0) & ~mask1) | mask2) & ~mask3);
      *J = Jr[j] / 2.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 0)
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hpmpm term");
        LogMessageChar("\n");
      }
      *new_state = ((((bitmap | mask1) & ~mask0) | mask3) & ~mask2);
      *J = Jr[j] / 2.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 0 && s1 == 0 && s2 == 1 && s3 == 1) /*Hzpm down down up up*/
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask1) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 1 && s1 == 1 && s2 == 0 && s3 == 0) /*Hzpm up up down down*/
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask3) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 0 && s1 == 1 && s2 == 1 && s3 == 0) /*Hzpm down up up down*/
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }
      *new_state = ((bitmap | mask3) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 1 && s1 == 0 && s2 == 0 && s3 == 1) /*Hzpm up down down  up*/
    {
      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }
      *new_state = ((bitmap | mask1) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 0) /*Hzpm up down up down */
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }
      *new_state = ((bitmap | mask3) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 1) /*Hzpm down up down up */
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
  } /* if m=0 ring */

  if (s0 + s1 + s2 + s3 == 1) /* if m=-1 ring */
  {
    if (s0 == 1) /*Hzpm up down down down*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s1 == 1) /*Hzpm down up down down*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s2 == 1) /*Hzpm down down up down*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }
      *new_state = ((bitmap | mask1) & ~mask2);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s3 == 1) /*Hzpm down down down up*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask0) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
  } /* if m=-1*/

  if (s0 + s1 + s2 + s3 == 3) /* if m=1 ring */
  {
    if (s0 == 0) /*Hzpm down up up up*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s1 == 0) /*Hzpm up down up up*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s2 == 0) /*Hzpm up up down up*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask2) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
    if (s3 == 0) /*Hzpm up up up down*/
    {

      if (input_flags->TEST_HAM4)
      {
        LogMessageChar("  Hzzpm term");
        LogMessageChar("\n");
      }

      *new_state = ((bitmap | mask3) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);

      *new_state = ((bitmap | mask3) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile, input_flags);
    }
  } /*for m=1 rings*/
}

void Eigenvector_test(long long k[NSYM], komplex *evec, komplex *tmp, struct FLAGS *input_flags)
{
  long long i;
  komplex product;

  LogMessageChar("Entering Eigenvector_test.\n");
  // Hamilton(evec, tmp, k);
  // DLC: Seems like there was an old Hamilton function no longer in use here.
  product = zero;
  LogMessageChar("Middle of Eigenvector_test.\n");
  for (i = 0; i < Nunique; i++)
    product += evec[i] * conj(tmp[i]);
  WriteGSEnergy(product);

  if (input_flags->write_states)
  {
    WriteState("Eigenvector:", evec);
    WriteState("H|e>:", tmp);
  }

  return;
}
