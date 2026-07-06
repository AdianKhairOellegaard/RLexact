/* Program file RLlancz.C -
* Performing the Lanczos diagonalization
* Last change: KL 06.04.15
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

#include <errno.h>
#include <string.h>
#include <RLexact.h>
#include "Functions.h"
#include <nr.h>
#include <cnr.h>

/* Functions defined in this file */

long long CrossLanczos(long long *);

void MakeSeed(komplex *);
void MakeSeedCross(komplex *, long long);
#ifdef FIND_MAG
double findmag(komplex *);
void findmaggs();
#endif /* FIND_MAG*/
double Normalize(komplex *);

/* Functions declared elsewhere */
extern void WriteEnergy(double);
extern void WriteState(const char *, double **);
extern void ApplySparse(komplex *vectin, komplex *vectout, long long *k);
extern void fatalerror(const char *, long long);
extern void Warning(const char *, long long);
extern void LogMessageChar(const char *);
extern void LogMessageInt(long long);
extern void LogMessageCharDouble(const char *, double);
extern void LogMessageCharInt(const char *, long long);
extern void LogMessageChar3Vector(const char *, double, double, double);
extern void OutMessageChar(const char *);
extern long long htqli(double *, double *, long long, komplex **);

/* Global variables declared in RLexact.c */
extern unsigned long long *unique;
extern long long Nunique;
extern long long *mag;
#ifdef FIND_MAG
extern double *magnetisation;
#endif // FIND_MAG
extern long long *Nocc;
extern komplex *vec1, *vec2, *vec3;
extern komplex *gs, *szxygs;
extern double *energies, *cross, *sz;
extern double cosine[], sine[], invsqrt[];
extern long long Ncoup, Nspins;
extern long long m;
extern long Nsym;
#ifdef FIND_MAG
extern double maggs;
#endif // end M_SYM

extern double Ritz_conv;

/* Regional variables in this file */
double length;
unsigned long long r, imin;
komplex **z, *first, *second, *third, *temp;
double diag[MAX_LANCZOS], subdiag[MAX_LANCZOS];
double diag_copy[MAX_LANCZOS], subdiag_copy[MAX_LANCZOS];
#ifdef FIND_MAG
double lanczmag[MAX_LANCZOS];
#endif // end FIND_MAG
/* LowestLanczos(k) constructs and diagonalizes the Hamilton operator for
k[] and given h/m and finds the ground state.
The routine returns the energy of the ground state.  */
double LowestLanczos(long long k[NSYM], komplex *resvect, long long *Nener,
                     long long flag, struct FLAGS *input_flags)
{
  double scale, emin;
  long long j, Nvec;
  long long i;
  char filename[150];
  FILE *seedfile;

  first = vec1;
  second = vec2;
  third = vec3;

  if (flag == CROSS)
    // flag bliver ikke brugt i MakeSeedCross() /ABP
    // Kopierer szxygs ind i first
    MakeSeedCross(first, CROSS);
  else
    MakeSeed(first);

  scale = Normalize(first);
  if (scale == 0.0)
  {
    if (Nener)
      *Nener = 0;
    return LARGE_NUMBER;
  }

#ifdef TEST_FINDGROUND
  if (flag == CROSS)
  {
    LogMessageChar("CALCULATION OF CROSSSECTION:   ");
  }
  LogMessageChar("Normalized seed vector: \n");
  for (i = 0; i < Nunique; i++)
  {
    LogMessageCharDouble("( ", real(first[i]));
    LogMessageCharDouble(" + i ", imag(first[i]));
    LogMessageChar(") \n");
  }
#endif

  if (input_flags->find_eigenstate)
  {
    if (flag == RECONSTRUCT)
    {
      sprintf(filename, "seed%lld.bin", Nspins);
      errno = 0;
      seedfile = fopen(filename, "w");
      if (seedfile == NULL)
      {
        fatalerror("Cannot open seedfile (writing), sorry!", errno);
        exit(-1);
      }

#ifdef TEST_FINDGROUND
      LogMessageChar("Seedvector before save to file:\n");
#endif /* TEST_FINDGROUND */
      // write seedvector to file for later use in finding the groundstate
      for (i = 0; i < Nunique; i++)
      {
#ifdef TEST_FINDGROUND
        LogMessageCharDouble("( ", real(first[i]));
        LogMessageCharDouble(" + i ", imag(first[i]));
        LogMessageChar(") \n");
#endif /* TEST_FINDGROUND */
        errno = 0;
        if (fwrite(&(first[i]), sizeof(komplex), 1, seedfile) <= 0)
        {
          fatalerror("Error in writing seedfile in lanczos for unique", errno);
          exit(1);
        }
      }
      fclose(seedfile);
    }
  }

  /* first[][] is now set with a proper seed */
  for (i = 0; i < Nunique; i++)
  {
    second[i] = zero;
    third[i] = zero;
  }

  z = kmatrix(1, MAX_LANCZOS, 1, MAX_LANCZOS);

  Nvec = LanczosLoop(0, k, NULL, input_flags);

  /* Copy magnetisation to arrays */
#ifdef FIND_MAG
  if (flag == NORMAL)
  {

    for (i = 0; i < Nvec; i++)
    {
#ifdef TEST_LANCZOS_VECTOR
      LogMessageCharInt("Vector number ", i);
      LogMessageChar(": (");
#endif
      magnetisation[i] = 0;
      for (int j = 0; j < Nvec; j++)
      {
        magnetisation[i] += ((real(z[j + 1][i + 1]) * real(z[j + 1][i + 1])) + (imag(z[j + 1][i + 1]) * imag(z[j + 1][i + 1]))) * lanczmag[j];
#ifdef TEST_LANCZOS_VECTOR
        LogMessageCharDouble("(", real(z[j + 1][i + 1]));
        LogMessageCharDouble(" i ) ,", imag(z[j + 1][i + 1]));
#endif
      }
#ifdef TEST_LANCZOS_VECTOR
      LogMessageChar(")\n");
#endif
    }
  }
#endif /* FIND_MAG */

  /* Copy eigenvalues and cross-sections to arrays */
  for (i = 0; i < Nvec; i++)
    energies[i] = diag_copy[i];

  if (input_flags->find_cross)
  {
    if (flag == CROSS)
    {
      for (i = 0; i < r; i++)
      {
        cross[i] = 2 * PI * sqrabs(scale * z[1][i + 1]); // TODO: check if this works!
#ifdef TEST_LANCCROSS
        LogMessageChar("k = [ ");
        for (j = 0; j < Nsym; j++)
          LogMessageInt(k[i]);
        LogMessageChar("] \n");
        LogMessageCharInt("Overlap from smpzq to exited state ", i);
        LogMessageCharDouble(" with energy ", energies[i]);
        LogMessageCharDouble(" is z =", real(z[1][i + 1]));
        LogMessageCharDouble(" + i ", imag(z[1][i + 1]));
        LogMessageCharInt("\n r =", r);
        LogMessageCharDouble(", scale =", scale);
        LogMessageCharDouble(" Cross section: ", cross[i]);
        LogMessageCharInt(", Nvec =", Nvec);
        LogMessageChar("\n");
#endif // TEST_LANCCROSS
      }
    }
  }

  /* Find index of smallest eigenvalue */
  emin = LARGE_NUMBER;
  for (i = 0; i < r; i++)
  {
#ifdef TEST_ENERGIES
    WriteEnergy(energies[i]);
#endif /* TEST_ENERGIES */
    if (energies[i] < emin)
    {
      emin = energies[i];
      imin = i;
    }
  }

#ifdef LANCZOS_MESSAGES
  LogMessageCharDouble(" Smallest eigenvalue:", energies[imin]);
  LogMessageCharInt(" no ", imin);
  LogMessageChar("\n");
#endif

  *Nener = r;

  if (input_flags->find_eigenstate)
  {
    if (flag == RECONSTRUCT)
    {
      /* Reconstruct Lanczos sequence to construct eigenvector */

      errno = 0;
      seedfile = fopen(filename, "r");
      if (seedfile == NULL)
      {
        fatalerror("Cannot open seed file (reading), sorry!", errno);
      }

      // read seedvector from seed file
#ifdef TEST_FINDGROUND
      LogMessageChar("Seed vector after load from file:\n");
#endif /* TEST_FINDGROUND */
      for (i = 0; i < Nunique; i++)
      {
        if (fread(&(first[i]), sizeof(komplex), 1, seedfile) <= 0)
        {
          LogMessageCharInt("Error in reading seedfile in lanczos for unique ", i);
          fatalerror("", errno);
        }

#ifdef TEST_FINDGROUND
        LogMessageCharDouble("", real(first[i]));
        LogMessageCharDouble(" + i ", imag(first[i]));
#endif /* TEST_FINDGROUND */
      }
#ifdef TEST_FINDGROUND
      LogMessageChar("Seed file read (while finding eigenvectors) \n");
#endif /* TEST_FINDGROUND */

      fclose(seedfile);

      for (i = 0; i < Nunique; i++)
        resvect[i] = zero;

#ifdef TEST_FINDGROUND
      LogMessageChar("Finding eigenvector\n");
#endif
      LanczosLoop(Nvec, k, resvect, input_flags);

#ifdef FIND_MAG
      for (i = 0; i < Nvec; i++)
      {
#ifdef TEST_FINDGROUND
        LogMessageCharInt("AFTER GS Vector number", i);
        LogMessageChar("\n");
#endif
        magnetisation[i] = 0;
        for (int j = 0; j < Nvec; j++)
        {
          magnetisation[i] += ((real(z[j + 1][i + 1]) * real(z[j + 1][i + 1])) + (imag(z[j + 1][i + 1]) * imag(z[j + 1][i + 1]))) * lanczmag[j];
#ifdef TEST_FINDGROUND
          LogMessageCharDouble("", real(z[j + 1][i + 1]));
          LogMessageCharDouble("+ i ", imag(z[j + 1][i + 1]));
          LogMessageCharDouble(" with ", lanczmag[j]);
          LogMessageChar("\n");
#endif
        }
#ifdef TEST_FINDGROUND
        LogMessageCharDouble("REAL MAG IS ", magnetisation[i]);
        LogMessageChar("\n");
#endif
      }
#endif // FIND_MAG

#ifdef TEST_EIG
      //  eigenvector_test(k,resvect,third);
#endif /* TEST_EIG */
    }
  }

  freekmatrix(z, 1, MAX_LANCZOS, 1, MAX_LANCZOS);

  return emin;
}

// LanczosLoop() makes one loop of the Lanczos algorithm. It
// returns the number of vectors generated. If called with
// Nvec > 0, it reconstructs the sequence without checking for
// convergence.
long long LanczosLoop(long long Nvec, long long k[NSYM], komplex *eigvect,
                      struct FLAGS *input_flags)
{
  long long j;
  long long i;
  double emin;
  double lastit = 10000;
  komplex product;
  long long iterations = 0;

  length = 1;
  /* Construct Lanczos sequence */
  r = 0;
#ifdef FIND_MAG
  lanczmag[0] = findmag(first);
#endif /* FIND_MAG */
#ifdef VERBOSE_LANCZOS
  LogMessageChar("Starting LanczosLoop()\n");
#endif
  while (length > ZERO_VEC_LENGTH)
  {
#ifdef VERBOSE_LANCZOS
    LogMessageCharInt("Nvec: ", Nvec);
    LogMessageCharInt("r: ", r);
    LogMessageCharDouble("length: ", length);
    LogMessageCharDouble("ZERO_VEC_LENGTH: ", ZERO_VEC_LENGTH);
    LogMessageChar("\n");
#endif
#ifdef TEST_FINDGROUND
    LogMessageCharDouble("Termination length = ", length);
    LogMessageChar("\n Lanczos vector: (");
    for (i = 0; i < Nunique; i++)
    {
      LogMessageCharDouble(",", real(first[i]));
      LogMessageCharDouble(" + i ", imag(first[i]));
    }
    LogMessageChar(")\n");
#endif // TEST_FINDGROUND
    // If reconstructing, reconstruct
    if (Nvec > 0)
    {
#ifdef TEST_FINDGROUND
      LogMessageChar("Before update: (");
      for (i = 0; i < Nunique; i++)
      {
        LogMessageCharDouble(",", real(eigvect[i]));
        LogMessageCharDouble(" + i ", imag(eigvect[i]));
      }
      LogMessageChar(")\n");
#endif // TEST_FINDGROUND

      for (i = 0; i < Nunique; i++)
      {
        eigvect[i] += z[r + 1][imin + 1] * first[i];
      }
#ifdef TEST_FINDGROUND
      LogMessageChar("After update:\n");
      for (i = 0; i < Nunique; i++)
      {
        LogMessageCharDouble(",", real(eigvect[i]));
        LogMessageCharDouble(" + i ", imag(eigvect[i]));
      }
      LogMessageChar(")\n");
#endif // TEST_FINDGROUND
    }

    length = NextLanczos(first, second, third, r, k, input_flags);
#ifdef TEST_FINDGROUND
    LogMessageCharInt("r = ", r);
    LogMessageCharDouble(", length = ", length);
#endif
    /* Swap vectors */
    temp = third;
    third = second;
    second = first;
    first = temp;

    r++;
#ifdef FIND_MAG
    lanczmag[r] = findmag(first);
#endif /* not M_SYM */

    // If not reconstructing, check for convergence
    if (Nvec == 0)
    {
      for (i = 1; i <= MAX_LANCZOS; i++)
        for (j = 1; j <= MAX_LANCZOS; j++)
          z[j][i] = zero;
      for (i = 1; i <= MAX_LANCZOS; i++)
        z[i][i] = one;
      for (i = 0; i < MAX_LANCZOS; i++)
        diag_copy[i] = diag[i];
      for (i = 0; i < MAX_LANCZOS; i++)
        subdiag_copy[i] = subdiag[i];

#ifdef TEST_FINDGROUND
      LogMessageCharDouble("diag[0]=", diag[0]);
      for (i = 1; i < r; i++)
      {
        LogMessageCharInt("i=", i);
        LogMessageCharDouble(", diag[i]=", diag[i]);
        LogMessageCharDouble(", subdiag[i]=", subdiag[i]);
        LogMessageChar("\n");
      }
#endif
      /*      LogMessageChar("Now starting htqli. "); */
      iterations = htqli(diag_copy - 1, subdiag_copy - 1, r, z);
      /*      LogMessageCharInt("Now at lanczos cycle ",r);
       * 	LogMessageCharInt(", used",iterations);
       *	LogMessageChar(" iterations \n");
       */

#ifdef TEST_EIGENVECTORS
      LogMessageCharInt("MAX_LANCZOS =", MAX_LANCZOS);
      LogMessageChar("Diagonalization done, resulting eigenvectors\n");
      for (i = 1; i <= MAX_LANCZOS; i++)
      {
        for (j = 1; j <= MAX_LANCZOS; j++)
        {
          LogMessageCharDouble(" ( ", real(z[i][j]));
          LogMessageCharDouble(" + i ", imag(z[i][j]));
          LogMessageChar(" ) ");
        }
        LogMessageChar("\n");
      }
#endif /* TEST_EIGENVECTORS */

#ifdef TEST_FINDGROUND
      for (i = 0; i < r; i++)
      {
        LogMessageCharInt("i = ", i);
        LogMessageCharDouble("diag_copy[i]=", diag_copy[i]);
        LogMessageChar("\n");
      }
#endif

      // Find index of smallest eigenvalue
      emin = LARGE_NUMBER;
      for (i = 0; i < r; i++)
      {
        if (diag_copy[i] < emin)
        {
          emin = diag_copy[i];
          imin = i;
        }
      }
      // Check magnitude of last element in Ritz vector
      // TODO: This criteria is very crude - an improvement could be made
      // using e.g. the backward error.
      if (subdiag[r] * abs(z[r][imin + 1]) < Ritz_conv)
      {
        length = 0;
        /*          LogMessageCharDouble("subdiag is too low: ",subdiag[r]);
                  LogMessageCharDouble("Ritz_conv is : ",Ritz_conv);
            LogMessageChar("\n");
        */
      }
      if (lastit == 0.0)
      {
        length = 0;
      }
      else
      {
        if (lastit < Ritz_conv)
        {
#ifdef LANCZOS_MESSAGES
          LogMessageCharDouble("STOPPING. Criterion ", subdiag[r] * abs(z[r][imin + 1]));
          LogMessageCharInt(", step ", r);
          LogMessageChar("\n");
#endif
          length = 0;
        }
        else
        {
          lastit = subdiag[r] * abs(z[r][imin + 1]);
        }
      }
    }
    // If specified number of vectors have been generated, terminate
    if (r == Nvec)
      length = 0;
    if ((Nvec == 0) && ((r % 1) == 0))
    {
#ifdef LANCZOS_MESSAGES
      LogMessageCharDouble(" Termination length is ", subdiag[r] * abs(z[r][imin + 1]));
      LogMessageCharDouble(" and lastit is ", lastit);
      LogMessageCharInt(" step ", r);
      LogMessageChar("\n");
#endif
    }
    if (r > MAX_LANCZOS - 1)
    {
      LogMessageCharDouble(" MAXLANCZ reached. Terminating at length ", length);
      length = 0;
    }
  }

#ifdef LANCZOS_MESSAGES
  LogMessageCharDouble("Termination length at termination = ", length);
  LogMessageCharInt("\n Lanzcos sequence terminated with r = ", r);
  LogMessageChar("\n");
#endif
#ifdef VERBOSE_LANCZOS
  LogMessageChar("Ending LanczosLoop()\n");
#endif

  return (r);
}

/* NextLanczos() constructs the next vector in the Lanczos sequence. It
returns the length of the vector before normalisation. The next vector
is returned in dummy; vectors must be swapped after calling.
vr = first, vr_1 = second, dummy = third */
double NextLanczos(komplex *vr, komplex *vr_1,
                   komplex *dummy,
                   unsigned long long r, long long k[NSYM],
                   struct FLAGS *input_flags)
{
  long long j;
  long long i;
  double lg, tmp;

  /* Construct |ur> = H|vr> - |vr_1><vr_1|H|vr> */
#ifdef TEST_APPLYSPARSE
  LogMessageCharInt("Starting Applysparse now, step ", r);
  LogMessageChar("\n");
#endif
  ApplySparse(vr, dummy, k, input_flags); // make dummy = H |vr>
#ifdef TEST_APPLYSPARSE
  LogMessageChar("Stopping Applysparse now\n");
#endif

  if (r != 0)
    for (i = 0; i < Nunique; i++)
      dummy[i] -= subdiag[r] * vr_1[i];

  /* Find <vr|H|vr> */
  lg = 0;
  for (i = 0; i < Nunique; i++)
    lg += real(vr[i]) * real(dummy[i]) + imag(vr[i]) * imag(dummy[i]);
  diag[r] = lg;

  /* Construct |wr> = |ur> - |vr><vr|H|vr> */
  for (i = 0; i < Nunique; i++)
    dummy[i] -= lg * vr[i];

  /* Set |v(r+1)> = |wr> / || |wr> || */
  lg = Normalize(dummy);

  /* Set subdiagonal element */
  if ((r + 1) < Nunique)
    subdiag[r + 1] = lg;

#ifdef DEBUG_LANCSTEP
  // output the vector and the subdiagonal element
  LogMessageCharDouble("next subdiagonal element (beta): ", lg);
  LogMessageChar("\n Next vector in trigonalization series: (");
  for (long long i = 0; i < Nunique; i++)
  {
    LogMessageCharDouble(" ", real(dummy[i]));
    LogMessageCharDouble(" + i ", imag(dummy[i]));
  }
  LogMessageCharDouble(") Length of this vector: ", lg);
  LogMessageChar("\n");
#endif /* DEBUG_LANCSTEP */

  return (lg);
}

/* MakeSeed() computes the seed to be used in the Lanzcos algorithm */
void MakeSeed(komplex *seed)
{
  unsigned long long i;

#ifdef TEST_SEED
  LogMessageChar("MakeSeed called \n");
#endif

  /* Make seed state with low energy */
  // FOR 1D : lowest Ising state |010101...>, has highest unique index -->  seed[0][Nunique-1] = 1; */
#ifdef DEBUG_SEED
  for (i = 0; i < Nunique; i++)
  {
    if (Nocc[i] > 0)
      seed[i] = zero;
    else
      seed[i] = zero;
  }
  seed[2] = 1.0 + I * 0.0;
  OutMessageChar(" DEBUG_SEED on. \n");

#else  /*not DEBUG_SEED */
  for (i = 0; i < Nunique; i++)
  {
    if (Nocc[i] > 0)
    {
      seed[i] = 2.0 * ((double)rand() / RAND_MAX) - 1.0 + I * 2.0 * ((double)rand() / RAND_MAX) - 1.0;
    }
    else
      seed[i] = zero;
  }
#endif /* DEBUG_SEED */

#ifdef TEST_SEED
  LogMessageChar("MakeSeed returning: (");
  for (i = 0; i < Nunique; i++)
  {
    LogMessageCharDouble(", ", real(seed[i]));
    LogMessageCharDouble(" + i ", imag(seed[i]));
  }
  LogMessageChar(") \n");
#endif
  return;
}

/* MakeSeedCross() creates the seed for the Lanczos algorithm used
   for cross section calculations */
void MakeSeedCross(komplex *first, long long flag)
{
  unsigned long long i;
  long long j;

  /* Use  Smp(k)|g>, Spm(k), or Szz(k)|g> as seed */
  for (i = 0; i < Nunique; i++)
    first[i] = szxygs[i];

#ifdef TEST_LANCCROSS
  for (int i = 0; i < Nunique; i++)
  {
    LogMessageCharDouble("\nfirst[i] =", real(first[i]));
    LogMessageCharDouble("+ i", imag(first[i]));
  }
#endif

  return;
}

// findmag finds the magnetisation of a given vector
#ifdef FIND_MAG
double findmag(komplex *vector)
{
  double result = 0;

#ifdef TEST_FINDMAG
  LogMessageChar("TESTING FINDMAG:\n");
#endif
  for (int n = 0; n < Nunique; n++)
  {
    result += ((real(vector[n]) * real(vector[n])) + (imag(vector[n]) * imag(vector[n]))) * mag[n];
    /*TODO: Magnetisation must be calculated directly from the unique-vectors and not, as now, a weighed average of the lanczos vectors */
#ifdef TEST_FINDMAG
    LogMessageCharDouble(" ", real(vector[n]));
    LogMessageCharDouble(" + i ", imag(vector[n]));
#endif
  }
#ifdef TEST_FINDMAG
  LogMessageCharDouble("GIVES: ", result);
  LogMessageChar("\n");
#endif
  return result;
}

void findmaggs()
{
  int i;
  double numSpins = Nspins;
  maggs = 0;
  for (i = 0; i < Nunique; i++)
  { // Nunique is number of uniques
    // LogMessageCharInt("Mag = ", (mag[i]));
    maggs += 1 / numSpins * (real(gs[i]) * real(gs[i]) + imag(gs[i]) * imag(gs[i])) * (mag[i]);
  }
}
#endif /* FIND_MAG */

// Normalize() normalizes a komplex vector, returning the length before normalization
double Normalize(komplex *vector1)
{
  long long i;
  double lg;

  lg = 0;

  for (i = 0; i < Nunique; i++)
    lg += SQR(real(vector1[i])) + SQR(imag(vector1[i]));

  lg = sqrt(lg);

  if (lg != 0)
    for (i = 0; i < Nunique; i++)
      vector1[i] /= lg;

  return (lg);
}