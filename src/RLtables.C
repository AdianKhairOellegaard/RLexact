/* Program file RLtables.C -
* Filling and using tables for lookup and symmetry management
* Last change: SJ 06.01.17
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
#include <complex>

#include <RLexact.h>
#include <nr.h>
#include <cnr.h>
#include <errno.h>
#include <string.h>
#include "Functions.h"

/* Functions defined in this file */


/* Functions defined elsewhere */
extern unsigned long long SymOp(long long, unsigned long long);


/* Variables defined in RLexact.c */
extern long long Nspins, Nsym, Nsymvalue[], symlist[];
extern long long Nu2, Nunique, Nuniq_k;
extern double cosine[], sine[], sqroot[];
extern unsigned long long *unique;
extern long long *uniq_k;
extern long long *Nocc, *Nocc_0;
extern char *infile_name;
extern long long *mag;
extern long long twom;

/* InvertMatrix () inverts the translation matrix. TODO: replace by the other matrix routine! */
/* The rational behind this is the following:
   The eigenvalues for the translation operator T_d is t_d = exp(i j_d d.q),
   where d is a translation vector (typically nearest neighbour),
   j_d is an integer characterizing the state, and q is the "real space"
   propagation vector of the state.
   As we have periodic boundary conditions, we must have
   t_d = exp(2 pi / n_d), where n_d is the periodicity. So d.q=2 pi j_d/n_d .
   There are D translation vectors, where D is the dimensionality
   of the problem. Writing the D eigenvalue equations in matrix form,
   we reach  (d).q = J  , where the elements in J are J_d=j_d/n_d,
   the matrix (d) is expressed in units of the lattice constant, a, and
   q is written in units of (2 pi/a). What we are calculating here is (d)^-1 */
void InvertMatrix(long long dim, double **matrix, double **inverse, struct FLAGS *input_flags)
{
  double det, d1, d2;
  long long i, j;

  inverse = (double **)malloc(4 * sizeof(double *));
  for (long long k = 0; k < 4; k++)
  {
    inverse[k] = (double *)malloc(4 * sizeof(double));
  }

  LogMessageCharInt("Dimension of the problem: ", dim);
  if (dim == 1)
  {
    inverse[0][0] = 1 / matrix[0][0];
  }
  if (dim == 2)
  {

    det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
    d1 = matrix[0][0] * matrix[1][1];
    d2 = matrix[1][0] * matrix[0][1];
    inverse[0][0] = matrix[1][1] / det;
    inverse[1][1] = matrix[0][0] / det;
    inverse[0][1] = -matrix[1][0] / det;
    inverse[1][0] = -matrix[0][1] / det;
  }

if (input_flags->TEST_INVERTMATRIX)
{
  LogMessageCharInt("InvertMatrix: dim = ", dim);
  LogMessageChar("\n InvertMatrix called with matrix:\n");
  for (i = 0; i < dim; i++)
  {
    LogMessageChar("(");
    for (j = 0; j < dim; j++)
      LogMessageCharDouble(" ", matrix[i][j]);
    LogMessageChar(")\n");
  }

  LogMessageChar("\n InvertMatrix returning matrix:\n");
  for (i = 0; i < dim; i++)
  {
    LogMessageChar("(");
    for (j = 0; j < dim; j++)
      LogMessageCharDouble(" ", inverse[i][j]);
    LogMessageChar(")\n");
  }
}

  if (dim > 2)
    fatalerror("In InvertMatrix: No of dimensions too large", dim);
  return;
}

/* BuildTables() makes tables of cos, sin, sqrt, and
   phases for structure factors */
void BuildTables(struct FLAGS *input_flags)
{
  long long i, j, q;
  double tmp;

  for (i = 0; i < Nspins; i++)
  {
    cosine[i] = cos(2 * PI * i / Nspins);
    // Beware of period lengths on lookup !!
    sine[i] = sin(2 * PI * i / Nspins);
  }
  if (Nspins % 2 == 0)      // If Nspins is even
    sine[Nspins / 2] = 0.0; // Avoid floating-point errors
  for (i = 0; i < 2 * (Nspins + 1) * Nspins; i++)
    sqroot[i] = sqrt((double)i);

if (input_flags->TEST_TABLES)
{
  for (i = 0; i < Nspins; i++)
    LogMessageCharInt(" i:", i);
  LogMessageChar3Vector(". (cos, sin, sqrt): ", cosine[i], sine[i], sqroot[i]);
}
}

/* Count() counts the number of up spins in a state */
/* TODO: modify for S>1/2 */
long long Count(unsigned long long state, struct FLAGS *input_flags)
{
  long long i, c = 0;

  for (i = 0; i < Nspins; i++)
  {
    c += ((state & (((unsigned long long)1) << i)) != 0);
if (input_flags->TEST_COUNT)
{
    LogMessageCharInt(" State: ", state);
    LogMessageCharInt("i = ", i);
    LogMessageCharInt(" count: ", c);
    LogMessageChar("\n");
}
  }

  return c;
}

/* ReadUniqueObservables() reads diagonal elements for all uniques from a file */
/* */
void ReadUniqueObservables(struct FLAGS *input_flags)
{
  char uniqueobs_name[30];
  FILE *uniobsfile;
  long long i;
  char m_name[2];

  // first, we construct the appropriate name of the unique file
  if (input_flags->m_sym)
  {                     // At present ReadUniqueObservables are not called by RLexact for M_SYM on, so this #ifdef does not make much sense.
    itoa(twom, m_name); // TT debugging 5/2-10
    strcpy(uniqueobs_name, infile_name);
    strcat(uniqueobs_name, ".m"); // The line strcat(uniqueobs_name,UNIOBSEND); was originally here, but it must have been a bug never discovered due to the function never being called in this #ifdef.
    strcat(uniqueobs_name, m_name);
    strcat(uniqueobs_name, UNIOBSEND);
  }
  else
  { /*M_SYM */
    strcpy(uniqueobs_name, infile_name);
    strcat(uniqueobs_name, UNIOBSEND);
  }

  // then, we open the file
  errno = 0;
  uniobsfile = fopen(uniqueobs_name, "r");
  if (uniobsfile == NULL)
  {
    printf("Expected uniqueobs_name %s \n", uniqueobs_name);
    printf("Expected infile_name %s \n", infile_name);
    fatalerror("Cannot open unique observables file, sorry", errno);
  }
  fflush(uniobsfile);

  if (!input_flags->m_sym)
  {
    // read the magnetisation and periodicity in q=0
    errno = 0;
    if (Nunique != fread(mag, sizeof(long long), Nunique, uniobsfile))
    {
      fatalerror("Cannot read mag array from uniques observables file, sorry!", errno);
    }

    errno = 0;
    if (Nunique != fread(Nocc_0, sizeof(long long), Nunique, uniobsfile))
    {
      fatalerror("Cannot read Nocc_0 array from uniques observables file, sorry!", errno);
    }
  }

  fclose(uniobsfile);
}

/* ReadUnique() reads the number of uniques and a list of the uniques from file */
long long ReadUnique(long long twom, int Nunique_only, struct FLAGS *input_flags)
{
  char unique_name[30], m_name[2];
  FILE *unifile;
  long long i;

  // first, we construct the appropriate name of the unique file
  strcpy(unique_name, infile_name); // TODO: generalize for other file names
  if (input_flags->m_sym)
  {
    /*  LogMessageCharInt("\n 2m value for the un file is:",twom);
      itoa(twom,m_name);
      strcat(unique_name,".m");
      strcat(unique_name,m_name); */
  }
  strcat(unique_name, UNIEND);

  // then, we open the file
  errno = 0;
  unifile = fopen(unique_name, "r");
  if (unifile == NULL)
  {
    LogMessageChar("\n In ReadUnique: Searching for Unique file ");
    LogMessageChar(unique_name);
    fatalerror("\n In ReadUnique: Cannot open uniques file", errno);
  }
  fflush(unifile);

  // now, read the number of uniques, Nu2, and the unique[] table
  if (1 != fread(&Nunique, sizeof(long long), 1, unifile))
    fatalerror("Cannot read Nunique from uniques file", errno);
  if (!Nunique_only) // Read also the other contents in file
  {
    if (1 != fread(&Nu2, sizeof(long long), 1, unifile))
      fatalerror("Cannot read Nu2 from uniques file", errno);
    if (Nunique != fread(unique, sizeof(unsigned long long), Nunique, unifile))
      fatalerror("Cannot read unique array from uniques file", errno);
  }

if (input_flags->TEST_READUNIQUE)
{
  LogMessageCharInt("Nunique found: ", Nunique);
  for (int i = 0; i < Nunique; i++)
  {
    LogMessageCharInt("unique no. ", i);
    if (!Nunique_only)
    {
      LogMessageCharInt("is ", unique[i]);
    }
    LogMessageChar("\n");
  }

}

  fclose(unifile);
  return (Nunique);
}

/* WriteUnique writes the number of uniques, followed by the unique table to the unique-file. */
void WriteUnique(long long twom, struct FLAGS *input_flags)
{
  char unique_name[30];
  char m_name[2];
  FILE *unifile;
  long long i;

  // First, construct the appropriate name of the unique file
  strcpy(unique_name, infile_name);
  if (input_flags->m_sym)
  {
    /* itoa(twom,m_name);
    strcat(unique_name,".m");
    strcat(unique_name,m_name);
    LogMessageChar("m_name : ");
    LogMessageChar(m_name); */
  }
  strcat(unique_name, UNIEND);
  LogMessageChar("unique_name : ");
  LogMessageChar(unique_name);

  // Then, open the file
  errno = 0;
  unifile = fopen(unique_name, "w");
  if (unifile == NULL)
  {
    LogMessageChar("Filename was:");
    LogMessageChar(unique_name);
    fatalerror("Cannot open uniques file", errno);
  }
  fflush(unifile);

  // Now, write in the number of uniques (beginning of file)
  errno = 0;
  if (1 != fwrite(&Nunique, sizeof(long long), 1, unifile))
    fatalerror("Cannot write Nunique to uniques file", errno);

  // then Nu2
  if (1 != fwrite(&Nu2, sizeof(long long), 1, unifile))
    fatalerror("Cannot write Nu2 to uniques file", errno);

  // and then the unique[] table
  if (Nunique != fwrite(unique, sizeof(unsigned long long), Nunique, unifile))
    fatalerror("Cannot write unique array to uniques file", errno);

  fclose(unifile);
}

/* WriteUniqueObservables() writes diagonal elements of the uniques to file */
void WriteUniqueObservables(struct FLAGS *input_flags)
{
  char uniqueobs_name[30];
  FILE *uniobsfile;
  long long i;

  // first, we construct the appropriate name of the unique file
  strcpy(uniqueobs_name, infile_name);
  strcat(uniqueobs_name, UNIOBSEND);
  LogMessageChar("filename:");
  LogMessageChar(uniqueobs_name);

  // then, we open the file
  errno = 0;
  uniobsfile = fopen(uniqueobs_name, "w");
  if (uniobsfile == NULL)
  {
    fatalerror("Cannot open unique observables file, sorry!", errno);
  }
  fflush(uniobsfile);

  if (!input_flags->m_sym)
  {
    // write the magnetisation
    errno = 0;
    if (Nunique != fwrite(mag, sizeof(long long), Nunique, uniobsfile))
    {
      fatalerror("Cannot write mag array to uniques observables file, sorry!", errno);
    }
    // write the periodicity in q=0 for use in cross sections
    if (Nunique != fwrite(Nocc_0, sizeof(long long), Nunique, uniobsfile))
    {
      fatalerror("Cannot write Nocc_0 array to uniques observables file, sorry!", errno);
    }
  }
  fclose(uniobsfile);
}

/* FillUnique fills the table of unique states if the flag CountOnly is not set.
The number of unique states are returned in any case. */
unsigned long long FillUnique(long long twom, int CountOnly, struct FLAGS *input_flags)
{
  unsigned long long bitmap, basis_c = 0;
  long long i, j, n, pos[NSPINS + 1];

if (input_flags->TEST_FILLUNIQUE)
{
  if (CountOnly == 1)
    LogMessageChar("\nCountOnly\n");
  LogMessageCharInt("Start of FillUnique. MAX_STATE = ", (unsigned long long)MAX_STATE);
  LogMessageChar("\n");
}

  if (input_flags->m_sym)
  {
    n = (Nspins + twom) / 2; // number of spin-ups
                             // Begin with all spins to the right
    for (j = 0; j < n; j++)
      pos[j] = j;
    bitmap = 0;
    for (j = 0; j < n; j++)
      bitmap += ((unsigned long long)1) << pos[j];
if (input_flags->TEST_FILLUNIQUE)
{
    LogMessageCharInt(" Loop over all states with ", n);
    LogMessageCharInt(" spins up starting with state ", bitmap);
    LogMessageChar("\n");
}
    while (bitmap < ((unsigned long long)MAX_STATE))
    {
if (input_flags->TEST_FILLUNIQUE)
{
      LogMessageCharInt("\nLooking at state ", bitmap);
}
      if (IsUnique(bitmap, input_flags))
      {
if (input_flags->TEST_FILLUNIQUE)
{
        LogMessageCharInt("\nEnter IsUnique criteria, basis_c=", basis_c);
}
        if (!CountOnly)
        {
          unique[basis_c] = (unsigned long long)bitmap;
        }
        basis_c++;
if (input_flags->TEST_FILLUNIQUE)
{
        LogMessageCharInt("\n basis_c=", basis_c);
        if (!CountOnly)
        {
          LogMessageCharInt("Its Unique, number ", basis_c - 1);
          LogMessageCharInt(": ", unique[basis_c - 1]);
          LogMessageCharInt("is ", bitmap);
          LogMessageChar("\n");
        }
        else
        {
          LogMessageCharInt("Its Unique, number ", basis_c - 1);
          LogMessageCharInt("is ", bitmap);
          LogMessageChar("\n");
        }
}
      }
      j = 0; // now, reposition the spins to cover all combinations of spin up
      while ((j < n - 1) && (pos[j + 1] == (pos[j] + 1)))
        j++;
      bitmap += (((unsigned long long)1) << (pos[j] + 1)) - (((unsigned long long)1) << pos[j]);
      pos[j]++;
      for (i = 0; i < j; i++)
      {
        bitmap += (((unsigned long long)1) << i) - (((unsigned long long)1) << pos[i]);
        pos[i] = i;
      }
    }
  }
  else
  { // M_SYM

    if (Nspins >= 48)
      fatalerror("In FillUnique: Cannot run through all states for Nspins>=", 48);
    if (Nspins > 32)
    {
      Warning("You are looping over many states in FillUnique:", exp(Nspins));
    }
    for (bitmap = 0; bitmap < MAX_STATE; bitmap++)
    {
if (input_flags->TEST_FILLUNIQUE)
{
      LogMessageCharInt(" basis_c = ", basis_c);
      LogMessageCharInt(" bitmap = ", bitmap);
}
      if (IsUnique(bitmap, input_flags))
      { // Unique found
        if (!CountOnly)
        {
          unique[basis_c] = bitmap;
if (input_flags->TEST_FILLUNIQUE_LIST)
{
          LogMessageCharInt("\n Unique no: ", basis_c);
          LogMessageCharInt(" , ", unique[basis_c]);
}
        }
        basis_c++;
      }
    }
  }

  Nu2 = ((unsigned long long)1) << (long long)ceil(log(basis_c + 1) / log(2)); // find smallest Nu2 = 2^^j > n

if (input_flags->TEST_FILLUNIQUE)
{
  LogMessageCharInt("\n Number of uniques found: ", basis_c);
  LogMessageChar("End of FillUnique  reached. \n");
}

  return (basis_c);
}

/* FillUniqueObservables fills the table of observables of unique states. */
void FillUniqueObservables(struct FLAGS *input_flags)
{
  LogMessageChar("IN FillUniqueObservables\n");
  unsigned long long state, index;
  long long count, q, j;
  double spin_j;
  komplex sum;

  for (index = 0; index < Nunique; index++)
  {
    state = unique[index];
    count = Count(state, input_flags);
if (input_flags->TEST_FILLUNIQUEOBSERVABLES)
{
    LogMessageCharInt("State no ", index);
    LogMessageCharInt("out of states numbers ", Nunique);
    LogMessageCharInt(" : ", state);
    LogMessageCharInt(", count: ", count);
}
    if (!input_flags->m_sym)
    {
      mag[index] = count - (Nspins + 1) / 2;
      ;
      /* Integer division, force integer m, even if Sz is half-integer */
if (input_flags->TEST_FILLUNIQUEOBSERVABLES)
{
      LogMessageCharInt(", m=", mag[index]);
}
    }

if (input_flags->TEST_FILLUNIQUEOBSERVABLES)
{
    LogMessageChar(", Done! \n");
}
  }

  return;
}

/* Build table of occurencies in a full cycle, for one q-value */
void BuildCycle(long long q[NSYM], struct FLAGS *input_flags)
{
  long long sym, s, j, count, tot_count, phi, uniq_count = 0;
  int T[NSYM];
  unsigned long long state;
  unsigned long long new_state;
  double p_i, p_r;
  double qlength;

if (input_flags->TEST_OCC)
{
  LogMessageChar("\n BuildCycle called for");
  LogMessageCharInt(" q =(", q[0]);
  LogMessageCharInt(",", q[1]);
  // LogMessageCharInt("), 2*m=",twom);
  LogMessageChar("\n");
}

  for (j = 0; j < Nunique; j++)
  {
    state = unique[j];
    count = 0;
    tot_count = 0;
    p_r = 0;
    p_i = 0;

    /* Run through symmetry operations on bitmap and count how
    many times it is encountered */
    TLOOP_BEGIN
    tot_count++;
    if (new_state == state)
    {
      count++;
      /* Calculate the phase, TODO: ADD MORE HERE ??? */
      for (phi = 0, s = 0; s < Nsym; s++)
      {
        phi += T[s] * q[s] * Nspins / Nsymvalue[s];
      }
      phi = phi % Nspins;
      p_r += cosine[phi];
      p_i += sine[phi];
    }
    TLOOP_END

    if (tot_count % count != 0) // re-occurence of a state in the symmetry loop should be a simple fraction of the total cycle length
    {
      LogMessageCharInt("Counts ", count);
      LogMessageCharInt(" of ", tot_count);
      fatalerror("In BuildCycle: Inconsistency in T-LOOP for state", unique[j]);
    }
    if (fabs(p_i) > SMALL_NUMBER) // summed phase factor over a cycle cannot be imaginary
    {
      LogMessageCharDouble("p: (", p_r);
      LogMessageCharDouble(" + i ", p_i);
      LogMessageChar(") \n");
      fatalerror("In BuildCycle: Imaginary value in BuildCycle for state", unique[j]);
    }
    Nocc[j] = (long long)floor(p_r + SMALL_NUMBER); // round off phase sum to a real integer

    // save the Nocc in q=0 for use in cross section calcs with Lanczos
    if (input_flags->use_lanczos)
    {
      qlength = 0;
      for (int k = 1; k < Nsym; k++)
      {
        qlength += q[k] * q[k];
      }
      if (qlength == 0)
      {
        Nocc_0[j] = Nocc[j];
      }
    }

    if (input_flags->use_exact_matrix)
    {
      if (Nocc[j] > 0)
        uniq_k[uniq_count++] = j;
    }

if (input_flags->TEST_OCC)
{
    LogMessageCharInt(" Unique", j);
    LogMessageCharInt(" : ", unique[j]);
    LogMessageCharInt(" , phi=", phi);
    LogMessageCharInt(" , Nocc[j]=", Nocc[j]);
    LogMessageCharInt(" , Nocc_0[j]=", Nocc_0[j]);
    LogMessageCharInt(" , uniq_count=", uniq_count);
    LogMessageChar("\n");
}
  }

  if (input_flags->use_exact_matrix)
  {
    Nuniq_k = uniq_count;
  }

if (input_flags->TEST_OCC)
{
  LogMessageChar(" Exit BuildCycle \n");
}
  return;
}

/* FindUnique() finds the unique state corresponding to
an arbitrary state by searching for the smallest number
found by applying symmetry operations. */

unsigned long long FindUnique(unsigned long long state, int *Tvec, struct FLAGS *input_flags)
{
  /* This routine is where most of the time is spent (for Lanczos). */

  long long sym, i;
  int T[NSYM];
  unsigned long long new_state, min;

if (input_flags->TEST_FINDUNIQUE)
{
  LogMessageCharInt(" FindUnique called with state ", state);
  LogMessageChar("\n");
}

  min = state;
  for (i = 0; i < Nsym; i++)
    Tvec[i] = 0;
  /* Run through symmetry operations on bitmap to find the
  lowest state encountered */
  TLOOP_BEGIN
  if (new_state < min)
  {
    min = new_state;
    for (i = 0; i < Nsym; i++)
      Tvec[i] = T[i];

if (input_flags->TEST_FINDUNIQUE_DETAIL)
{
    LogMessageCharInt("sym =", sym);
    LogMessageChar("\n T=(");
    for (sym = 0; sym < Nsym; sym++)
      LogMessageInt(T[sym]);
    LogMessageChar(") Tvec=(");
    for (sym = 0; sym < Nsym; sym++)
      LogMessageInt(Tvec[sym]);
    LogMessageCharInt(") State: ", new_state);
    LogMessageChar("\n");
}
  }
  TLOOP_END

if (input_flags->TEST_FINDUNIQUE)
{
  LogMessageCharInt(" Found u = ", min);
  LogMessageChar("\n");
}

  return min;
}

/* Check if a state is a unique state
   Code similar to FindUnique */
long long IsUnique(unsigned long long state, struct FLAGS *input_flags)
{
  int j;

  long long sym;
  int T[NSYM];
  unsigned long long new_state;

  /* Run through symmetry operations on state */
  TLOOP_BEGIN
if (input_flags->TEST_ISUNIQUE)
{
  LogMessageChar(" T =");
  for (j = 0; j < Nsym; j++)
    LogMessageInt(T[j]);
  LogMessageCharInt(" led to state ", new_state);
}
  if (new_state < state)
  {
    // LogMessageChar("\n");
    return 0; /* state is not unique */
  }
  TLOOP_END

if (input_flags->TEST_ISUNIQUE)
{
  LogMessageChar("Is unique\n");
}

  return 1;
}

/* LookUpU() looks up the index of a unique state
             in the sorted unique table. */
long long LookUpU(unsigned long long u, struct FLAGS *input_flags)
{
  unsigned long long entry;
  long long index = Nu2 >> 1, Nb2 = Nu2 >> 2;

  while (1)
  {
    entry = unique[index - 1];
if (input_flags->TEST_LOOKUP)
{
    LogMessageCharInt("Lookup searching for", u);
    LogMessageCharInt(", reached ", entry);
}
    if (entry > u)
      index -= Nb2;
    else if (entry < u)
    {
      index += Nb2;
      if (index > Nunique)
        index = Nunique;
    }
    else                  // entry = u
      return (index - 1); // Correct for offset of the table
    if (Nb2 == 0)
      fatalerror("In LookUpU: Unique not found:", u);
    Nb2 = Nb2 >> 1;
  }
}


