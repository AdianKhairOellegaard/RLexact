/* Program file RLsymm.C -
* Initializing and performing the symmetry operations
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex>
#include <RLexact.h>
#include <cnr.h>
#include "Functions.h"

/* Functions defined in this file */


/* Global variables defined in RLexact.c */
extern long long Nspins, Ncoup, Nsym, Nsymvalue[NSYM], Nsymadd;
extern long long hamil_coup[NCOUP][2], symlist[NSYM];
extern long long **symadd;
extern double Jxy[NCOUP], Jzz[NCOUP], Janis[NCOUP];
extern double Jr[NRING];
extern double Jdip[NCOUP], geom_13[NCOUP], r_vector[NCOUP][3];

/* Regional variables defined in this file */
unsigned long long flipmask;

/* Perform symmetry operations. Containing spinflip, linear, N=32 fcc.
   TODO: Should probably even be coded with explicit operations and no arrays.  */
unsigned long long SymOp(long long sym_index, unsigned long long bitmap)
{
  /* Reflection in plane */
  static long long M0[32] = {0, 19, 2, 17, 4, 23, 6, 21, 8, 27, 10, 25, 12, 31, 14, 29,
                             16, 3, 18, 1, 20, 7, 22, 5, 24, 11, 26, 9, 28, 15, 30, 13};
  /* Rotation around 3-fold axis */
  static long long R3[32] = {0, 7, 10, 13, 3, 6, 9, 12, 2, 5, 8, 15, 1, 4, 11, 14,
                             18, 21, 24, 31, 17, 20, 27, 30, 16, 23, 26, 29, 19, 22, 25, 28};
  /* Rotation around 4-fold axis */
  static long long R4[32] = {0, 13, 10, 7, 28, 25, 22, 19, 8, 5, 2, 15, 20, 17, 30, 27,
                             24, 21, 18, 31, 4, 1, 14, 11, 16, 29, 26, 23, 12, 9, 6, 3};
  /* Rotation around 2-fold axis */
  static long long R2[32] = {0, 7, 10, 13, 4, 11, 14, 1, 8, 15, 2, 5, 12, 3, 6, 9,
                             24, 31, 18, 21, 28, 19, 22, 25, 16, 23, 26, 29, 20, 27, 30, 17};

  unsigned long long new_state;
  long long j, s, sa;

  s = symlist[sym_index];
  switch (s)
  {
  case IDENTITY:
    return bitmap;
  case SPIN_FLIP:
    return (~bitmap) & flipmask;
  case LINEAR_T:
    new_state = right(bitmap);
    return new_state;
  case FCC32Tx: /* Translation along line */
                /* NB! NB! "j" changed by "new_state". Possible error in FCC32 calc. ??! */
    new_state = (bitmap & MASK_L2) >> 3;
    new_state += ((bitmap & MASK_L1) << 1);
    new_state += ((bitmap & MASK_L2) >> 3) + ((bitmap & MASK_L1) << 1);
    return new_state;
  case FCC32Ty: /* Translation between lines */
    return ((bitmap & MASK_L3) >> 12) + ((bitmap & MASK_L4) << 4);
  case FCC32Tz: /* Translation between planes */
    new_state = ((bitmap & MASK_P1) << 16) + ((bitmap & MASK_P2) >> 16);
    new_state = ((new_state & MASK_L2) >> 03) + ((new_state & MASK_L1) << 01);
    return ((new_state & MASK_L3) >> 12) + ((new_state & MASK_L4) << 4);
  case FCC32MIRROR: /* Mirror reflection */
    new_state = 0;
    for (j = 0; j < 32; j++)
      if (bitmap & (1 << j))
        new_state += 1 << M0[j];
    return new_state;
  case FCC32R2: /* Two-fold rotation */
    new_state = 0;
    for (j = 0; j < 32; j++)
      if (bitmap & (1 << j))
        new_state += 1 << R2[j];
    return (new_state);
  case FCC32R4: /* Four-fold rotation */
    new_state = 0;
    for (j = 0; j < 32; j++)
      if (bitmap & (1 << j))
        new_state += 1 << R4[j];
    return (new_state);
  case FCC32R3: /* Three-fold rotation */
    new_state = 0;
    for (j = 0; j < 32; j++)
      if (bitmap & (1 << j))
        new_state += 1 << R3[j];
    return (new_state);
  default: /* interpret added operators */
    new_state = 0;
    sa = s - NSYM;
    // Safely delete
    // LogMessageCharInt(" SymOp ",sa);
    // LogMessageCharInt(", state ",bitmap);
    //
    if ((sa >= 0) && sa < Nsymadd)
    {
      for (j = 0; j < Nspins; j++)
      {
        if (bitmap & (((unsigned long long)1) << j))
          new_state += ((unsigned long long)1) << symadd[sa][j];
      }
      // Safely delete
      // LogMessageCharInt("new state, ",new_state);
      // LogMessageChar("\n");
      //
      return (new_state);
    }
    else
      fatalerror(" Unknown symmetry operation ", s);
  }
  return 0; /* Never get here */
}

/* Symmetry operations for point symmetries on site indices.
 TODO: Should be updated to be parallel to SymOp */
unsigned long long SymOpSite(long long i, long long site)
{
  /* Reflection in plane */
  static long long M0[32] = {0, 19, 2, 17, 4, 23, 6, 21, 8, 27, 10, 25, 12, 31, 14, 29,
                             16, 3, 18, 1, 20, 7, 22, 5, 24, 11, 26, 9, 28, 15, 30, 13};
  /* Rotation around 3-fold axis */
  static long long R3[32] = {0, 7, 10, 13, 3, 6, 9, 12, 2, 5, 8, 15, 1, 4, 11, 14,
                             18, 21, 24, 31, 17, 20, 27, 30, 16, 23, 26, 29, 19, 22, 25, 28};
  /* Rotation around 4-fold axis */
  static long long R4[32] = {0, 13, 10, 7, 28, 25, 22, 19, 8, 5, 2, 15, 20, 17, 30, 27,
                             24, 21, 18, 31, 4, 1, 14, 11, 16, 29, 26, 23, 12, 9, 6, 3};
  /* Rotation around 2-fold axis */
  static long long R2[32] = {0, 7, 10, 13, 4, 11, 14, 1, 8, 15, 2, 5, 12, 3, 6, 9,
                             24, 31, 18, 21, 28, 19, 22, 25, 16, 23, 26, 29, 20, 27, 30, 17};

  long long j, s, sa, new_site;

  /*  i is an index in the symlist array */
  s = symlist[i];
  switch (s)
  {
  case IDENTITY:
    return site;
  case SPIN_FLIP:
    return (unsigned long long)-1; // Not a point symmetry
  case LINEAR_T:
    new_site = site + 1;
    if (new_site == Nspins)
      new_site = 0;
    return new_site;
  case FCC32Tx:                    /* Translation along line */
                                   /* NB! NB! "j" changed by "new_state". Possible error in FCC32 calc. ??! */
    return (unsigned long long)-1; // Not implemented
  case FCC32Ty:                    /* Translation between lines */
    return (unsigned long long)-1; // Not implemented
  case FCC32Tz:                    /* Translation between planes */
    return (unsigned long long)-1; // Not implemented
  case FCC32MIRROR:                /* Mirror reflection */
    return (unsigned long long)-1; // Not implemented
  case FCC32R2:                    /* Two-fold rotation */
    return (unsigned long long)-1; // Not implemented
  case FCC32R4:                    /* Four-fold rotation */
    return (unsigned long long)-1; // Not implemented
  case FCC32R3:                    /* Three-fold rotation */
    return (unsigned long long)-1; // Not implemented
  default:                         /* interpret added operators */
    sa = s - NSYM;
    if ((sa >= 0) && sa < Nsymadd)
    {
      return symadd[sa][site];
    }
    else
      fatalerror(" Unknown point symmetry operation ", s);
  }
  return 0; /* Never get here */
}

/* Initialise symmetry operations and check if requested */
void InitSym(struct FLAGS *input_flags)
{

  long long sym, Nsv[NSYM + NSYMADD], i, j;
  unsigned long long tmp, start = SPIN_0_UP;

  flipmask = ~(~((unsigned long long)0) << Nspins);
  if (flipmask == 0)
    flipmask = ~flipmask;
  if (input_flags->TEST_SYM)
  {
    LogMessageCharInt("flipmask= ", flipmask);
    LogMessageChar("\n");
  }
  Nsv[SPIN_FLIP] = 2;
  Nsv[FCC32Tx] = 4;
  Nsv[FCC32Ty] = 4;
  Nsv[FCC32Tz] = 2;
  Nsv[FCC32MIRROR] = 2;
  Nsv[FCC32R2] = 2;
  Nsv[FCC32R4] = 4;
  Nsv[FCC32R3] = 3;
  Nsv[LINEAR_T] = Nspins;
  Nsv[IDENTITY] = 1;
  for (i = 0; i < Nsymadd; i++)
  { /* Find translation period for added symmetry */
    tmp = start;
    LogMessageCharInt("Added symmetry : ", i);
    LogMessageCharInt(", round : ", 0);
    LogMessageCharInt(", state : ", tmp);
    LogMessageChar("\n");
    for (j = 1; (tmp = SymOp(Nsym + i, tmp)) != start && j < Nspins; j++)
    {
      LogMessageCharInt("Added symmetry : ", i);
      LogMessageCharInt(", round : ", j);
      LogMessageCharInt(", state : ", tmp);
      LogMessageChar("\n");
    }
    Nsv[NSYM + i] = j;
  }

  Nsym += Nsymadd;
  for (sym = 0; sym < Nsym; sym++)
  {
    Nsymvalue[sym] = Nsv[symlist[sym]];
    if (input_flags->TEST_SYM)
    {
      LogMessageChar3Vector("sym,symlist,Nsv : ", sym, symlist[sym], Nsv[symlist[sym]]);
      LogMessageChar("\n");
    }
  }

  if (input_flags->TEST_SYM)
  {
    TestSym(input_flags);
  }

  return;
}

// Used for debugging the symmetry operations
void TestSym(struct FLAGS *input_flags)
{
  long long sym, i, j, s, T[NSYM];
  unsigned long long new_state, state;

  LogMessageChar("Nsymvalue=(");
  for (i = 0; i < Nsym; i++)
    LogMessageInt(Nsymvalue[i]);
  LogMessageChar(" ) \n");

  if (Nsymadd > 0)
  {
    for (j = 0; j < Nspins; j++)
      LogMessageInt(symadd[0][j]);
    LogMessageChar("\n");
  }
  for (i = 0; i < Nsym; i++)
  {
    LogMessageCharInt(" Testing symmetry number ", symlist[i]);
    LogMessageChar(" (makes sense only for translations) \n");
    for (j = 0; j < Nspins; j++)
    {
      LogMessageCharInt(" Translates spin ", j);
      LogMessageCharInt(" into ", SYMNUM(i, j));
      LogMessageCharInt(", state ", ((unsigned long long)1) << j);
      LogMessageCharInt(" -> ", SymOp(i, ((unsigned long long)1) << j));
      LogMessageChar("\n");
    }
  }

  state = (0xABCDABCD & flipmask);
  LogMessageCharInt(" Operating on state ", state);
  LogMessageChar("\n");
  for (i = 0; i < Nsym; i++)
  {
    LogMessageInt(Nsymvalue[i]);
    LogMessageCharInt(" applications of symmetry ", symlist[i]);
    LogMessageChar(" gave \n (");
    for (j = 0; j < Nsymvalue[i]; j++)
      LogMessageInt(state = SymOp(i, state));
    LogMessageChar(" ) \n");
  }

  state = (0x3456789A & flipmask);
  /* Run through symmetry operations on state and notice
  if state reappears */
  TLOOP_BEGIN
  if (new_state == state)
  {
    LogMessageChar(" state encountered for T= (");
    for (j = 0; j < Nsym; j++)
      LogMessageInt(T[j]);
    LogMessageChar(") \n");
  }
  TLOOP_END
  if (new_state != state)
  {
    LogMessageCharInt(" Returned to ", new_state);
    LogMessageCharInt(" rather than to ", state);
    LogMessageChar(" \n ");
  }
  return;
}

void MakeSymCoup(struct FLAGS *input_flags)
/* MakeSymCoup constructs the full Hamiltonian from a general
   pattern and knowledge of the symmetries */
/* IN TESTING. TODO: Finish this function */
{
  long long i, j, Nc = Ncoup, T[NSYM];
  long long s0, s1, str, sym, new_s0, new_s1;

  for (i = 0; i < Ncoup; i++)
  {
    s0 = hamil_coup[i][0];
    s1 = hamil_coup[i][1];
    for (sym = 0; sym < Nsym; sym++)
      T[sym] = 0;
    if (input_flags->TEST_SYMCOUPLING)
    {
      LogMessageCharInt("Initial coupling ", i);
      LogMessageCharInt(" : ", s0);
      LogMessageCharInt(" -> ", s1);
      LogMessageChar("\n");
    }

    sym = 0;
    T[sym] = 1; /* don't copy the initial coupling pattern */
    new_s0 = s0;
    new_s1 = s1;

    while (T[Nsym - 1] < Nsymvalue[Nsym - 1])
    {
      new_s0 = SymOpSite(sym, new_s0);
      new_s1 = SymOpSite(sym, new_s1);

      if (input_flags->TEST_SYMCOUPLING)
      {
        LogMessageCharInt("coupling ", new_s0);
        LogMessageCharInt(" -> ", new_s1);
        LogMessageChar("\n");
      }
      hamil_coup[Nc][0] = new_s0;
      hamil_coup[Nc][1] = new_s1;
      Jzz[Nc] = Jzz[i];
      Jxy[Nc] = Jxy[i];
      Janis[Nc] = Janis[i];
      if (input_flags->dipole)
      {
        Jdip[Nc] = Jdip[i];
        r_vector[Nc][X] = r_vector[i][X];
        r_vector[Nc][Y] = r_vector[i][Y];
        r_vector[Nc][Z] = r_vector[i][Z];
        geom_13[Nc] = geom_13[i];
      }

      for (sym = 0; ++T[sym] == Nsymvalue[sym];)
        if (sym < Nsym - 1)
        {
          T[sym++] = 0;
          new_s0 = SymOpSite(sym, new_s0);
          new_s1 = SymOpSite(sym, new_s1);
        }
      Nc++;
    } /* while */
  }
  Ncoup = Nc;

  return;
}
