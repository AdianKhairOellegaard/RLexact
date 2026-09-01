/* Include file RLexact.h -
* Main include file for diagonalization program RLexact
* Last change: SJ 20.11.2017
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94
* Version 4.0, September 2017
*
============================================
*/
#include "Max_problem.h"
#ifndef HEADER // Firstly make sure to only include all this once.

#define HEADER
// Struct definition for loading file
struct FLAGS
{
  // The following long long ints are used as booleans.
  long long use_lanczos; // long long is here to allow for use of matchlines function
  long long use_exact_matrix;
  long long m_sym;           // M_symmetry present = 1, absent = 0
  long long dipole;          // DLC: This is defineable, but needs testing
  long long ring_exchange;   // DLC: This is defineable, but needs testing
  long long find_eigenstate; //
  long long find_cross;      // Find cross in terms of S^zz (q,w), S^xx (q,w) and S^yy(q,w)
  long long find_cross_pm;   // Find cross in S^+-(q,w) and S^-+(q,w)
                             //  instead of S^xx (q,w) and S^yy(q,w).
                             // Requires MSYM and FIND_CROSS.
  long long motive;          // DLC TODO: We need a description of the MOTIVE. Maybe
  long long find_mag;        // DLC TODO: Needs description. Debugging required! Should only be used WITHOUT MSYM SJ 20/11/17
  long long cross;           // DLC: I have literally no idea what this is

  // OUTPUT SPECIFIERS
  long long write_energies;
  long long write_states;        // Prints groundstate in dat-file
                                 // and all eigenstates if MATRIX;
  long long write_magnetisation; // DLC TODO: Needs description
                                 // Should only be used WITHOUT MSYM, Works only for MATRIX-mode, SJ 31/5/16

  // ============================================================
  // DEBUG / TEST FLAGS (converted from historical #define switches)
  // ============================================================

  // --- From RLexact.C / RLexact.c ---
  long long random_number_seed;
  long long TEST_GS_SEARCH;
  long long TEST_ALLOCATE;
  long long MAIN_LOOP_MESSAGES;

  // --- From RLhamil.c ---
  long long TEST_HAMILTON; // only in use in hamil.dipole
  long long TEST_MAT_ELEM;
  long long TEST_HAM2;
  long long TEST_HAMZEE;
  long long TEST_FILLHAM;

  // --- From RLsparse.C ---
  long long TEST_MAKESPARSE;
  long long TEST_APPLYSPARSE;
  long long TEST_APPLYSPARSE_READ;
  long long TEST_APPLYSPARSE_DEEP;
  long long DEBUG_APPLYSPARSE;
  long long TEST_HAM4;

  // --- From RLmatrix.c ---
  long long MATRIX_MESSAGES;
  long long TEST_EIG;
  long long TEST_EVEC;
  long long TEST_MATRIXMAG;
  long long TEST_CALC_M;
  long long TEST_CALC_SZZQ;

  // --- From RLlanczos.c ---
  long long TEST_SEED;
  long long TEST_EIGENVECTORS;
  long long TEST_FINDGROUND;
  long long TEST_STATES;
  long long TEST_ENERGIES;
  long long TEST_LANCZOS_VECTOR;
  long long LANCZOS_MESSAGES;
  long long VERBOSE_LANCZOS;
  long long TEST_FINDMAG;
  long long DEBUG_SEED; // WARNING: May cause premature stopping of lanczos algorithm
  long long DEBUG_LANCSTEP;
  long long TEST_LANCCROSS;

  // --- From RLtables.c ---
  long long TEST_TABLES;
  long long TEST_INVERTMATRIX;
  long long TEST_OCC;
  long long TEST_ISUNIQUE;
  long long TEST_FINDUNIQUE;
  long long TEST_FINDUNIQUE_DETAIL;
  long long TEST_FILLUNIQUE;
  long long TEST_READUNIQUE;
  long long TEST_FILLUNIQUE_LIST;
  long long TEST_FILLUNIQUEOBSERVABLES;
  long long TEST_COUNT;
  long long TEST_LOOKUP;
  long long TEST_CHECKPHASE;

  // --- From RLsymm.c ---
  long long TEST_SYM;
  long long TEST_SYMCOUPLING;

  // --- From RLio.c and RLutil.c ---
  long long TEST_WRITECROSS;
  long long TEST_SPINFLIP;
  long long TEST_INPUT;
  long long TEST_ROTATION;
  long long TEST_WRITEHMQ;

  // --- From RLcross.c ---
  long long TEST_APPLYSMP;
  long long TEST_APPLYSZQ;
  long long TEST_CROSS;

  // --- From Diagonal.C ---
  long long TEST_DIAGONALIZE;

  // --- From regcc.cpp ---
  long long TEST_MULTIMATCH;
  long long TEST_FILEREAD;

  // ============================================================
  // Existing verbosity controls (already present in your struct)
  // ============================================================
  long long VERBOSE_TIME_LV1;
  long long VERBOSE_TIME_LV2;
  long long VERBOSE;
  long long FILEREAD_VERBOSE;
};

// Highest possible state
//  #define MAX_STATE ((((unsigned long) 1)<<(Nspins-1)) -1)
// #else
#define MAX_STATE ((unsigned long long)1 << (Nspins))
// #endif

// Generated vectors with squared lengths below the following
// cutoff are considered to be zero in Lanczos calculation because of numerical precision.
#define ZERO_VEC_LENGTH 1e-9
#define SHORT_VEC_LENGTH (double)1 / NUNIQUE
#define MAX_LANCZOS 600
#define right(i) (((i & 1) << (Nspins - 1)) + (i >> 1))
#define inverse(i) (i ^ (MAX_STATE - 1))

// Symmetry related definitions
#define MASK_L1 0x77777777
#define MASK_L2 0x88888888
#define MASK_L3 0xF000F000
#define MASK_L4 0x0FFF0FFF
#define MASK_P1 0x0000FFFF
#define MASK_P2 0xFFFF0000

// Hard coded symmetry flags
#define SPIN_FLIP 0
#define FCC32Tx 1
#define FCC32Ty 2
#define FCC32Tz 3
#define FCC32MIRROR 4
#define FCC32R2 5
#define FCC32R4 6
#define FCC32R3 7
#define IDENTITY 9
#define LINEAR_T 10

#define SYMNUM(i, j) ((long long)(log((double)SymOp((long long)i, ((unsigned long long)1) << j) + 1E-6) / log((double)2)))

// Deal with complex numbers
#define komplex std::complex<double>
#define kvector cvector
#define freekvector freecvector
#define kmatrix cmatrix
#define freekmatrix freecmatrix
#define I komplex(0.0, 1.0) //_Complex_I
#define zero komplex(0.0, 0.0)
#define one komplex(1.0, 0.0)
#define Arg(a) (a == zero ? 0.0 : arg(a))
#define arg(a) atan2(imag(a), real(a))
#define abs(a) sqrt(real(a) * real(a) + imag(a) * imag(a))
#define sqrabs(a) (real(a) * real(a) + imag(a) * imag(a))
#define eksp(a) exp(real(a)) * (cos(imag(a)) + I * sin(imag(a)))
#define conj(a) (real(a) - I * imag(a))
#define skrt(a) (sqrt(abs(a)) * (cos(Arg(a) / 2.0) + I * sin(Arg(a) / 2.0)))

// Various definitions
#define SPIN_0_UP ((unsigned long long)1)
#define SQR(a) ((a) * (a))
#define RANDOM ((double)2 * rand() / RAND_MAX - 1)
#define PI 3.14159265358979323846
#define LARGE_NUMBER 9.999E+99
#define SMALL_NUMBER 1E-6
#define float double

/* Output definitions */
#define CSVOUT // output data to s** as comma separated values
#define FILEEND ".dat"
#define LOGFILEEND ".log"
#define SZZEND ".szz" // only used of FIND_CROSS
#define SXXEND ".sxx" // only used of FIND_CROSS
#define SYYEND ".syy" // only used of FIND_CROSS
#define SPMEND ".spm" // only used of FIND_CROSS and FIND_CROSS_PM
#define SMPEND ".smp" // only used of FIND_CROSS and FIND_CROSS_PM

#define MATRIXFILENAME "SqMat"
#define COEND ".gs"
#define GSCOEND ".gs"
#define UNIEND ".un"
#define UNIOBSEND ".uno"

// Flags
#define START 1
#define STOP 0
#define X 0
#define Y 1
#define Z 2
// #define REAL 0
// #define IMAG 1 //Not used and clashed with OpenMPI
#define NORMAL 0
#define RECONSTRUCT 1
#define CROSS 2
#define SZZ 0
#define SPM 1
#define SMP 2
// parallelization flags
#define MODEN 0    // Normal mode
#define MODEW 1    // Write Hamiltonian to file
#define MODEGS 2   // Groundstate mode
#define MODERC 3   // Reconstruct mode
#define MODEQ 4    // Cross section mode
#define UNIMODEN 0 // Normal mode
#define UNIMODEW 1 // Write unique
#define UNIMODER 2 // Read unique

// Program pieces
#define TLOOP_BEGIN                         \
  for (sym = 0; sym < Nsym; sym++)          \
    T[sym] = 0;                             \
  sym = 0, new_state = state;               \
  while (T[Nsym - 1] < Nsymvalue[Nsym - 1]) \
  {
#define TLOOP_END                            \
  new_state = SymOp(0, new_state);           \
  for (sym = 0; ++T[sym] == Nsymvalue[sym];) \
    if (sym < Nsym - 1)                      \
    {                                        \
      T[sym++] = 0;                          \
      new_state = SymOp(sym, new_state);     \
    }                                        \
  } /* while */
/*TRANSLOOP: for 1D this will be a loop over Nspins. statrs at sym=1 to avoid spin-flip/identity symmetry*/
#define TRANSLOOP_BEGIN                     \
  T[0] = 1;                                 \
  for (sym = 1; sym < Nsym; sym++)          \
    T[sym] = 0;                             \
  sym = 1, new_state = state;               \
  while (T[Nsym - 1] < Nsymvalue[Nsym - 1]) \
  {
#define TRANSLOOP_END                        \
  new_state = SymOp(1, new_state);           \
  for (sym = 1; ++T[sym] == Nsymvalue[sym];) \
    if (sym < Nsym - 1)                      \
    {                                        \
      T[sym++] = 0;                          \
      new_state = SymOp(sym, new_state);     \
    }                                        \
  } /* while */

#define QLOOP_BEGIN                         \
  for (sym = 0; sym < Nsym; sym++)          \
    q[sym] = 0;                             \
  sym = 0;                                  \
  while (q[Nsym - 1] < Nsymvalue[Nsym - 1]) \
  {
#define QLOOP_END                            \
  for (sym = 0; ++q[sym] == Nsymvalue[sym];) \
  {                                          \
    if (sym < Nsym - 1)                      \
      q[sym++] = 0;                          \
  }                                          \
  } /* while */

#endif
