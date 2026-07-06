/* Program file RLio.c -
* Processes all user input/output
* and some geometric conversions
*
* Last change: SJ 19.12.16
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
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex>
#include <cnr.h>
#include "Functions.h"
#include <RLexact.h>

#include <mpi.h>
extern int rank, nprocs;

/* Functions defined in this file */

void TransformCoup(long long);
#ifdef FIND_MAG
void WriteMaggs(long long *);
#endif // FIND_MAG
void time_stamp(time_t *, long long, const char *);
void Warning(const char *, long long);
void LogMessageChar(const char *);
void OutMessageChar(const char *);
void LogMessageInt(long long);
void LogMessageImag(long long);
void LogMessageCharDouble(const char *, double);
void LogMessageCharInt(const char *, long long);
void OutMessageCharInt(const char *, long long);
void LogMessageChar3Vector(const char *, double, double, double);
void WriteState(const char *, komplex *);
void WriteStates(komplex **);

void WriteGSEnergy(komplex);
void WriteEnergy(double);
void WriteQvalue(long long *);
void WriteGSdata(double, long long *);
void WriteGSstate(komplex *);
void ReadGSdata(double *, long long *, komplex *);
void ReadGSenergy(double *, long long *);

/* Functions in RLutils.C */
extern void NormalizeVector(double *);
extern void RotateVector(double *);
extern void Bubblesort(double *, double *, long long);
extern void FillRotationMatrix(double *);

/* Functions defined elsewhere */
extern void MakeSymCoup();
extern void itoa(long long, char[]);
extern void filereader(char *, char *, long long);
extern long long filesizer(char *);
extern long long multimatch(char *, long long, const char *, double **, long long *, long long);
extern long long multimatch(char *, long long, const char *, long long **, long long *, long long);
extern long long matchlines(char *, const char *, double *, bool);
extern long long matchlines_wrapper(char *, const char *, long long *, bool);
extern double *dvector(long long, long long);

/* Global variables defined in RLexact.c */
extern long long Nspins, Nsym, Nsymadd, Nunique, Nuniq_k;
extern long long **symadd;
extern long long Ndimensions;
extern long long *TransIds;
extern long long Trans_Qmax[3];
extern long long hamil_coup[NCOUP][2], Ncoup;
extern long long ring_coup[NCOUP][4], Nring;
extern komplex *gs;
extern long long twom;
extern double mstart, mend;
extern double h, hstart, hend, hstep;
extern double field[3];
extern double Jzz[NCOUP], Jxy[NCOUP], Janis[NCOUP];
extern double Jr[NCOUP];
extern double Jdip[NCOUP], geom_13[NCOUP], r_vector[NCOUP][3];
extern long long symlist[NSYM];
extern double *energies;
#ifdef FIND_MAG
extern double *magnetisation;
extern double maggs;
#endif /* FIND_MAG */
extern long long Nq_choice;
extern long long **q_choice;
extern FILE *outfilezz;
extern FILE *outfilexx, *outfileyy;
extern FILE *outfilepm, *outfilemp;
extern double *cross;
#ifdef MOTIVE
extern long long Nspins_in_uc;
extern float **spin_positions;
#endif // MOTIVE
extern FILE *gscoinfile;
extern FILE *outfile;
extern FILE *logfile;
extern long long q_gs[NSYM];
extern double gs_energy;
extern char *infile_name;
extern bool name_on_commandline;
extern double szqlength;
extern long long mode;
extern long long unimode;
extern int spinflip_number;
extern int spinflip_present;
/* Lanczos variables defines in RLlanz.C */
extern double Ritz_conv;
extern unsigned long long *unique;

/* Regional variables defined here */
FILE *infile;
FILE *gscofile;
char *filedata;
long long filesize;

/************************************************/
long long intro(struct FLAGS *input_flags)
{
  /* Introduce, read input and perform file handling */
  char outfile_name[256];
  char logfile_name[256];
  char gscofile_name[256];
  char gscoinfile_name[256];
  char qstr[2];

  if (!name_on_commandline)
  {
    int namelen;
    if (rank == 0)
    {
      OutMessageChar(" Please type the name of the input file : ");
      scanf("%s", infile_name);
      OutMessageChar(" ... \n");
      namelen = strlen(infile_name) + 1;
    }
    MPI_Bcast(&namelen, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(infile_name, namelen, MPI_CHAR, 0, MPI_COMM_WORLD);
  }

  // strcpy(outfile_name,infile_name);
  // strcat(outfile_name,FILEEND);
  errno = 1;
  if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, FILEEND) >= 256)
  {
    fatalerror("infile_name too large", errno);
    return -1;
  }
  outfile = fopen(outfile_name, "w");
  if (outfile == NULL)
  {
    fatalerror("Cannot open general output file, sorry!", errno);
    return -1;
  }
  fflush(outfile);
  // strcpy(logfile_name,infile_name);
  // strcat(logfile_name,LOGFILEEND);
  errno = 1;
  if (snprintf(logfile_name, 255, "%s-%d%s", infile_name, rank, LOGFILEEND) >= 256)
  {
    fatalerror("infile_name too large", errno);
    return -1;
  }
  logfile = fopen(logfile_name, "w");
  if (logfile == NULL)
  {
    fatalerror("Cannot open logfile, sorry!", errno);
    return -1;
  }
  ReadInputFlags(infile_name, input_flags);
  if (rank == 0)
  {

    OutMessageChar("Welcome to the Exact Diagonalization Program, RLexact \n");
    if (input_flags->use_exact_matrix)
      OutMessageChar("Using Matrix");
    else if (input_flags->use_lanczos)
      OutMessageChar("Using Lanczos");
    OutMessageCharInt("diagonalization. BUFFERSIZE =", BUFFERSIZE);
    if (input_flags->m_sym)
      OutMessageChar("M-symmetry present. \n");
    else
      OutMessageChar(" M-symmetry absent. \n");
    OutMessageChar(" Observables:");
#ifdef FIND_MAG
    OutMessageChar(" Magnetization,");
#endif /* FIND_MAG */
#ifdef CROSS
    OutMessageChar(" S^zz(q,w),");
    if (!input_flags->m_sym)
    {
      if (!input_flags->find_cross_pm)
        OutMessageChar(" S^xx(q,w), S^yy(q,w),");
      else
        OutMessageChar(" S^+-(q,w), S^-+(q,w),");
#endif
    }
    OutMessageChar(" Energy.\n");
    OutMessageChar(" For more information, see the manual.\n");
  }

  ReadCoupPattern(infile_name, input_flags); // this function does the actual file-input handling

  if (input_flags->find_cross)
  {
    if (rank == 0)
    {
      if (mode == MODEQ || mode == MODERC)
      {
        strcpy(gscoinfile_name, infile_name);
        strcat(gscoinfile_name, GSCOEND);
        errno = 1;
        gscoinfile = fopen(gscoinfile_name, "r");
        if (gscoinfile == NULL)
        {
          fatalerror("Cannot open gs file, sorry!", errno);
          return -1;
        }
        fflush(gscoinfile);
      }

      if (mode == MODEGS)
      {
        strcpy(gscofile_name, infile_name);
        strcat(gscofile_name, COEND);
        errno = 1;
        gscofile = fopen(gscofile_name, "w");
        if (gscofile == NULL)
        {
          fatalerror("Cannot open carry over file, sorry!", errno);
          return -1;
        }
        fflush(gscofile);
      }
      if (mode == MODERC)
      {
        strcpy(gscofile_name, infile_name);
        strcat(gscofile_name, COEND);
        errno = 1;
        gscofile = fopen(gscofile_name, "a");
        if (gscofile == NULL)
        {
          fatalerror("Cannot open carry over file, sorry!", errno);
          return -1;
        }
        fflush(gscofile);
      }
    } // Not MODEN - Perl scripts are lost so I don't know how/if these should be
      // parallelised. ABP - 2025/03/12

    // strcpy(outfile_name,infile_name);
    // strcat(outfile_name,SZZEND);
    errno = 1;
    if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, SZZEND) >= 256)
    {
      fatalerror("infile_name too large", errno);
      return -1;
    }
    if ((outfilezz = fopen(outfile_name, "w")) == NULL)
    {
      fatalerror("Cannot open output file for Szz, sorry!", errno);
      return -1;
    }
    fflush(outfilezz);

    if (!input_flags->find_cross_pm)
    {
      // strcpy(outfile_name,infile_name);
      // strcat(outfile_name,SXXEND);
      errno = 1;
      if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, SXXEND) >= 256)
      {
        fatalerror("infile_name too large", errno);
        return -1;
      }
      if ((outfilexx = fopen(outfile_name, "w")) == NULL)
      {
        fatalerror("Cannot open output file for S+-, sorry!", errno);
        return -1;
      }
      fflush(outfilexx);

      // strcpy(outfile_name,infile_name);
      // strcat(outfile_name,SYYEND);
      errno = 1;
      if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, SYYEND) >= 256)
      {
        fatalerror("infile_name too large", errno);
        return -1;
      }
      if ((outfileyy = fopen(outfile_name, "w")) == NULL)
      {
        fatalerror("Cannot open output file for S-+, sorry!", errno);
        return -1;
      }
      fflush(outfileyy);
    }

    if (input_flags->find_cross_pm)
    {
      // strcpy(outfile_name,infile_name);
      // strcat(outfile_name,SPMEND);
      errno = 1;
      if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, SPMEND) >= 256)
      {
        fatalerror("infile_name too large", errno);
        return -1;
      }
      if ((outfilepm = fopen(outfile_name, "w")) == NULL)
      {
        fatalerror("Cannot open output file for S+-, sorry!", errno);
        return -1;
      }
      fflush(outfilepm);

      // strcpy(outfile_name,infile_name);
      // strcat(outfile_name,SMPEND);
      if (snprintf(outfile_name, 255, "%s-%d%s", infile_name, rank, SMPEND) >= 256)
      {
        fatalerror("infile_name too large", errno);
        return -1;
      }
      errno = 1;
      if ((outfilemp = fopen(outfile_name, "w")) == NULL)
      {
        fatalerror("Cannot open output file for S-+, sorry!", errno);
        return -1;
      }
      fflush(outfilemp);
    }
  }
#ifdef TEST_INPUT
  LogMessageChar("All filenames are OK. \n");
#endif /* TEST_INPUT */

  return 0; /* All is OK */
}

/* ----------------------------------------------------------------------- */
void ReadInputFlags(char *filename, struct FLAGS *input_flags)
{
  // This function simply reads in the input file, and captures all the flag
  // settings.
  filesize = filesizer(filename);
  filedata = (char *)malloc(filesize * sizeof(char));
  if (filedata == NULL)
  {
    printf("\nERROR: Filedata not allocated");
    exit(1);
  }
  filereader(filename, filedata, filesize); // the entire file is now in filedata
  input_flags->use_lanczos = 1;             // Using lanczos is default.
  input_flags->use_exact_matrix = 0;
  input_flags->m_sym = 1; // Use m_sym as default
  input_flags->dipole = 0;
  input_flags->ring_exchange = 0;
  input_flags->find_eigenstate = 1; // As default, find eigenstates
  input_flags->find_cross = 1;
  input_flags->find_cross_pm = 0;

  input_flags->write_energies = 1; // Output energies and states as default
  input_flags->write_states = 1;

  matchlines_wrapper(filedata, "Use_Lanczos", &input_flags->use_lanczos, true);
  matchlines_wrapper(filedata, "Use_Exact_Matrix", &input_flags->use_exact_matrix, true);
  matchlines_wrapper(filedata, "M_Symmetry", &input_flags->m_sym, true);
  matchlines_wrapper(filedata, "Dipole", &input_flags->dipole, true);
  matchlines_wrapper(filedata, "Ring_exchange", &input_flags->ring_exchange, true);
  matchlines_wrapper(filedata, "Find_Eigenstate", &input_flags->find_eigenstate, true);
  matchlines_wrapper(filedata, "Find_cross", &input_flags->find_cross, true);
  matchlines_wrapper(filedata, "Find_cross_pm", &input_flags->find_cross_pm, true);

  matchlines_wrapper(filedata, "Write_Energies", &input_flags->write_energies, true);
  matchlines_wrapper(filedata, "Write_States", &input_flags->write_states, true);

  matchlines_wrapper(filedata, "VERBOSE_TIME_LV1", &input_flags->VERBOSE_TIME_LV1, true);
  matchlines_wrapper(filedata, "VERBOSE_TIME_LV2", &input_flags->VERBOSE_TIME_LV2, true);
  matchlines_wrapper(filedata, "VERBOSE", &input_flags->VERBOSE, true);
}

/* ----------------------------------------------------------------------- */

long long ReadCoupPattern(char *filename, struct FLAGS *input_flags)
{
  /* Read all couplings (pairs, types and strengths)
     and symmetry numbers from file */
  long long Ncoupstr; /* Number of different coupling strengths */
  long long Nringstr;
  double hamring[NRINGSTR];
  long long n1, n2, str, i, j, symconstruct;
  double hamzz[NCOUPSTR], hamxy[NCOUPSTR], hamanis[NCOUPSTR];
  double r, rx, ry, rz;
  double Dip[NCOUPSTR];
  double B2;
  char test[80];
  long long filesize; // of inputfile
  char *filedata;

  if (rank == 0)
  {
    filesize = filesizer(filename);
    filedata = (char *)malloc(filesize * sizeof(char));
    if (filedata == NULL)
    {
      printf("\nERROR: Filedata not allocated");
      exit(1);
    }
    filereader(filename, filedata, filesize); // the entire file is now in filedata

#ifdef TEST_INPUT
    LogMessageChar("Input file opened...\n");
#endif /* TEST_INPUT */

    matchlines_wrapper(filedata, "Number of spins", &Nspins, true);
#ifdef TEST_INPUT
    LogMessageCharInt(" Nspins:", Nspins);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }
  MPI_Bcast(&Nspins, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  // *********** Mandatory input: Mode *************
  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Mode", &mode, true);
#ifdef TEST_INPUT
    switch (mode)
    {
    case MODEN:
      LogMessageChar("Mode is: NORMAL MODE. \n");
      break;
    case MODEGS:
      LogMessageChar("Mode is: GS MODE. \n");
      break;
    case MODEQ:
      LogMessageChar("Mode is: Q MODE. \n");
      break;
    default:
      fatalerror("Unknown mode", mode);
    }
#endif /* TEST_INPUT */
  }
  MPI_Bcast(&mode, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    // *********** Mandatory input: Unique mode *************
    matchlines_wrapper(filedata, "Unimode", &unimode, true);
#ifdef TEST_INPUT
    switch (unimode)
    {
    case UNIMODEN:
      LogMessageChar("Unique mode is: NORMAL MODE. \n");
      break;
    case UNIMODEW:
      LogMessageChar("Unique mode is: WRITE MODE. \n");
      break;
    case UNIMODER:
      LogMessageChar("Unique mode is: READ MODE. \n");
      break;
    default:
      fatalerror("Unknown unique mode", unimode);
    }
#endif /* TEST_INPUT */
  }
  MPI_Bcast(&unimode, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  // ******** Mandatory input: Read number of lines to input **********
  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Number of couplings", &Ncoup, true);
    matchlines_wrapper(filedata, "Number of coupling strengths", &Ncoupstr, true);
  }
  MPI_Bcast(&Ncoup, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ncoupstr, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  if (input_flags->ring_exchange)
  {
    if (rank == 0)
    {
      matchlines_wrapper(filedata, "Number of rings", &Nring, true);
      matchlines_wrapper(filedata, "Number of ringstrength", &Nringstr, true);
    }
    MPI_Bcast(&Nring, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Nringstr, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  }
#ifdef TEST_INPUT
  LogMessageCharInt("Ncoup:", Ncoup);
  LogMessageCharInt("Ncoupstr:", Ncoupstr);
  LogMessageChar("\n");
  if (input_flags->ring_exchange)
  {
    LogMessageCharInt("Nring:", Nring);
    LogMessageCharInt("Nringstr:", Nringstr);
    LogMessageChar("\n");
  }
#endif /* TEST_INPUT */

  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Number of hardcoded symmetries", &Nsym, true);
    matchlines_wrapper(filedata, "Number of custom symmetries", &Nsymadd, true);
    matchlines_wrapper(filedata, "Construct symmetries", &symconstruct, true);
    if ((Nsym + Nsymadd) <= 0)
    {
      printf("\nNsym=%lld, Nsymadd=%lld\n", Nsym, Nsymadd);
      fatalerror("At least one symmetry must be defined; use for instance IDENTITY", IDENTITY);
    }
  }
  MPI_Bcast(&Nsym, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nsymadd, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&symconstruct, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  LogMessageCharInt("Number of hardcoded symmetries: Nsym=", Nsym);
  LogMessageCharInt("\nNumber of custom symmetries: Nsymadd=", Nsymadd);
  LogMessageCharInt("\nsymconstruct:", symconstruct);
  LogMessageChar("\n");
#endif /* TEST_INPUT */

  // ********* Input symmetry info ********************
  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Hardcoded symmetries", symlist, false);
  }
  MPI_Bcast(&symlist, Nsym, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  LogMessageCharInt("Number of hardcoded symmetries scanned:", Nsym);
  LogMessageChar("\n");
#endif /* TEST_INPUT */

  symadd = (long long **)malloc(Nsymadd * sizeof(long long *));
  for (int i = 0; i < Nsymadd; i++)
    symadd[i] = (long long *)malloc(Nspins * sizeof(long long));

  if (rank == 0)
  {
    long long *dummy = (long long *)malloc(Nsymadd * sizeof(long long));
    multimatch(filedata, filesize, "Custom symmetry", symadd, dummy, Nsymadd);
    free(dummy);
  }
  for (int i = 0; i < Nsymadd; i++)
    MPI_Bcast(symadd[i], Nspins, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Number of dimensions", &Ndimensions, true);
  }
  MPI_Bcast(&Ndimensions, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  TransIds = (long long *)malloc(3 * sizeof(long long));
  if (rank == 0)
  {
    for (int i = 0; i < 3; i++)
      TransIds[i] = 0; // Ugly,but needed for loop in RLcross.
    matchlines_wrapper(filedata, "Translation indices", TransIds, 1);
  }

  MPI_Bcast(TransIds, 3, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  LogMessageCharInt("Number of dimensions:", Ndimensions);
  LogMessageChar("\n");
  LogMessageChar("Translations at indices:");
  for (int i = 0; i < Ndimensions; i++)
    LogMessageInt(TransIds[i]);
  LogMessageChar("\n");
#endif

#ifdef TEST_SPINFLIP
  LogMessageCharInt("Testing for spin flip sym. Matching symlist with SPIN_FLIP=", SPIN_FLIP);
  LogMessageChar("\n symlist=");
#endif

  for (i = 0; i < Nsymadd; i++) // Give numbers to custom symmetries
  {
    symlist[Nsym + i] = NSYM + i;

#ifdef TEST_SPINFLIP
    LogMessageCharInt("\nsymlist[Nsym+i]=", symlist[Nsym + i]);
#endif
  }

#ifdef TEST_SPINFLIP
  LogMessageChar("\n");
#endif

  spinflip_number = -1;
  spinflip_present = 0;
  for (i = 0; i < Nsym; i++) // Look for spin flip
  {
#ifdef TEST_SPINFLIP
    LogMessageCharInt("\nMatch SPIN_FLIP with symlist[Nsym+i]=", symlist[Nsym + i]);
#endif
    if (symlist[i] == SPIN_FLIP)
    {
      spinflip_present = 1;
      spinflip_number = i;
#ifdef TEST_SPINFLIP
      LogMessageCharInt("Spinflip found! at i =", i);
      LogMessageChar("\n");
#endif
    }
  }
#ifdef TEST_SPINFLIP
  LogMessageCharInt("\nRLio: spinflip_number=", spinflip_number);
#endif

  // ************ Input GS q-values ******************
  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Number of chosen GS q-values", &Nq_choice, true);
  }
  MPI_Bcast(&Nq_choice, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
#ifdef TEST_INPUT
  LogMessageCharInt("Nq_choice:", Nq_choice);
#endif /* TEST_INPUT */

  q_choice = (long long **)malloc(Nq_choice * sizeof(long long *));
  if (rank == 0)
  {
    long long *dummy = (long long *)malloc(Nq_choice * sizeof(long long));
    multimatch(filedata, filesize, "Chosen GS q-value", q_choice, dummy, Nq_choice);
    free(dummy);
  }
  MPI_Bcast(q_choice, Nq_choice, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  for (i = 0; i < Nq_choice; i++)
  {
    LogMessageChar("q_choice: ( ");
    for (j = 0; j < Nsym + Nsymadd; j++)
    {
      LogMessageCharInt(" ", q_choice[i][j]);
    }
    LogMessageChar(") \n");
  }
#endif /* TEST_INPUT */

#ifdef MOTIVE
  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Number of spins in unit cell", &Nspins_in_uc, true);
  }
  MPI_Bcast(&Nspins_in_uc, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  LogMessageCharInt("Number of spins in unit cell: ", Nspins_in_uc);
  LogMessageChar("\n");
#endif // TEST_INPUT

  spin_positions = (double **)malloc(Nspins_in_uc * sizeof(double *));
  for (int i = 0; i < Nspins_in_uc; i++)
  {
    spin_positions[i] = (double *)malloc(3 * sizeof(double));
  }
  if (rank == 0)
  {
    long long *dummy = (long long *)malloc(Nspins_in_uc * sizeof(long long));
    multimatch(filedata, filesize, "Relative position", spin_positions,
               dummy, Nspins_in_uc);
    free(dummy);
  }

  for (int i = 0; i < Nspins_in_uc; i++)
  {
    MPI_Bcast(spin_positions[i], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

#ifdef TEST_INPUT
  for (i = 0; i < Nspins_in_uc; i++)
  {
    LogMessageChar3Vector("Spin pos.", spin_positions[i][X],
                          spin_positions[i][Y],
                          spin_positions[i][Z]);
  }
#endif // TEST_INPUT

  if (rank == 0)
  {
    matchlines_wrapper(filedata, "Qmax translation", Trans_Qmax, true);
  }
  MPI_Bcast(Trans_Qmax, 3, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
#ifdef TEST_INPUT
  LogMessageChar3Vector("Qmax translation", Trans_Qmax[X],
                        Trans_Qmax[Y],
                        Trans_Qmax[Z]);

#endif

#endif // MOTIVE

  if (input_flags->m_sym)
  {

    if (rank == 0)
    {
      matchlines(filedata, "M start", &mstart, true);
      matchlines(filedata, "M end", &mend, true);
    }
    MPI_Bcast(&mstart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
    LogMessageCharDouble("mstart:", mstart);
    LogMessageCharDouble("mend:", mend);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }
  else
  {

    if (rank == 0)
    {
      matchlines(filedata, "H start", &hstart, true);
      matchlines(filedata, "H end", &hend, true);
      matchlines(filedata, "H step", &hstep, true);
    }
    MPI_Bcast(&hstart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hstep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
    LogMessageCharDouble("hstart:", hstart);
    LogMessageCharDouble("hend:", hend);
    LogMessageCharDouble("hstep:", hstep);
    LogMessageChar("\n");
#endif /* TEST_INPUT */

    if (rank == 0)
    {
      matchlines(filedata, "Hx", &field[X], true);
      matchlines(filedata, "Hy", &field[Y], true);
      matchlines(filedata, "Hz", &field[Z], true);
    }
    MPI_Bcast(field, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
    LogMessageChar3Vector("Field direction:", field[X], field[Y], field[Z]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
    for (B2 = 0, j = X; j <= Z; j++)
      B2 += SQR(field[j]);
    if (B2 > 0)
      for (j = X; j <= Z; field[j++] /= sqrt(B2))
        ;
    else
      fatalerror(" Field direction not well defined", 0);
#ifdef TEST_INPUT
    LogMessageChar3Vector("Normalized field direction:", field[X],
                          field[Y],
                          field[Z]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
    // FillRotationMatrix(field); //doesnt work and unnecessary, SJ 20.02.17
  }
  if (input_flags->use_lanczos)
  {
    // ********* Lanzcos-Mandatory input : Read Lanzcos numbers *********
    if (rank == 0)
    {
      matchlines(filedata, "Ritz_conv", &Ritz_conv, true);
    }
    MPI_Bcast(&Ritz_conv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
    LogMessageCharDouble("Ritz_conv:", Ritz_conv);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }

  if (rank == 0)
  {
    long long *dummy = (long long *)malloc(Ncoupstr * sizeof(long long));
    double **dummyresdouble = (double **)malloc(Ncoupstr * sizeof(double *));
    multimatch(filedata, filesize, "Coupling strength vector", dummyresdouble,
               dummy, Ncoupstr);

    for (long long k = 0; k < Ncoupstr; k++)
    {
      hamzz[k] = dummyresdouble[k][0];
      hamxy[k] = dummyresdouble[k][1];
      hamanis[k] = dummyresdouble[k][2];
      if (input_flags->dipole)
      {
        Dip[k] = dummyresdouble[k][3];
        printf("Coupling strengths: %g, %g, %g, %g",
               hamzz[k], hamxy[k],
               hamanis[k], Dip[k]);
      }
    }
    free(dummy);
    free(dummyresdouble);
  }
  for (long long k = 0; k < Ncoupstr; k++)
  {
    MPI_Bcast(&hamzz[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hamxy[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hamanis[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (input_flags->dipole)
      MPI_Bcast(&Dip[k], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
    LogMessageChar3Vector(" Jzz, Jxy, Janis:", hamzz[k], hamxy[k], hamanis[k]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }

  if (input_flags->ring_exchange)
  {
    if (rank == 0)
    {
      long long *dummy = (long long *)malloc(Nringstr * sizeof(long long));
      double **dummyresdouble1 = (double **)malloc(Nringstr * sizeof(double *));
      multimatch(filedata, filesize, "Ring strength", dummyresdouble1, dummy,
                 Nringstr);

      for (long long k = 0; k < Nringstr; k++)
      {
        hamring[k] = dummyresdouble1[k][0];
#ifdef TEST_INPUT
        LogMessageCharDouble(" Jr:", hamring[k]);
        LogMessageChar("\n");
#endif /* TEST_INPUT */
      }
      free(dummy);
      free(dummyresdouble1);
    }
    MPI_Bcast(hamring, Nringstr, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  // MPI Not tested, since documentation is sparse on this functionality - ABP

  if (rank == 0)
  {
    long long *dummy = (long long *)malloc(Ncoup * sizeof(long long));
    long long **dummyres = (long long **)malloc(Ncoup * sizeof(long long *));
    multimatch(filedata, filesize, "Coupling vector", dummyres, dummy, Ncoup);

    for (long long k = 0; k < Ncoup; k++)
    {
      hamil_coup[k][0] = dummyres[k][0];
      hamil_coup[k][1] = dummyres[k][1];
      Jzz[k] = hamzz[dummyres[k][2]];
      Jxy[k] = hamxy[dummyres[k][2]];
      Janis[k] = hamanis[dummyres[k][2]];
    }
    if (input_flags->dipole)
    {
      for (long long k = 0; k < Ncoup; k++)
      {
        Jdip[k] = Dip[dummyres[k][2]];
        r_vector[k][X] = dummyres[k][3];
        r_vector[k][Y] = dummyres[k][4];
        r_vector[k][Z] = dummyres[k][5];
        NormalizeVector(r_vector[k]);
        RotateVector(r_vector[k]);
        geom_13[k] = 1.0 - 3.0 * SQR(r_vector[k][Z]);
#ifdef TEST_ROTATION
        LogMessageChar3Vector("   transformed direction: ", r_vector[k][X],
                              r_vector[k][Y],
                              r_vector[k][Z]);
        LogMessageChar("\n");
#endif /* TEST_ROTATION */
      }
      // Not parallelised Dipole is missing documentation
    }
    free(dummy);
    free(dummyres);
  }

  for (long long k = 0; k < Ncoup; k++)
  {
    MPI_Bcast(hamil_coup[k], 2, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(Jzz, Ncoup, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(Jxy, Ncoup, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(Janis, Ncoup, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef TEST_INPUT
  for (long long k = 0; k < Ncoup; k++)
  {
    LogMessageCharInt(" Coupling from ", hamil_coup[k][0]);
    LogMessageCharInt(", to ", hamil_coup[k][1]);
    LogMessageChar(", with strength: ");
    LogMessageChar("\n\t");
    LogMessageCharDouble(" Jzz ", Jzz[k]);
    LogMessageCharDouble(", Jxy ", Jxy[k]);
    LogMessageCharDouble(", Janis ", Janis[k]);
    LogMessageChar("\n");
  }
#endif /* TEST_INPUT */

  if (input_flags->ring_exchange)
  {
    if (rank == 0)
    {
      long long *dummy = (long long *)malloc(Nring * sizeof(long long));
      long long **dummyresring = (long long **)malloc(Nring * sizeof(long long *));
      multimatch(filedata, filesize, "Coupling ring vector", dummyresring, dummy,
                 Nring);

      for (long long k = 0; k < Nring; k++)
      {
        ring_coup[k][0] = dummyresring[k][0];
        ring_coup[k][1] = dummyresring[k][1];
        ring_coup[k][2] = dummyresring[k][2];
        ring_coup[k][3] = dummyresring[k][3];
        Jr[k] = hamring[dummyresring[k][4]];
#ifdef TEST_INPUT
        LogMessageCharInt(" Coupling from:", ring_coup[k][0]);
        LogMessageCharInt(", over:", ring_coup[k][1]);
        LogMessageCharInt(", and:", ring_coup[k][2]);
        LogMessageCharInt(", to:", ring_coup[k][3]);
        LogMessageCharInt(", and back to:", ring_coup[k][0]);
        LogMessageCharDouble(" Jr: ", Jr[k]);
        LogMessageChar("\n");
#endif /* TEST_INPUT */
      }
      free(dummy);
      free(dummyresring);
    } // Not parallelised Ring_exchange is missing documentation
  }

  if (rank == 0)
  {
    if (symconstruct == 1)
      MakeSymCoup(input_flags);

    LogMessageChar("RLio.C successful. \n");
    free(filedata);
  }
  return 0;
}

#ifndef NEVER // it doesnt work currently: rotating J is the wrong way!
void TransformCoup(long long j)
{
  double vec[3];

  vec[Z] = Jzz[j];
  vec[X] = Jxy[j] + Janis[j];
  vec[Y] = Jxy[j] - Janis[j];
  RotateVector(vec);
  Jzz[j] = vec[Z];
  Jxy[j] = (vec[X] + vec[Y]) / 2.0;
  Janis[j] = (vec[X] - vec[Y]) / 2.0;
  LogMessageCharDouble("Now Jz =", Jzz[j]);
  LogMessageCharDouble(", Jxy =", Jxy[j]);
  LogMessageCharDouble(", Janis =", Janis[j]);
  LogMessageChar("\n");
  return;
}
#endif /* NEVER */

long long SortCross(long long Nener)
{
  /* Finds and sums crosssections for same energies, and sorts energies and
   * cross by energy.
   * Could be optimized greatly, but MAX_LANCZOS is small: no point */
  long long XS = 0; // number of unique crosssections
  bool newener;

  double *tmpener = dvector(0, MAX_LANCZOS - 1);
  double *tmpcross = dvector(0, MAX_LANCZOS - 1);

  for (long long i = 0; i < Nener; i++)
  {
    newener = true;
    for (long long j = 0; j < XS; j++)
    {
      if (energies[i] == tmpener[j])
      {
        newener = false;
        tmpcross[j] += cross[i];
      }
    }
    if (newener)
    {
      tmpener[XS] = energies[i];
      tmpcross[XS++] = cross[i];
    }
  }

  Bubblesort(tmpener, tmpcross, XS); // Maybe not necessary, but nice

  for (long long i = 0; i < XS; i++)
  {
    cross[i] = tmpcross[i];
    energies[i] = tmpener[i];
  }
  freedvector(tmpener, 0, MAX_LANCZOS - 1);
  freedvector(tmpcross, 0, MAX_LANCZOS - 1);
  return XS;
}

void time_stamp(time_t *tim, long long flag, const char *string)
{
  if (flag == START)
  {
    *tim = -clock();
    fprintf(logfile, "Starting %s ... ", string);
  }
  if (flag == STOP)
  {
    *tim += clock();
    fprintf(logfile, " %s done, %lg seconds used\n",
            string, *tim / (double)CLOCKS_PER_SEC);
  }

  return;
}

void outro(struct FLAGS *input_flags)
{
  fclose(outfile);
  if (input_flags->find_cross)
  {
    fclose(outfilezz);
    if (!input_flags->find_cross_pm)
    {
      fclose(outfilexx);
      fclose(outfileyy);
    }
    else
    {
      fclose(outfilepm);
      fclose(outfilemp);
    }
  }
  LogMessageChar("\n End of diagonalization program RLexact.\n");

  return;
}

void fatalerror(const char *str, long long i)
{
  fprintf(logfile, "*** FATAL ERROR *** %s %lld %s \n", str, i, strerror(i));
  printf("*** FATAL ERROR *** %s %lld %s \n", str, i, strerror(i));
  exit(-1);
}

void Warning(const char *str, long long i)
{
  fprintf(logfile, "*** WARNING *** %s %lld \n", str, i);
  fflush(logfile);
  return;
}

void LogMessageInt(const long long i)
{
  fprintf(logfile, " %lld ", i);
  fflush(logfile);
  return;
}

void LogMessageImag(const double a, const double b)
{
  fprintf(logfile, " (  %lg + %lg i )", a, b);
  fflush(logfile);
  return;
}

void LogMessageKomplex(const komplex z)
{
  fprintf(logfile, " (  %lg + %lg i )", real(z), imag(z));
  fflush(logfile);
  return;
}

void LogMessageChar(const char *str)
{
  fprintf(logfile, " %s ", str);
  fflush(logfile);
  return;
}
void OutMessageChar(const char *str)
{
  fprintf(stdout, " %s ", str);
  return;
}

void LogMessageCharDouble(const char *str, double d)
{
  fprintf(logfile, " %s %g ", str, d);
  fflush(logfile);
  return;
}

void LogMessageCharKomplex(const char *str, komplex z)
{
  fprintf(logfile, " %s (  %lg + %lg i )", str, real(z), imag(z));
  fflush(logfile);
}

void LogMessageCharInt(const char *str, long long i)
{
  fprintf(logfile, " %s %lld ", str, i);
  fflush(logfile);
  return;
}
void OutMessageCharInt(const char *str, long long i)
{
  fprintf(stdout, " %s %lld ", str, i);
  return;
}

void LogMessageChar3Vector(const char *str, double d1, double d2, double d3)
{
  fprintf(logfile, " %s (%g, %g, %g ) \n", str, d1, d2, d3);
  fflush(logfile);
  return;
}

/* Overload for longs */
void LogMessageChar3Vector(const char *str,
                           long long l1, long long l2, long long l3)
{
  fprintf(logfile, " %s (%lld, %lld, %lld ) \n", str, l1, l2, l3);
  fflush(logfile);
  return;
}

void WriteState(const char *msg, komplex *state)
{
  /* Output one state vector -- with a text message */
  long long i;

  fprintf(outfile, "%s \n", msg);
  for (i = 0; i < Nunique; i++)
  {
#ifdef CSVOUT
    fprintf(outfile, "%g,%g,%lld\n", real(state[i]), imag(state[i]), unique[i]);
#else
    fprintf(outfile, "  %g + i * %g\n", real(state[i]), imag(state[i]), unique[i]);
#endif
  }
  return;
}

void WriteStates(komplex **hamil)
{
  long long i, j;

  fprintf(outfile, "{\n");
  for (i = 0; i < Nuniq_k; i++)
  {
    fprintf(outfile, "( ");
    for (j = 0; j < Nuniq_k; j++)
      fprintf(outfile, "( %lg +i %lg ), ", real(hamil[j + 1][i + 1]),
              imag(hamil[j + 1][i + 1]));
    fprintf(outfile, ") \n");
  }

  fprintf(outfile, "}\n");
  return;
}

void WriteGSdata(double energy, long long *symmetry)
{
  fwrite(&energy, sizeof(double), 1, gscofile);
  fwrite(symmetry, sizeof(long long), Nsym, gscofile);
  return;
}
void WriteGSstate(komplex *state)
{
  fwrite(state, sizeof(komplex), Nunique, gscofile);
  return;
}

void ReadGSdata(double *energy, long long *symmetry, komplex *state)
{
  //  gscoinfile=fopen("./12-0.gs","r");
  fread(energy, sizeof(double), 1, gscoinfile);
  fread(symmetry, sizeof(long long), Nsym, gscoinfile);
  fread(state, sizeof(komplex), Nunique, gscoinfile);
  return;
}

void ReadGSenergy(double *energy, long long *symmetry)
{
  //  gscoinfile=fopen("./12-0.gs","r");
  fread(energy, sizeof(double), 1, gscoinfile);
  fread(symmetry, sizeof(long long), Nsym, gscoinfile);
  return;
}

void WritehmQ(long long *q, struct FLAGS *input_flags)
{
  long long i;

  fprintf(outfile, "!");
  if (input_flags->m_sym)
    fprintf(outfile, " m= %g ; ", double(twom) / 2);
  else
    fprintf(outfile, " h= %g ; ", h);
#ifdef TEST_WRITEHMQ
  printf("\n no of q values in dat files, Nsym=%d\n", Nsym);
#endif
  fprintf(outfile, "q= ( ");
  for (i = 0; i < Nsym; i++)
    fprintf(outfile, "%lld ", q[i]);
  fprintf(outfile, " ) \n");

  return;
}

void WriteResults(long long N, struct FLAGS *input_flags)
{
  /* Output the energies and other observables of the eigenstates */
  long long i, q;
#ifndef WRITE_MAGNETISATION
  // otherwise the energy and magnetisation pairs get mixed up
  Bubblesort(energies, NULL, N);
#endif

  fprintf(outfile, "[ \n");
  for (i = 0; i < N; i++)
  {

    // throw away rounding errors
    if (fabs(energies[i]) < SMALL_NUMBER)
    {
      energies[i] = 0;
    }

    fprintf(outfile, "energy= %9.6g ", energies[i]);
#ifdef FIND_MAG
    // throw away rounding errors
    if (fabs(magnetisation[i]) < SMALL_NUMBER)
    {
      magnetisation[i] = 0;
    }
#ifdef WRITE_MAGNETISATION
    if (input_flags->use_exact_matrix)
    {
      fprintf(outfile, ", mag_z= %9.6g ", magnetisation[i]);
    }
#endif // WRITE_MAGNETISATION
#endif /* FIND_MAG */

    fprintf(outfile, "\n");
  }
  fprintf(outfile, "] \n");
  fflush(outfile);

  return;
}

void WriteCross(long long Nener, long long *symvalue, long long flag, struct FLAGS *input_flags)
{
  /* Output the cross sections of the ground state */

  FILE *crossfile;
  if (flag == 0)
  {
    crossfile = outfilezz;
  } // SZZ
  else if (!input_flags->find_cross_pm)
  {
     if (flag == 1)
    {
      crossfile = outfilexx;
    } // SXX
    else if (flag == 2)
    {
      crossfile = outfileyy;
    } // SYY
  }
  else
  {
    if (flag == 1)
    {
      crossfile = outfilemp;
    } // SMP
    else if (flag == 2)
    {
      crossfile = outfilepm;
    } // SPM
  }
  if (flag < 0 || flag > 2)
  {
    LogMessageChar("\nError in WriteCross flag!\n");
  }
  if (Nener <= 0)
  {
    return;
  }

#ifdef TEST_WRITECROSS
  LogMessageChar("In WriteCross, ");
  if (input_flags->m_sym)
    LogMessageCharInt(", m =", twom / 2);
  LogMessageCharDouble(", gs_energy=", gs_energy);

  LogMessageCharInt(", q_gs = (", q_gs[0]);
  for (long long i = 1; i < Nsym; i++)
    LogMessageCharInt(",", q_gs[i]);
  LogMessageChar(")");
  for (int i = 0; i < Nsym; i++)
  {
    LogMessageCharInt("\n i =", i);
    LogMessageCharInt(", symvalue[i] =", symvalue[i]);
  }
  LogMessageChar("\n");
#endif

#ifdef CSVOUT
  // Print header
  if (input_flags->m_sym)
    fprintf(crossfile, "m,");
  else
    fprintf(crossfile, "h,");

  for (int i = 0; i < Nsym; i++)
    fprintf(crossfile, "q%d,", i);
  fprintf(crossfile, "E,S\n");

  // This doesn't sort the array as below, but this shouldn't be needed as
  // CSV-out is mainly for loading into other software and plotting/looking
  // at data there.
  for (long long i = 0; i < Nener; i++)
  {
    if (input_flags->m_sym)
      fprintf(crossfile, "%g,", double(twom) / 2);
    else
      fprintf(crossfile, "%g,", double(h));
    for (int i = 0; i < Nsym; i++)
      fprintf(crossfile, "%lld,", symvalue[i] - q_gs[i]);

    fprintf(crossfile, "%g,%g\n", energies[i] - gs_energy, cross[i]);
  }
#else

  if (input_flags->m_sym)
  {
    fprintf(crossfile, "[\n m= %g, q= (%lld", double(twom) / 2, symvalue[0] - q_gs[0]);
    // fprintf(crossfile, "[ \n m= %g, q= (%lld",double(twom)/2,symvalue[0]);
  }
  else
  {
    fprintf(crossfile, "[ \n h= %g, q= (%lld", double(h), symvalue[0] - q_gs[0]);
    // fprintf(crossfile, "[ \n h= %g, q= (%lld",double(h),symvalue[0]-q_gs[0]);
  }

  for (long long i = 1; i < Nsym; i++) // print spin-wave q: q'=symvalue[i]-q_gs[i]
  {
    fprintf(crossfile, " %lld", symvalue[i] - q_gs[i]); //
    // fprintf(crossfile, ",%2lld",symvalue[i]); //
  }

  fprintf(crossfile, ")\n length of Szpm(q)|gs> = %g", szqlength);
  fprintf(crossfile, " S(q) = %g\n\n", 2 * PI * szqlength);

#ifdef TEST_WRITECROSS
  for (int i = 0; i < Nener + 1; i++)
  {
    LogMessageCharDouble("\nenergies[i] =", energies[i]);
    LogMessageCharDouble(", cross[i] =", cross[i]);
  }
#endif

  long long Xsections = SortCross(Nener);
  // destroys *energies and *Cross. Maybe not so smart?

  for (long long i = 0; i < Nener; i++)
  {
    // for (long long i = 0;i<Nener;i++) {
    // fprintf (outfilezz," omega = %9.6g Szz = %9.6g\n",energies[i]-gs_energy,cross[i]);
    fprintf(crossfile, " %9.6g %9.6g\n", energies[i] - gs_energy, cross[i]);
    // fprintf (crossfile," %9.6g %9.6g\n",energies[i],cross[i]);
  }
  fprintf(crossfile, "]\n\n");
#endif // CSVOUT
  return;
}

void WriteGSEnergy(komplex E)
{
  printf("<e|H|e> = %g +i %g\n", real(E), imag(E));

  return;
}

void WriteEnergy(double re)
{
  printf(" E = %g \n", re);

  return;
}

#ifdef FIND_MAG
void WriteMaggs(long long *q)
{
  int i;
  WriteQvalue(q);
  fprintf(outfile, "\nMagnetization of groundstate %g\n", maggs);
  fprintf(outfile, "\n");
}
#endif // FIND_MAG

void WriteQvalue(long long *qvec)
{
  int i;
  fprintf(outfile, "q-vector: (");
  fprintf(outfile, " %lld", qvec[0]);
  for (i = 0; i < Nsym; i++)
  {
    fprintf(outfile, ", %lld", qvec[i]);
  }
  fprintf(outfile, ")\n");
}
