
#include "RLexact.h"

#ifndef FUNCTIONS_H // Make sure to only include this once
#define FUNCTIONS_H

// =============================================================================
// Functions in RLexact.C
// =============================================================================
void Solve_Lanczos(struct FLAGS *);
void Solve_Matrix(struct FLAGS *);
void allocate(struct FLAGS *);
void deallocate(struct FLAGS *);

// =============================================================================
// Functions in RLio
// =============================================================================
void WritehmQ(long long *, struct FLAGS *);
void WriteResults(long long, struct FLAGS *);
void WriteCross(long long, long long *, long long, struct FLAGS *);
void WriteMaggs(long long *);
void ReadInputFlags(char *, struct FLAGS *);
long long intro(struct FLAGS *);
long long ReadCoupPattern(char *, struct FLAGS *);
void fatalerror(const char *str, long long i);
void time_stamp(time_t *tim, long long flag, const char *string);
void outro(struct FLAGS *);
void TransformCoup(long long);
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

// =============================================================================
// Functions in regc
// =============================================================================

double atod(char *);
long long regexperr(long long, const char *, bool);
long long multimatch(char *, long long, const char *, double **, long long *, long long, struct FLAGS *);
long long multimatch(char *, long long, const char *, long long **, long long *, long long, struct FLAGS *);
long long matchlines_wrapper(char *, const char *, long long *, bool, struct FLAGS *);
long long matchlines(char *, const char *, double *, bool, struct FLAGS *);

void filereader(char *, char *, long long, struct FLAGS *);
long long filesizer(char *);

// =============================================================================
// Functions in RLcross
// =============================================================================
void CrossLanczos(long long *, struct FLAGS *);
void ApplySzq(long long *, struct FLAGS *);
double lengthofvector(komplex *);
void ApplySmp(long long *, long long, komplex *, struct FLAGS *);

// =============================================================================
// Functions in RLsparse
// =============================================================================
void MakeSparse(struct FLAGS *);
void ApplySparse(komplex *vectin, komplex *vectout, long long *k, struct FLAGS *);
void FillHamilSparse(komplex **hamil, long long *k, struct FLAGS *);
void WriteCouplingFiles(unsigned long long, unsigned long long, int *,
                        komplex, long long, long long, int *, long long *,
                        FILE *, FILE *, FILE *, struct FLAGS *);

// =============================================================================
// Functions in RLlancz
// =============================================================================
double NextLanczos(komplex *, komplex *, komplex *,
                   unsigned long long, long long *, struct FLAGS *);
long long LanczosLoop(long long, long long *, komplex *, struct FLAGS *);
double LowestLanczos(long long *, komplex *, long long *, long long, struct FLAGS *);
double findmag(komplex *, struct FLAGS *);
void findmaggs();
void MakeSeed(komplex *, struct FLAGS *);
void MakeSeedCross(komplex *, long long, struct FLAGS *);

// =============================================================================
// Functions in RLmatrix
// =============================================================================
double Matrix_gs(komplex **, long long *, long long *, komplex *, struct FLAGS *);
void CalculateMatrixM(komplex **, double *, struct FLAGS *);
double CalculateM(komplex *, struct FLAGS *);

// =============================================================================
// Functions in RLhamil
// =============================================================================
void Hamil2_sparse(unsigned long long, unsigned long long *, long long, long long,
                   int *, long long *, int *, komplex *,
                   double *, FILE *, FILE *, FILE *, struct FLAGS *);
void Hamil_Zeeman(unsigned long long, unsigned long long *, long long,
                  int *, long long *, int *, komplex *, double *,
                  FILE *, FILE *, FILE *, struct FLAGS *);
void Eigenvector_test(long long *, komplex *, komplex *, struct FLAGS *);

void Hamil4_sparse(unsigned long long, unsigned long long *,
                   long long, long long, int *, long long *,
                   int *, komplex *, double *, FILE *, FILE *, FILE *,
                   struct FLAGS *);
// =============================================================================
// Functions in RLhamil.dipole
// =============================================================================
void FillHamilton(int k[], komplex **, struct FLAGS *);
double HamDiag(struct FLAGS *);
void Hamil2(int *, komplex, komplex *, struct FLAGS *);
void Hamilton(komplex *, komplex *, int k[], struct FLAGS *);
void matrixelement(komplex, int *, komplex, komplex *, struct FLAGS*);

// =============================================================================
// Functions in Diagonal
// =============================================================================
void htred2(komplex **, long long, komplex *, komplex *, struct FLAGS *);
long long htqli(double *, double *, long long, komplex **);
void Diagonalize(komplex **, long long, double *, struct FLAGS *);

// =============================================================================
// Functions in RLtables
// =============================================================================
void BuildTables(struct FLAGS*);
// Fill the tables of often used math functions and complex phases
long long Count(unsigned long long, struct FLAGS*);
// Count the number of up-spins in a state (represented by a bitmap)
unsigned long long FillUnique(long long, int, struct FLAGS *);
// Fill the table of unique states
void FillUniqueObservables(struct FLAGS *);
// Write diagonal values of the unique states to file
void BuildCycle(long long *, struct FLAGS *);
// Make table of number of occurences of a particular unique in a particular symmetry cycle
unsigned long long FindUnique(unsigned long long, int *, struct FLAGS*);
// Find the unique corresponding to a particular state
long long IsUnique(unsigned long long, struct FLAGS*);
// Test if a given state is a unique
long long LookUpU(unsigned long long, struct FLAGS*);
// Find the index of a given unique state
void InvertMatrix(long long, double[4][4], double[4][4], struct FLAGS*);
// Inverts 1x1 and 2x2 matrices. TODO: test when this is used and if it should be generalized
void WriteUnique(long long, struct FLAGS *);
// Write list of uniques to file
long long ReadUnique(long long, int, struct FLAGS *);
// Read list of uniques from file
void WriteUniqueObservables(struct FLAGS *);
// Write the diagonal elements of the uniques to file
void ReadUniqueObservables(struct FLAGS *);
// Read the diagonal elements of the uniques from file

// =============================================================================
// Functions in RLsymm
// =============================================================================

// Construct the Hamiltonian from the symmetry. TODO: Finish this function
unsigned long long SymOp(long long, unsigned long long); // Perform the actual symmetry operation on the states
unsigned long long SymOpSite(long long, long long);      // Similar to SymOp. TODO: Merge the two functions
void MakeSymCoup(struct FLAGS *);                        // TODO: Beskriv
void TestSym(struct FLAGS *);                            // Test symmetry operations (during initialization)
void InitSym(struct FLAGS *);                            // Initialize the symmetry operations and corresponding tables
// =============================================================================
// Functions in RLutil
// =============================================================================

void reverse(char *);
void itoa(long long, char *);
void FillRotationMatrix(double *, struct FLAGS*);
void NormalizeVector(double *);
void RotateVector(double *);
void Bubblesort(double *, double *, long long);

#endif