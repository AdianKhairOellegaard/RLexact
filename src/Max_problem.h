// Written by DLC, Jan 2026.
#ifndef MAX_PROBLEM
#define MAX_PROBLEM
// The following parameters define the maximal problem size that is
// allowed when compiling the program. This is for optimization reasons.

#define NCOUP 400          // Max number of couplings
#define NCOUPSTR 10        // Maximal number of different coupling strengths
#define NSPINS 100         // Maximum number of spins
#define NSYM 20            // Maximum number of symmetries
#define NSYMADD 10         // Maximum needs a description
#define NUNIQUE 4500       // Needs a descriotion
#define MAXARRAYSIZE 200   // max number of entries on a given line in the input file
#define NRING 100          // Needs a description
#define NRINGSTR 10        // Maximum number of different ring strengths
#define BUFFERSIZE 4194304 // Set up size of buffer for sparse matrix read/write

#endif