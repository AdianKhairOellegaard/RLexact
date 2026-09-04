/* Program file RLcross.C -
 * Contains functions for calculating local spin operator expectation values
 * Last change: AKOE 04.09.2026
 *
 ============================================
 *
 * RLexact: The exact diagonalization package
 * Christian Rischel & Kim Lefmann, 26.02.94
 * Version 4.0, September 2017
 *
 ============================================
 */
#include <complex>
#include <RLexact.h>
#include "Functions.h"

/* Functions defined in this file */

/* Global variables defined in RLexact.c */
extern long long Nspins;

/* Regional variables defined in this file */



/* Local spin operation functions */

//S^z_p operator function
int ApplySz(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        unsigned long long mask = ((unsigned long long)1) << p;
        int spin_value = 2 * ((*bitmap & mask) != 0) - 1;
        return spin_value;
    }
    else
    {
        fatalerror("Position exceeds spin positions", p);
        return 0;
    }
}

//S^+_p operator function
int ApplySp(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        unsigned long long mask = ((unsigned long long)1) << p;

        if ((*bitmap & mask) == 0)
        {
            *bitmap = *bitmap | mask;
            return 1;
        }
        return 0;
    }
    else
    {
        fatalerror("Position exceeds spin positions", p);
        return 0;
    }
}

//S^m_p operator function
int ApplySm(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        unsigned long long mask = ((unsigned long long)1) << p;

        if ((*bitmap & mask) != 0)
        {
            *bitmap = *bitmap & ~mask;
            return 1;
        }
        return 0;
    }
    else
    {
        fatalerror("Position exceeds spin positions", p);
        return 0;
    }
}