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
#include <strings.h>
#include <stdio.h>
#include "Functions.h"

/* Functions defined in this file */

/* Global variables defined in RLexact.c */
extern long long Nspins, Nsymops, Nsym, Nsymvalue[NSYM];
extern long long *Nocc;
extern int numcanonical;
extern CanonicalRep *canonical;
extern OrbitTable *otable;
/* Regional variables defined in this file */
unsigned long long maskN;
unsigned long long maskp;
/* Builds symmetry orbits used to simplify local expectation value calculations */
void BuildOrbitTable(OrbitTable *otable, CanonicalRep *canonical, struct FLAGS *input_flags)
{
    long long symcom_count, spin_count, sym = 0;
    long long p_new, alpha_new;
    int T[NSYM];
    komplex lambda;
    unsigned long long state, new_state;
    numcanonical = 0;
    maskN = (((unsigned long long)1) << Nspins) - 1;
    int *single_visit = (int *)calloc(Nspins * 3, sizeof(*single_visit));

    if (single_visit == NULL)
    {
        fatalerror("Could not allocate single_visit", Nspins);
    }

    for (int p = 0; p < Nspins; ++p)
    {
        for (int alpha = 0; alpha < 3; ++alpha)
        {
            // If this state has already been visited by a previous orbit, skip it!
            if (single_visit[p * 3 + alpha]) continue;

            // Otherwise it is a new canonical
            canonical[numcanonical++] = (CanonicalRep){p, alpha};
            
            state = ((unsigned long long)1) << p;
            symcom_count = 0;
            // Now, loop over all symmetry operations to mark the canonicals entire orbit.
            TLOOP_BEGIN
            spin_count = Count(new_state, input_flags);
            // Check if a spin-flip occurred by counting the number of set bits
            if (spin_count == 1)
            {
                p_new = ffsll(new_state) - 1;
                alpha_new = alpha;
                lambda = 1;
            }
            else
            {
                p_new = ffsll((~new_state) & maskN) - 1;
                if (alpha == 2)
                {
                    alpha_new = alpha;
                    lambda = -1;
                }
                else
                {
                    alpha_new = 1 - alpha;
                    lambda = 1;
                }
            }
            //Save symmetrical values and lambda.
            store_orbit_element(otable, p, alpha, symcom_count++, p_new, alpha_new, lambda);
            single_visit[p_new * 3 + alpha_new] = 1;
            TLOOP_END
        }
    }
    free(single_visit);
}

//Function for retrieving array index in the flat orbit table
inline long long get_orbit_idx(int p, int alpha, long long symcom_id)
{
    if (p < 0 || p >= Nspins)
        fatalerror("Spin position index is outside the range of number of spins", (long long)p);
    if (alpha < 0 || alpha >= 3)
        fatalerror("Spin operator index is outisde the range of 3 (+,-,z)", (long long)alpha);
    if (symcom_id < 0 || symcom_id >= Nsymops)
        fatalerror("Spin operator index is outside the range of symmetry operator combinations", Nsymops);  
    return (p * 3 + alpha) * Nsymops + symcom_id;
}
//A storing function for the orbit table
inline void store_orbit_element(OrbitTable* table, int p, int alpha, long long symcom_id, 
                                       int p_prime, int alpha_prime, komplex lambda)
{   
    if (p_prime < 0 || p_prime >= Nspins)
        fatalerror("Spin position index is outside the range of number of spins", (long long)p_prime);
    if (alpha_prime < 0 || alpha_prime >= 3)
        fatalerror("Spin operator index is outisde the range of 3 (+,-,z)", (long long)alpha_prime);
    long long idx = get_orbit_idx(p, alpha, symcom_id);
    table[idx] = (OrbitTable){p_prime, alpha_prime, lambda};
}
//A getter function
inline const OrbitTable* get_orbit_element(const OrbitTable* table, int p, int alpha, long long symcom_id) {
    long long idx = get_orbit_idx(p, alpha, symcom_id);
    return &table[idx];
}

/* Expectation value function */

/* Local spin operation functions */

//S^z_p operator function
int ApplySz(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        maskp = ((unsigned long long)1) << p;
        int spin_value = 2 * ((*bitmap & maskp) != 0) - 1;
        return spin_value;
    }
    else
    {
        fatalerror("Spin position index is outside the range of number of spins", (long long)p);
        return 0;
    }
}

//S^+_p operator function
int ApplySp(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        maskp = ((unsigned long long)1) << p;

        if ((*bitmap & maskp) == 0)
        {
            *bitmap = *bitmap | maskp;
            return 1;
        }
        return 0;
    }
    else
    {
        fatalerror("Spin position index is outside the range of number of spins", (long long)p);
        return 0;
    }
}

//S^m_p operator function
int ApplySm(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        maskp = ((unsigned long long)1) << p;

        if ((*bitmap & maskp) != 0)
        {
            *bitmap = *bitmap & ~maskp;
            return 1;
        }
        return 0;
    }
    else
    {
        fatalerror("Spin position index is outside the range of number of spins", (long long)p);
        return 0;
    }
}