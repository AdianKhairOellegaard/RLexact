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
#include <cnr.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include "Functions.h"

/* Functions defined in this file */

/* Global variables defined in RLexact.c */
extern long long Nspins, Nsymops, Nsym, Nunique, Nsymvalue[NSYM];
extern long long *Nocc, q_gs[NSYM];
extern unsigned long long *unique;
extern komplex *gs;
extern OrbitTable *otable;
extern komplex *S1exp;
extern komplex **S2exp;
extern double cosine[], sine[], sqroot[];
extern double *one_tangle, **concurrence, *two_tangle, *QFI;
/*Only relevant for QFI calculation using expectation values*/
extern long long q_T_count;
extern long long *TransIds;
extern long long Ndimensions, Nspins_in_uc;
extern double **spin_positions;
extern long long Nsymvalue[NSYM];
extern long long Trans_Qmax[3];

/* Output files */
extern FILE *outfileexp;
extern FILE *outfilewit;

/* Regional variables defined in this file */
unsigned long long maskN, maskp, outvec;

/* Finds 1 spin Expectation values of the system*/
void FindExpectationValues(OrbitTable *otable, CanonicalRep *canonical, int numcanonical, struct FLAGS *input_flags)
{
    komplex expectation_value;
    long long p, alpha;
    const OrbitTable *el;
    komplex lambda;
    komplex *dummyS;
    double (*op)(int, unsigned long long *) = NULL;
    int m_zero_exp;

    for (int i = 0; i < numcanonical; i++)
    {
        p = canonical[i].p;
        alpha = canonical[i].alpha;

        switch(alpha)
        {
            case (PLUS):
                op = ApplySp;
                m_zero_exp = 1;
                break;
            case (MINUS):
                op = ApplySm;
                m_zero_exp = 1;
                break;
            case (Z):
                op = ApplySz;
                m_zero_exp = 0;
                break;
            default:
                fatalerror("The canonical does not have a spin value in range", alpha);
        }

        if (input_flags->m_sym && m_zero_exp != 0)
            expectation_value = 0;
        else
            expectation_value = expect_value(op, p, input_flags);

        for (long long symcom_count = 0; symcom_count < Nsymops; symcom_count++)
        {
            el = get_orbit_element(otable, p, alpha, symcom_count);
            S1exp[el->target_site*3 + el->target_alpha] = el->lambda * expectation_value;
        }
    }
    
    /*Turn from S^+, S^- to S^x, S^y unless flagged otherwise*/
    if ((!input_flags->find_expect_pm || input_flags->find_witness_exp) && (!input_flags->m_sym))
    {   
        dummyS = kvector(0, Nspins * 3 - 1);
        memcpy(dummyS, S1exp, Nspins * 3 * sizeof(komplex));
        for (int i = 0; i < Nspins; i++)
        {
            S1exp[i*3 + 0] = (dummyS[i*3 + 0] + dummyS[i*3 + 1]) / 2.0;
            S1exp[i*3 + 1] = (dummyS[i*3 + 0] - dummyS[i*3 + 1]) / (2.0 * I);
        }
        freekvector(dummyS, 0, Nspins * 3 - 1);
    }
    
    if ((!input_flags->find_expect_pm || input_flags->find_witness_exp) && (input_flags->m_sym))
    {   
        for (int i = 0; i < Nspins; i++)
        {
            S1exp[i*3 + 0] = S1exp[i*3 + 2];
            S1exp[i*3 + 1] = S1exp[i*3 + 2];
        }
    }   
}

/* Finds 2 spin expectation values of the system*/
void FindExpectationValues2(OrbitTable *otable, CanonicalPair *canonical2, int numcanonical2, struct FLAGS *input_flags)
{
    komplex expectation_value;
    long long p1, alpha1, p2, alpha2;
    const OrbitTable *el1, *el2;
    double (*op1)(int, unsigned long long *) = NULL;
    double (*op2)(int, unsigned long long *) = NULL;
    komplex *dummyS;
    int m_zero_exp; //Ensures that expecvalue is called minimally if m-symmetry is present.
    int exp_wit; //Ensures that expectvalue is called minimally if only entanglement witnesses are required.
    int no_exp_calc;

    for (int i = 0; i < numcanonical2; i++)
    {
        p1 = canonical2[i].first.p;
        alpha1 = canonical2[i].first.alpha;
        p2 = canonical2[i].second.p;
        alpha2 = canonical2[i].second.alpha;

        switch (3 * alpha1 + alpha2)
        {
            case 3 * PLUS + PLUS:
                op1 = ApplySp;
                op2 = ApplySp;
                m_zero_exp = 1;
                exp_wit = 0;
                break;
            case 3 * PLUS + MINUS:
                op1 = ApplySp;
                op2 = ApplySm;
                m_zero_exp = 0;
                exp_wit = 0;
                break;
            case 3 * PLUS + Z:
                op1 = ApplySp;
                op2 = ApplySz;
                m_zero_exp = 1;
                exp_wit = 1;
                break;
            case 3 * MINUS + PLUS:
                op1 = ApplySm;
                op2 = ApplySp;
                m_zero_exp = 0;
                exp_wit = 0;
                break;
            case 3 * MINUS + MINUS:
                op1 = ApplySm;
                op2 = ApplySm;
                m_zero_exp = 1;
                exp_wit = 0;
                break;
            case 3 * MINUS + Z:
                op1 = ApplySm;
                op2 = ApplySz;
                m_zero_exp = 1;
                exp_wit = 1;
                break;
            case 3 * Z + PLUS:
                op1 = ApplySz;
                op2 = ApplySp;
                m_zero_exp = 1;
                exp_wit = 1;
                break;
            case 3 * Z + MINUS:
                op1 = ApplySz;
                op2 = ApplySm;
                m_zero_exp = 1;
                exp_wit = 1;
                break;
            case 3 * Z + Z:
                op1 = ApplySz;
                op2 = ApplySz;
                m_zero_exp = 0;
                exp_wit = 0;
                break;
            default:
                fatalerror("The canonical pair does not have valid spin values", 3 * alpha1 + alpha2);
        }
        //Expectation value calculation logic
        no_exp_calc = (input_flags->find_witness_exp && !input_flags->find_expect && exp_wit != 0);
        if (input_flags->m_sym)
        {
            if (m_zero_exp != 0)
                expectation_value = 0;
            else
                expectation_value = expect_value2(op1, p1, op2, p2, input_flags);
        }
        else
        {
            if (no_exp_calc)
                expectation_value = 0; //avoids unecessary expectation value calculations if only the entanglement witnesses are required. Supposes that no canonicals orbit from z,(+,-) to (+-),(-+) which is valid.
            else
                expectation_value = expect_value2(op1, p1, op2, p2, input_flags);
        }
        for (long long symcom_count = 0; symcom_count < Nsymops; symcom_count++)
        {
            el1 = get_orbit_element(otable, p1, alpha1, symcom_count);
            el2 = get_orbit_element(otable, p2, alpha2, symcom_count);
            S2exp[el1->target_site * 3 + el1->target_alpha][el2->target_site * 3 + el2->target_alpha] = el1->lambda * el2->lambda * expectation_value;
        }
    }

    /*Turn from S^+, S^- to S^x, S^y unless flagged otherwise*/
    if (!input_flags->find_expect_pm || input_flags->find_witness_exp)
    {   
        dummyS = kvector(0, 8);
        for (int i = 0; i < Nspins; i++)
        {
            for (int j = 0; j < Nspins; j++)
            {
                for (int a1 = 0; a1 < 3; a1++)
                    for (int a2 = 0; a2 < 3; a2++)
                        dummyS[a1*3 + a2] = S2exp[i*3 + a1][j*3 + a2];
                //S^xS^x
                S2exp[i*3 + 0][j*3 + 0] = (dummyS[0*3 + 0] + dummyS[0*3 + 1] + dummyS[1*3 + 0] + dummyS[1*3 + 1]) / 4.0;
                //S^xS^y
                S2exp[i*3 + 0][j*3 + 1] = (dummyS[0*3 + 0] - dummyS[0*3 + 1] + dummyS[1*3 + 0] - dummyS[1*3 + 1]) / (4.0 * I);
                //S^xS^z
                S2exp[i*3 + 0][j*3 + 2] = (dummyS[0*3 + 2] + dummyS[1*3 + 2]) / 2.0;
                //S^yS^x
                S2exp[i*3 + 1][j*3 + 0] = (dummyS[0*3 + 0] + dummyS[0*3 + 1] - dummyS[1*3 + 0] - dummyS[1*3 + 1]) / (4.0 * I);
                //S^yS^y
                S2exp[i*3 + 1][j*3 + 1] = - (dummyS[0*3 + 0] - dummyS[0*3 + 1] - dummyS[1*3 + 0] + dummyS[1*3 + 1]) / 4.0;
                //S^yS^z
                S2exp[i*3 + 1][j*3 + 2] = (dummyS[0*3 + 2] - dummyS[1*3 + 2]) / (2.0 * I);
                //S^zS^x
                S2exp[i*3 + 2][j*3 + 0] = (dummyS[2*3 + 0] + dummyS[2*3 + 1]) / 2.0;
                //S^zS^y
                S2exp[i*3 + 2][j*3 + 1] = (dummyS[2*3 + 0] - dummyS[2*3 + 1]) / (2.0 * I);
                //S^zS^z no change
            }
        }
        freekvector(dummyS, 0, 8);
    }   

}


/* Builds symmetry orbits used to simplify local expectation value calculations, set to X, Y and Z but can easily be changed to +,-,Z if required  */
void BuildOrbitTable(OrbitTable *otable, CanonicalRep *canonical, int *numcanonical, struct FLAGS *input_flags)
{
    long long symcom_count, spin_count, sym = 0;
    long long p_new, alpha_new;
    int T[NSYM];
    komplex lambda;
    unsigned long long state, new_state;
    *numcanonical = 0;
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
            // Keep canonical representatives unique, but populate orbit entries for every operator.
            if (!single_visit[p * 3 + alpha])
                canonical[(*numcanonical)++] = (CanonicalRep){p, alpha};
            
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

void BuildCanonicalPairList(CanonicalPair *canonical2, int *numcanonical2, struct FLAGS *input_flags)
{
    const long long operator_count = Nspins * 3;
    const long long pair_count = operator_count * operator_count;
    bool *visited = (bool *)calloc(pair_count, sizeof(bool));

    *numcanonical2 = 0;
    if (visited == NULL)
        fatalerror("Could not allocate pair-orbit visit table", pair_count);

    for (int p1 = 0; p1 < Nspins; ++p1)
    {
        for (int alpha1 = 0; alpha1 < 3; ++alpha1)
        {
            for (int p2 = 0; p2 < Nspins; ++p2)
            {
                for (int alpha2 = 0; alpha2 < 3; ++alpha2)
                {
                    long long pair_idx = (p1 * 3 + alpha1) * operator_count
                                       + p2 * 3 + alpha2;
                    if (visited[pair_idx])
                        continue;

                    canonical2[(*numcanonical2)++] = (CanonicalPair){{p1, alpha1}, {p2, alpha2}};

                    for (long long symcom_count = 0; symcom_count < Nsymops; symcom_count++)
                    {
                        const OrbitTable *first = get_orbit_element(otable, p1, alpha1, symcom_count);
                        const OrbitTable *second = get_orbit_element(otable, p2, alpha2, symcom_count);
                        long long target_idx =
                        (first->target_site * 3 + first->target_alpha) * operator_count
                        + second->target_site * 3 + second->target_alpha;
                        //if (first->target_site < 0 || first->target_site >= Nspins ||
                        //    first->target_alpha < 0 || first->target_alpha >= 3 ||
                        //    second->target_site < 0 || second->target_site >= Nspins ||
                        //    second->target_alpha < 0 || second->target_alpha >= 3)
                        //{
                        //    fatalerror("Invalid transformed pair in orbit table", target_idx);
                        //}
                        visited[target_idx] = true;
                    }
                }
            }
        }
    }
    free(visited);
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

/* Expectation value function for any local operator that mutates bitmaps cleanly <gs| Op_p |gs>*/
komplex expect_value(double (*Op)(int, unsigned long long*), int p, struct FLAGS *input_flags)
{
    long long symcom_count, sym = 0;
    long long l;
    double operation_value;
    unsigned long long state, new_state, op_state, u;
    int T[NSYM], d_T[NSYM];
    komplex base_fac, fac, sum;
    double phase;

    sum = 0.0;
    for (long long i = 0; i < Nunique; ++i)
    {
        if (abs(gs[i]) == 0)//(abs(gs[i]) < SMALL_NUMBER)
            continue;
        
        state = unique[i];
        base_fac = gs[i] / (sqroot[Nocc[i]] * Nsymops);
        fac = 0;
        TLOOP_BEGIN
        op_state = new_state;
        operation_value = Op(p, &op_state);
        if (operation_value != 0)
        {
            u = FindUnique(op_state, d_T, input_flags);
            l = LookUpU(u, input_flags);

            if (abs(gs[l]) != 0)//(abs(gs[l]) >= SMALL_NUMBER)
            {
                phase = 0.0;
                for (int j = 0; j < Nsym; ++j)
                {
                    phase += 1.0 * q_gs[j] * (d_T[j] - T[j]) / Nsymvalue[j];
                }
                fac += ((komplex)operation_value) * conj(gs[l]) * sqroot[Nocc[l]] * exp(2.0 * I * PI * phase);
            }
        }
        TLOOP_END
        sum += base_fac*fac;
    }
    return sum;
}

/* Expectation value of a product of two local operators: <gs| Op1_p1 Op2_p2 |gs> */
komplex expect_value2(double (*Op1)(int, unsigned long long*), int p1,
                      double (*Op2)(int, unsigned long long*), int p2,
                      struct FLAGS *input_flags)
{
    long long symcom_count, sym = 0;
    long long l;
    double val1, val2, operation_value;
    unsigned long long state, new_state, op_state, u;
    int T[NSYM], d_T[NSYM];
    komplex base_fac, fac, sum;
    double phase;

    sum = 0.0;
    for (long long i = 0; i < Nunique; ++i)
    {
        if (abs(gs[i]) == 0)
            continue;

        state = unique[i];
        base_fac = gs[i] / (sqroot[Nocc[i]] * Nsymops);
        fac = 0;
        TLOOP_BEGIN
        op_state = new_state; // Ops mutate their bitmap argument; don't corrupt the TLOOP state
        val2 = Op2(p2, &op_state);
        if (val2 != 0)
        {
            val1 = Op1(p1, &op_state);
            operation_value = val1 * val2;
            if (operation_value != 0)
            {
                u = FindUnique(op_state, d_T, input_flags);
                l = LookUpU(u, input_flags);

                if (abs(gs[l]) != 0)
                {
                    phase = 0.0;
                    for (int j = 0; j < Nsym; ++j)
                    {
                        phase += 1.0 * q_gs[j] * (d_T[j] - T[j]) / Nsymvalue[j];
                    }
                    fac += ((komplex)operation_value) * conj(gs[l]) * sqroot[Nocc[l]] * exp(2.0 * I * PI * phase);
                }
            }
        }
        TLOOP_END
        sum += base_fac*fac;
    }
    return sum;
}

/* Local spin operation functions */

//S^z_p operator function
double ApplySz(int p, unsigned long long* bitmap)
{
    if (p >= 0 && p < Nspins)
    {
        maskp = ((unsigned long long)1) << p;
        double spin_value = (2 * ((*bitmap & maskp) != 0) - 1) / 2.0;
        return spin_value;
    }
    else
    {
        fatalerror("Spin position index is outside the range of number of spins", (long long)p);
        return 0;
    }
}

//S^+_p operator function
double ApplySp(int p, unsigned long long* bitmap)
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
double ApplySm(int p, unsigned long long* bitmap)
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

//Entanglement witness functions based on expectation values.
void one_tangle_exp(komplex *S1)
{
    for (int p = 0; p < Nspins; p++)
    { 
        if ((imag(S1[p*3 + 0]) > SMALL_NUMBER) || (imag(S1[p*3 + 1]) > SMALL_NUMBER) || (imag(S1[p*3 + 2]) > SMALL_NUMBER))
            fatalerror("x,y and/or z expectation values have a significant imaginary component", p);
        
        one_tangle[p] = 1.0 - 4.0*(sqrabs(S1[p*3 + 0]) + sqrabs(S1[p*3 + 1]) + sqrabs(S1[p*3 + 2]));
    }
}

void concurrence_exp(komplex *S1, komplex **S2)
{
    double res1;
    double res2;
    for (int i = 0; i < Nspins; i++)
    {
        for (int j = 0; j < Nspins; j++)
        { 
            res1 = abs(S2[i*3 + 0][j*3 + 0] - S2[i*3 + 1][j*3 + 1]) - (1.0 / 4.0) + real(S2[i*3 + 2][j*3 + 2]);
            res2 = abs(S2[i*3 + 0][j*3 + 0] + S2[i*3 + 1][j*3 + 1]) - sqrt(sqrabs((1.0 / 4.0) + S2[i*3 + 2][j*3 + 2]) - sqrabs((S1[i*3 + 2] + S1[j*3 + 2]) / 2.0));

            if ((res1 < 0) && (res2 < 0))
                concurrence[i][j] = 0.0;
            else if (res1 < res2)
                concurrence[i][j] = 2.0 * res2;
            else
                concurrence[i][j] = 2.0 * res1;
        }
    }
}

void two_tangle_exp(double **C)
{
    for (int i = 0; i < Nspins; i++)
    {
        two_tangle[i] = 0;
        for (int j = 0; j < Nspins; j++)
        { 
            if (i != j)
                two_tangle[i] += C[i][j] * C[i][j];
        }
    }
}

// This function is really only for comparison and is a slow way of getting QFI. 
// Instead use the dynamical correlation functions. But if find_cross is not used and find_expect is used, then this function is handy.
void QFI_exp(komplex *S1, komplex **S2)
{

    if (S1 == NULL || S2 == NULL || QFI == NULL)
        fatalerror("Invalid arguments passed to QFI_exp", 0);

    long long q[NSYM];
    long long nqvalue[NSYM];
    for (int symmetry = 0; symmetry < Nsym; ++symmetry)
    {
        q[symmetry] = 0;
        nqvalue[symmetry] = 1;
    }
    for (int dimension = 0; dimension < Ndimensions; ++dimension)
    {
        const int symmetry = TransIds[dimension];
        nqvalue[symmetry] = Nsymvalue[symmetry] * Trans_Qmax[dimension];
    }

    komplex *phase_factor = (komplex *)malloc(Nspins * sizeof(komplex));
    if (phase_factor == NULL)
        fatalerror("Could not allocate QFI phase factors", Nspins);

    for (long long q_index = 0; q_index < q_T_count; ++q_index)
    {
        bool zero_q = true;
        for (int dimension = 0; dimension < Ndimensions; ++dimension)
            if (q[TransIds[dimension]] != 0)
                zero_q = false;

        long long spin_index = 0;
        for (long long x = 0; x < Nsymvalue[TransIds[X]]; ++x)
        {
            for (long long y = 0; y < Nsymvalue[TransIds[Y]]; ++y)
            {
                for (long long z = 0; z < Nsymvalue[TransIds[Z]]; ++z)
                {
                    for (long long unit_cell_spin = 0;
                         unit_cell_spin < Nspins_in_uc;
                         ++unit_cell_spin)
                    {
                        double phase = 0.0;
                        const long long cell[3] = {x, y, z};
                        for (int dimension = 0; dimension < Ndimensions; ++dimension)
                        {
                            const long long symmetry = TransIds[dimension];
                            phase += q[symmetry] *
                                     (cell[dimension] +
                                      spin_positions[unit_cell_spin][dimension]) /
                                     Nsymvalue[symmetry];
                        }
                        phase_factor[spin_index++] = exp(I * 2.0 * PI * phase);
                    }
                }
            }
        }

        komplex sq[3] = {zero, zero, zero};
        if (zero_q)
        {
            for (long long i = 0; i < Nspins; ++i)
                for (int alpha = 0; alpha < 3; ++alpha)
                    sq[alpha] += S1[i * 3 + alpha];
        }

        komplex sq_corr[3] = {zero, zero, zero};
        for (long long i = 0; i < Nspins; ++i)
        {
            for (long long j = 0; j < Nspins; ++j)
            {
                const komplex phase = conj(phase_factor[i]) * phase_factor[j];
                for (int alpha = 0; alpha < 3; ++alpha)
                    sq_corr[alpha] += phase * S2[i * 3 + alpha][j * 3 + alpha];
            }
        }

        for (int alpha = 0; alpha < 3; ++alpha)
        {
            sq[alpha] /= sqrt((double)Nspins);
            sq_corr[alpha] /= (double)Nspins;
            QFI[q_index * 3 + alpha] =
                4.0 * (real(sq_corr[alpha]) - sqrabs(sq[alpha]));
        }

        for (int dimension = 0; dimension < Ndimensions; ++dimension)
        {
            const int symmetry = TransIds[dimension];
            if (++q[symmetry] < nqvalue[symmetry])
                break;
            q[symmetry] = 0;
        }
    }

    free(phase_factor);
}