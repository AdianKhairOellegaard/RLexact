// File for calculating neutron scattering intensity for Azurite
// as used by Asbjørn B. Preuss in his thesis (2025).

Ritz_conv 0.000000001 //Lanczos precision
random_number_seed 42 //Number that seeds the random number generator. Should be decommented in real life!
Mode 0 //Perl magic for old cluster routines
Unimode 0
Number of spins 18
Number of spins in unit cell 3
// Coordinates of each spin in the unit cell 
Relative position 0 0 0
Relative position 0.2516 0.49776 0.0835
Relative position -0.2516 0.50224 -0.0835

Number of couplings 36
Number of coupling strengths 4
Coupling strength vector 0.47 0.47 0 //J1
Coupling strength vector 1 1 0 //J2
Coupling strength vector 0.209 0.209 0 //J3
Coupling strength vector 0.139 0.139 0 //Jm
//     10   7   4   1          |x
//      |\ /|\ /|\ /|\    y____|
//  ... |-9-|-6-|-3-|-0
//      |/ \|/ \|/ \|/
//     11   8   5   2
Coupling vector 0 1 0
Coupling vector 0 2 2
Coupling vector 1 2 1
Coupling vector 1 3 2
Coupling vector 2 3 0
Coupling vector 0 3 3
Coupling vector 3 4 0
Coupling vector 3 5 2
Coupling vector 4 5 1
Coupling vector 4 6 2
Coupling vector 5 6 0
Coupling vector 3 6 3
Coupling vector 6 7 0
Coupling vector 6 8 2
Coupling vector 7 8 1
Coupling vector 7 9 2
Coupling vector 8 9 0
Coupling vector 6 9 3
Coupling vector 9 10 0
Coupling vector 9 11 2
Coupling vector 10 11 1
Coupling vector 10 12 2
Coupling vector 11 12 0
Coupling vector 9 12 3
Coupling vector 12 13 0
Coupling vector 12 14 2
Coupling vector 13 14 1
Coupling vector 13 15 2
Coupling vector 14 15 0
Coupling vector 12 15 3
Coupling vector 15 16 0
Coupling vector 15 17 2
Coupling vector 16 17 1
Coupling vector 16 0 2
Coupling vector 17 0 0
Coupling vector 15 0 3
Number of hardcoded symmetries 0
Number of custom symmetries 3
Custom symmetry 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17//Identity
Custom symmetry 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17//TransX
Custom symmetry 15 16 17 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14//TransY
Translation indices 1 2 //Which symmetries are translations?
Qmax translation 1 2 1 //How many extra translations in each direction?
Number of dimensions 2
Number of Chosen GS q-values 0 //0 means search for all
//Chosen GS q-value 0 1 1
Construct symmetries 0 //Leftover from MakeSymCoup(). Needed for program to run
                       //Ask Kim if should be removed from RLio.C/RLsymm.C/...
M start 0
M end 3
