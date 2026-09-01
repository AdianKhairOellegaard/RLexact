// File for calculating neutron scattering intensity for alternating
// J1-J2 chain as used by Asbjørn B. Preuss in his thesis (2025).

Ritz_conv 0.000000001 //Lanczos precision
random_number_seed 42
Mode 0 //Perl magic for old cluster routines
Unimode 0
Number of spins 24
Number of spins in unit cell 2
// Coordinates of each spin in the unit cell 
Relative position 0.0958 0 0
Relative position 0.4057 0 0

Number of couplings 24
Number of coupling strengths 2
Coupling strength vector 1 1 0
Coupling strength vector 0.773 0.773 0
//   23=22-21=20-...
Coupling vector 0 1 0
Coupling vector 1 2 1
Coupling vector 2 3 0
Coupling vector 3 4 1
Coupling vector 4 5 0
Coupling vector 5 6 1
Coupling vector 6 7 0
Coupling vector 7 8 1
Coupling vector 8 9 0
Coupling vector 9 10 1
Coupling vector 10 11 0
Coupling vector 11 12 1
Coupling vector 12 13 0
Coupling vector 13 14 1
Coupling vector 14 15 0
Coupling vector 15 16 1
Coupling vector 16 17 0
Coupling vector 17 18 1
Coupling vector 18 19 0
Coupling vector 19 20 1
Coupling vector 20 21 0
Coupling vector 21 22 1
Coupling vector 22 23 0
Coupling vector 23 0 1
Number of hardcoded symmetries 0
Number of custom symmetries 2
Custom symmetry 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23//Identity
Custom symmetry 22 23 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21//TransX
Translation indices 1 //Which symmetries are translations?
Qmax translation 4 1 1 //How many extra translations in each direction?
Number of dimensions 1
Number of Chosen GS q-values 0 //0 means search for all
//Chosen GS q-value 0 1 1
Construct symmetries 0 //Leftover from MakeSymCoup(). Needed for program to run
                       //Ask Kim if should be removed from RLio.C/RLsymm.C/...
M start 0
M end 0
