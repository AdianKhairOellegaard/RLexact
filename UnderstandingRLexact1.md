# RLexact Study Notes (Conversation Summary)

Date: 2026-07-07
Purpose: Clean summary of today's walkthrough of the RLexact codebase, concepts, and outputs.

## Quick Mental Model
RLexact follows this pipeline:

raw spin bitstates -> symmetry-reduced unique basis -> q-sector projection -> Hamiltonian action (sparse) -> Lanczos ground state -> operator-applied seeds -> spectra S(q,w)

If this chain makes sense, most functions in the project will also make sense.

## 1) What RLexact does
RLexact is an exact-diagonalization code for spin models. In practical terms, it:
- reads model and run settings from input files,
- reduces Hilbert space with symmetries,
- builds Hamiltonian representation (primarily sparse),
- computes energies and eigenvectors (mostly Lanczos),
- computes dynamical cross sections like S(q,w).

## 2) Main runtime flow (src/RLexact.C)
Typical flow:
1. Initialize MPI (`MPI_Init`, rank, process count).
2. Read input and initialize files (`intro`).
3. Build lookup tables (`BuildTables`) and symmetry metadata (`InitSym`).
4. Build/read unique basis (`FillUnique` or `ReadUnique`).
5. Allocate vectors/matrices.
6. Loop over sectors (magnetization sectors or field values).
7. Solve with Lanczos or exact matrix mode.
8. Optionally reconstruct GS and compute cross sections.
9. Finalize output (`outro`), free memory, finalize MPI.

## 3) C/C++ and runtime basics clarified
- `#include`: inserts header content before compile.
- `void`: function has no return value.
- `extern ...`: declaration of symbol defined in another file.
- `#ifndef`: compile code only if macro is not defined.
- `long long`: large integer type for indices/counters/state ids.
- `main(int argc, char *argv[])`: command-line entry point.
- `malloc(...)`: heap allocation.
- buffering: data is first accumulated in memory for efficiency.
- `fflush(file)`: force buffered data to be written now.

## 4) Important flags and modes
### Compile-time flags
- `MOTIVE`: enables unit-cell position aware phase handling and related q logic.
- `FIND_MAG`: magnetization-related path; mostly legacy and typically off.

### Runtime mode flags
- `MODEN`: normal full workflow.
- `MODEGS`: find and store GS metadata.
- `MODERC`: reconstruct GS from stored metadata.
- `MODEQ`: evaluate cross section for chosen q.
- `UNIMODEW`: write unique tables.
- `UNIMODER`: read unique tables.

## 5) MPI behavior that matters
- Lanczos path is MPI-parallel over q sectors.
- Matrix path is effectively rank-0 centric in current code.
- `MPI_Barrier`: wait until all ranks reach same point.
- `MPI_Bcast`: copy root data to all ranks.
- `MPI_Allreduce(..., MPI_MINLOC)`: find global minimum value and the rank that found it.

## 6) Symmetry-reduced basis (core idea)
### `FillUnique`
Builds canonical representatives of symmetry orbits.
- Huge basis-size reduction.
- Can run in count-only mode for memory sizing.
- Returns and/or populates `Nunique` and `unique[]`.

### `ReadUnique` / `WriteUnique`
Reuse precomputed unique-basis tables from disk.

### `FillUniqueObservables`
Prepares per-unique diagonal helper arrays (notably `mag[]` in current paths).

### `BuildCycle(q)`
For one symmetry sector q, computes orbit consistency and normalization data (`Nocc`).
- If projection phases cancel, that projected basis component is zero in that sector.

### `FindUnique`
Maps any state bitmap to canonical representative and records symmetry-translation bookkeeping (`Tvec`).

## 7) Symmetry operations (src/RLsymm.C)
### `SymOp(sym_index, bitmap)`
Applies one symmetry from `symlist`:
- identity,
- spin flip,
- translations,
- FCC32-specific point/space operations,
- user-added operators.

### `InitSym()`
Initializes symmetry periods/ranges used in sector loops.

## 8) Sparse Hamiltonian workflow (src/RLsparse.C)
### `MakeSparse`
Writes sparse Hamiltonian data to binary files (index/coupling/translation/diag/count/total).

### `WriteCouplingFiles`
Maps generated transitions to unique basis and writes one sparse off-diagonal entry.

### `ApplySparse`
Computes matvec `v_out = H * v_in` from sparse files.
- Includes q-dependent phase factors and `Nocc` normalization.
- Main runtime hotspot for Lanczos.

### `FillHamilSparse`
Builds dense matrix from sparse representation for exact matrix mode.

## 9) Hamiltonian term builders (src/RLhamil.C)
### `Hamil_Zeeman`
- z-field contributes diagonal term.
- x/y-field contributes spin-flip off-diagonal terms.

### `Hamil2_sparse`
Pairwise exchange terms:
- diagonal Ising-like contributions,
- off-diagonal XY/anisotropic transitions.

### `Hamil4_sparse`
Ring exchange (four-spin) contributions with pattern-dependent transitions.

## 10) Lanczos workflow
### `LowestLanczos`
Sector-level Lanczos routine:
- finds lowest Ritz value,
- optionally reconstructs eigenvector in RECONSTRUCT mode.

### Ground-state logic in `Solve_Lanczos`
- scan chosen/all q sectors,
- reduce to global minimum with MPI,
- broadcast `q_gs` and GS vector,
- save metadata via `WriteGSdata(gs_energy, q_gs)`.

## 11) Cross-section workflow (src/RLcross.C)
### `CrossLanczos`
Top-level spectral routine for one q:
- build operator-applied seed,
- run Lanczos in CROSS mode,
- write spectral channels.

### `ApplySzq`
Builds seed vector `S^z_q |GS>`.

### `ApplySmp`
Builds `S^+_q |GS>` or `S^-_q |GS>` by:
- local spin operation on bitmaps,
- canonical remapping (`FindUnique` and `LookUpU`),
- symmetry/momentum phase and normalization factors.

### General operator question
Yes, in principle. Any operator can be added if you provide:
- its local action rules,
- transition amplitudes,
- mapping back to unique basis,
- phase/normalization handling consistent with current basis conventions.

## 12) Output format interpretation
### Groundstate lines in `.dat`
From `WriteState`, each line is:

`real_part,imag_part,unique_state_id`

Example (`RLexamples/test/J_4x1x1.h-0.dat`):
- `0.267573,0.511603,3` -> coefficient `0.267573 + i*0.511603` on unique state id 3.
- `-0.378405,-0.723516,5` -> coefficient `-0.378405 - i*0.723516` on unique state id 5.

Note: the third value is a unique-state label, not a translation index.

## 13) RLio utility behaviors (src/RLio.C)
- `time_stamp(start/stop, label)`: block timing logger.
- `outro`: closes all output/cross files safely.
- `LogMessage*`: structured logging helpers.
- `SortCross(Nener)`: merge duplicate-energy bins, sum weights, sort by energy, return compacted count `XS`.
- `Nener`: number of raw energy entries before compaction.

## 14) Final understanding from this session
By the end of today, the full conceptual route became clear:

state representation -> symmetry reduction -> sector projection -> Hamiltonian application -> eigen solving -> operator application -> spectral output.

## Suggested next step
Do one concrete, line-by-line run trace for a small case (for example `RLexamples/test/J_4x1x1.h`) and record:
- sector values,
- generated files,
- where each major function enters and exits.

That is usually the fastest way to make this understanding fully operational.