# DSA5208 Group Project 3 — Parallel Poisson Solver

This project contains three C++17 Jacobi implementations for the 2D and 3D Poisson equation: **serial**, **pure MPI**, and **hybrid MPI + OpenMP**. Numerical validation covers both 2D and 3D. The recorded scaling benchmarks cover **2D only**, using exactly 1,000 Jacobi iterations per run. The performance figures therefore measure iteration-loop time, **not** time to numerical convergence.

## Code layout

- `src/serial.cpp`, `include/solver.hpp`: serial Jacobi solver and numerical diagnostics.
- `src/mpi.cpp`, `include/mpi_solver.hpp`: MPI solver with one-dimensional slab decomposition and blocking halo exchange.
- `src/hybrid.cpp`, `include/hybrid_solver.hpp`: MPI + OpenMP version with thread-parallel local updates.
- `include/options.hpp`: CLI parser shared across solvers.
- `Makefile`: builds all three executables and provides test targets for all three implementations.
- `tests/`: stencil, field-comparison, convergence, and CLI tests.
- `scripts/`: HPC validation helper and a Slurm template (the recorded performance jobs used PBS).
- `*.pbs`: PBS performance and validation job scripts.
- `results/`: raw CSV files, logs, and validation records. Keep these as evidence for the report.

## Build

The Atlas8 experiments used GCC 9.2.0, Open MPI 4.1.2, and C++17 with `-O3`; OpenMP builds use `-fopenmp`. Builds require a C++17 compiler, an MPI compiler wrapper and launcher for the parallel targets, and OpenMP support for the hybrid target. The tests use Python 3. Load the matching compiler and MPI modules on Atlas8. From the project directory:

```bash
make clean
make CXX=g++ MPICXX=mpicxx
make mpi CXX=g++ MPICXX=mpicxx
make hybrid CXX=g++ MPICXX=mpicxx
```

The executables are `poisson_serial`, `poisson_mpi`, and `poisson_hybrid`. Compiled binaries are not included; rebuild them on the target system.

## Validation

Run tests inside a suitable multi-process allocation (at least four CPU slots for the MPI tests):

```bash
mkdir -p results
make test
make test-mpi
make test-hybrid
```

The recorded Atlas8 validation job was **1037997**. The completed tests included 14 serial numerical/stopping cases, seven invalid-input cases, 30 MPI full-field comparisons, and 48 hybrid full-field comparisons across rank/thread settings. The manufactured exact solution is the product of sine functions. The observed maximum-error convergence order is approximately two for both 2D and 3D. Small-grid tests do not establish correctness for every possible problem size or rank count.

### Example runs with Open MPI

The hybrid example requests six CPU cores.

```bash
./poisson_serial --dim 2 --n 64 --mode benchmark --iterations 1000 --output results/own_serial.csv
mpirun --bind-to core -np 4 ./poisson_mpi --dim 2 --n 64 --mode benchmark --iterations 1000 --output results/own_mpi.csv
OMP_NUM_THREADS=2 mpirun --map-by slot:PE=2 --bind-to core -np 3 ./poisson_hybrid --dim 2 --n 64 --mode benchmark --iterations 1000 --output results/own_hybrid.csv
```

Create the output directory before running. Validation tests can regenerate CSV files under `results/`. The original Atlas8 validation records are in [results/validation_success](results/validation_success). For convergence tests, use `--mode validation` with appropriate `--tol`, `--max-iters`, and `--check-every` settings.

## Recorded performance jobs

| Job | PBS script | Scope | Host noted in log |
| --- | --- | --- | --- |
| [1038177](results/serial_scaling/serial_1038177.csv) | `serial_scaling.pbs` | Serial N=256,512,1024,2048 | `cnode-33-43-28` |
| [1038009](results/strong_scaling/strong_1038009.csv) | `strong_scaling_12.pbs` | Serial + MPI 1–12 ranks | `cnode-33-43-28` |
| [1038129](results/strong_scaling/final_1038129.csv) | `strong_scaling_final.pbs` | Controlled MPI 1–96 ranks in one 96-CPU job | `cnode-33-43-52` |
| [1038022](results/weak_scaling/weak_1038022.csv) | `weak_scaling_12.pbs` | MPI weak scaling, 1–12 ranks | `cnode-33-43-28` |
| [1038015](results/hybrid_scaling/hybrid_1038015.csv) | `hybrid_scaling_12.pbs` | Six rank/thread mixes totaling twelve CPUs | `cnode-33-43-28` |
| [1038027](results/strong_scaling/strong96_1038027.csv) | `strong_scaling_96.pbs` | Separate 24/48/96-rank check | `cnode-33-43-52` |

Each setting was repeated three times; the report uses medians. The **strong-scaling speedup and efficiency use only job 1038129**, with its own one-rank MPI baseline. Other jobs had different allocations, and at least two nodes were involved. Do not combine runtimes from these jobs to claim a controlled serial-to-MPI or MPI-to-hybrid speedup. In particular, the pure-MPI 12-rank median was **0.5005 s** in job 1038129 but **2.1012 s** in job 1038009. The latter ran on a different node with a 12-CPU rather than 96-CPU allocation.

The original PBS output for these jobs is in [results/logs](results/logs). On Atlas8, build the executables before submitting the corresponding PBS script from the repository root. For example, run `qsub strong_scaling_final.pbs` to repeat the main MPI strong-scaling benchmark.

### Known limitations

- Fixed-iteration timing is **not** time to convergence. At N=2048, the relative residual after 1,000 iterations was approximately 0.9988.
- All performance runs were on one node at a time. No multi-node scaling data were collected.
- Weak scaling was tested up to twelve ranks because the recorded weak-scaling PBS script requested twelve CPUs. This does not explain the slowdown; further profiling is needed.
- The `compute_sec` and `halo_sec` CSV columns are separately reduced maximum rank times; **do not add them** to reconstruct total runtime. Blocking `halo_sec` includes waiting, not only network transfer.
- Separate jobs were not matched for node, resource reservation, affinity, and other system conditions. Reported runtime gaps are observations, not evidence for a specific hardware cause.

## AI usage disclosure

We used ChatGPT to assist with project planning, code and PBS script development and debugging, review of benchmark records, and report editing. The reported runtimes and job IDs come from the archived Atlas8 CSV files and PBS logs. Speedups and efficiencies were calculated from the recorded runtimes.
