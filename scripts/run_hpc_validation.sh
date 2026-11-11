#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
unset MPI_TEST_SINGLE
# Load the site's matching compiler and MPI modules BEFORE invoking this script.
# Run within an allocation with >=4 CPU slots, never as a login-node benchmark.
: "${MPIEXEC:=mpirun}"
: "${MPI_NP_FLAG:=-np}"
export MPIEXEC MPI_NP_FLAG
run_dir="results/hpc_$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir -p "$run_dir"
exec > >(tee "$run_dir/workflow.log") 2>&1
printf 'Output directory: %s\n' "$run_dir"
{
  date -u
  hostname
  uname -srm
  command -v "${CXX:-g++}"
  "${CXX:-g++}" --version
  command -v "${MPICXX:-mpicxx}"
  "${MPICXX:-mpicxx}" --version
  python3 --version
  printf 'MPIEXEC=%s\nMPI_NP_FLAG=%s\n' "$MPIEXEC" "$MPI_NP_FLAG"
  if type module >/dev/null 2>&1; then module list; fi
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    printf 'SLURM_JOB_ID=%s\nSLURM_NTASKS=%s\nSLURM_CPUS_PER_TASK=%s\n' "$SLURM_JOB_ID" "${SLURM_NTASKS:-}" "${SLURM_CPUS_PER_TASK:-}"
  fi
} > "$run_dir/environment.txt" 2>&1
make clean
make CXX="${CXX:-g++}" MPICXX="${MPICXX:-mpicxx}"
make mpi CXX="${CXX:-g++}" MPICXX="${MPICXX:-mpicxx}"
make test CXX="${CXX:-g++}" MPICXX="${MPICXX:-mpicxx}"
cp results/local_validation.csv "$run_dir/serial_validation.csv"
make test-mpi CXX="${CXX:-g++}" MPICXX="${MPICXX:-mpicxx}"
cp results/mpi_field_comparison.csv "$run_dir/mpi_field_comparison.csv"
printf 'PASS: HPC validation complete. Return this directory: %s\n' "$run_dir"
