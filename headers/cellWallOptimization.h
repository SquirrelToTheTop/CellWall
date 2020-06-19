#ifndef CELLWALL_OPTIMISATION_H
#define CELLWALL_OPTIMISATION_H

#include "math.h"

#include "cellWallParams.h"

void benchmarks_energy_cw(CellWallMonolayer *cwl, int mpi_rank, int mpi_size);

void analyze_cpu_time_energy(CellWallMonolayer *cwl, CellWallLipidLayer *ll, int);

// Simulated annealing technic or also called "recuit simulé" in sexy french
void optimize_simulated_annealing(CellWallMonolayer *, CellWallLipidLayer *, int);

// Simulated annealing technic or also called "recuit simulé" in sexy french
void optimize_simulated_annealing_force(CellWallMonolayer *, CellWallLipidLayer *, int);

// conjugate gradient
void conjugate_gradient(CellWallMonolayer *, int, int);

// test MC
void MC_simulated_annealing(CellWallMonolayer *cwl, int mpi_rank, int mpi_size);

#endif
