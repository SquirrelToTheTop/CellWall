#include <stdio.h>

#include "cellWallConfig.h"
#include "cellWallUtils.h"

#include <mpi.h>
#include <omp.h>

// Sexy and entertaining welcome message
void welcome_message(){

  int mpi_size;
  int omp_size, omp_max;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  #pragma omp parallel
  {
    omp_size = omp_get_num_threads();
  }

  omp_max = omp_get_max_threads();

  fprintf(stdout, "\n    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    \n");
  fprintf(stdout, "   @----------------------------------@   \n");
  fprintf(stdout, "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@  \n");
  fprintf(stdout, " @---@                              @---@  \n");
  fprintf(stdout, "@---@  CW Bacteria Simulator %d.%d    @---@ \n", 
          cellWall_VERSION_MAJOR, cellWall_VERSION_MINOR);
  fprintf(stdout, "@---@       Gram + / Gram -          @---@\n");                           
  fprintf(stdout, " @---@                              @---@  \n");
  fprintf(stdout,  "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@   \n");
  fprintf(stdout,  "   @----------------------------------@    \n");
  fprintf(stdout,  "    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     \n");

  fprintf(stdout, "\nParallel information : \n");
  fprintf(stdout, "\t> Number of MPI process involved : %d\n", mpi_size);
  fprintf(stdout, "\t> Number of OMP threads involved : %d / %d\n", omp_size, omp_max);

}

/*
 * Calcul la norme d'un vecteur stocké au format x1,y1,z1, ...
 */
double norm_vect(double *arr, int n){

  int i;
  double norm_f;

  norm_f = 0.0f;
  for(i=0; i<n; ++i){
    norm_f += arr[i]*arr[i];
  }

  norm_f = sqrt(norm_f);

  return norm_f;

}