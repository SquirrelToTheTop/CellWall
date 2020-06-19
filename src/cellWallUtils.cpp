#include <stdio.h>

#include "cellWallConfig.h"
#include "cellWallUtils.h"

// Sexy and entertaining welcome message
void welcome_message(){

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