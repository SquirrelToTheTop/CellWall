
#include "cellWallParams.h"

#include "cellWallForces.h"

/*
 * Compute the total energy of glycosidic springs according to
 * the equation of spring : e = 0.5 * k_g * (||r2-r1||-d_g)^2
 * (euclidian norme used)
 * 
 * Parameters:
 *             cellwall object, use coordinate array and glycosidic bonds
 * 
 */
double compute_energy_gbond(CellWallMonolayer *cwl){
  
  int i, mi, mj;
  double d, x, y , z;
  double g_energy = 0.0f;

  for(i=0; i<cwl->get_total_glycobonds()*2; i+=2){
    mi = cwl->glyco_bonds[i];
    mj = cwl->glyco_bonds[i+1];

    x = cwl->coordinate_xyz[mj] - cwl->coordinate_xyz[mi];
    y = cwl->coordinate_xyz[mj+1] - cwl->coordinate_xyz[mi+1];
    z = cwl->coordinate_xyz[mj+2] - cwl->coordinate_xyz[mi+2];

    d = sqrt(x*x + y*y + z*z) - d0_g;

    g_energy += d*d;
  }

  g_energy = g_energy * 0.5f * stiffness_g;

  return g_energy;

}

/*
 * Compute the total energy of peptidic springs according to
 * the equation of spring : e = 0.5 * k_p * (||r2-r1||-d_p)^2
 * (euclidian norme used)
 * 
 * Same function as compute_energy_gbond but could be changed 
 * according to the evolution of the model so better it be 
 * separed
 * 
 * Parameters:
 *             cellwall object, use coordinate array and peptidic bonds
 * 
 */
double compute_energy_pbond(CellWallMonolayer *cwl){
  
  int i, mi, mj;
  double d, x, y , z;
  double p_energy = 0.0f;

  for(i=0; i<cwl->get_total_peptibonds()*2; i+=2){
    mi = cwl->pepti_bonds[i];
    mj = cwl->pepti_bonds[i+1];

    x = cwl->coordinate_xyz[mj] - cwl->coordinate_xyz[mi];
    y = cwl->coordinate_xyz[mj+1] - cwl->coordinate_xyz[mi+1];
    z = cwl->coordinate_xyz[mj+2] - cwl->coordinate_xyz[mi+2];

    d = sqrt(x*x + y*y + z*z) - d0_p;

    p_energy += d*d;
  }

  p_energy = p_energy * 0.5f * stiffness_p;

  return p_energy;

}