
#include "cellWallParams.h"

#include "cellWallForces.h"

#include <stdio.h>

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
 * Compute the total energy of glycosidic - glycosidic angles according to
 * the equation  : e = 0.5 * k_gg * (Theta - Theta_0)^2
 * (euclidian norme used)
 * 
 * Parameters:
 *             cellwall object, use coordinate array and glycosidic bonds
 * 
 */
double compute_energy_gg_angles(CellWallMonolayer *cwl){
  
  int i, mi, mj, mk;
  double x_ij, y_ij , z_ij, norm_ij;
  double x_jk, y_jk , z_jk, norm_jk;
  double gg_energy = 0.0f, pscal, theta;

  for(i=0; i<cwl->get_total_gg_angles()*2; i+=2){

    // the two involved bonds
    mi = cwl->glyco_bonds[cwl->gg_angles[i]*2];
    mj = cwl->glyco_bonds[(cwl->gg_angles[i]*2)+1]; // middle mass #emoticon_lunette
    mk = cwl->glyco_bonds[(cwl->gg_angles[i+1]*2)+1];

    // mi, mj, mk are the index in the coordinate array of the "x" coordinate of the mass
    x_ij = cwl->coordinate_xyz[mi] - cwl->coordinate_xyz[mj];
    y_ij = cwl->coordinate_xyz[mi+1] - cwl->coordinate_xyz[mj+1];
    z_ij = cwl->coordinate_xyz[mi+2] - cwl->coordinate_xyz[mj+2];

    x_jk = cwl->coordinate_xyz[mk] - cwl->coordinate_xyz[mj];
    y_jk = cwl->coordinate_xyz[mk+1] - cwl->coordinate_xyz[mj+1];
    z_jk = cwl->coordinate_xyz[mk+2] - cwl->coordinate_xyz[mj+2];

    norm_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
    norm_jk = sqrt(x_jk*x_jk + y_jk*y_jk + z_jk*z_jk);

    // cos_theta = (mimj . mjmk)/||mimj||*||mjmk||
    pscal = x_ij*x_jk + y_ij*y_jk + z_ij*z_jk;
    theta = acos( pscal/(norm_ij*norm_jk) );

    gg_energy += (theta-alpha0_gg)*(theta-alpha0_gg);
  }

  gg_energy = gg_energy * 0.5f * stiffness_gg;

  return gg_energy;

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