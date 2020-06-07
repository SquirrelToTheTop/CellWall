
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
  double dij, x, y , z, tmp;
  double g_energy = 0.0f;

  for(i=0; i<cwl->get_total_glycobonds()*2; i+=2){
    mi = cwl->glyco_bonds[i];
    mj = cwl->glyco_bonds[i+1];

    x = cwl->coordinate_xyz[mj] - cwl->coordinate_xyz[mi];
    y = cwl->coordinate_xyz[mj+1] - cwl->coordinate_xyz[mi+1];
    z = cwl->coordinate_xyz[mj+2] - cwl->coordinate_xyz[mi+2];

    dij = sqrt(x*x + y*y + z*z);

    tmp = stiffness_g * (dij - d0_g) / dij;
    dij = dij - d0_g;

    cwl->forces_xyz[mi] += tmp * x;
    cwl->forces_xyz[mi+1] += tmp * y;
    cwl->forces_xyz[mi+2] += tmp * z;

    cwl->forces_xyz[mj] -= tmp * x;
    cwl->forces_xyz[mj+1] -= tmp * y;
    cwl->forces_xyz[mj+2] -= tmp * z;

    g_energy += dij * dij;
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
  double dij, x, y , z, tmp;
  double p_energy = 0.0f;

  for(i=0; i<cwl->get_total_peptibonds()*2; i+=2){
    mi = cwl->pepti_bonds[i];
    mj = cwl->pepti_bonds[i+1];

    x = cwl->coordinate_xyz[mj] - cwl->coordinate_xyz[mi];
    y = cwl->coordinate_xyz[mj+1] - cwl->coordinate_xyz[mi+1];
    z = cwl->coordinate_xyz[mj+2] - cwl->coordinate_xyz[mi+2];

    dij = sqrt(x*x + y*y + z*z);

    tmp = stiffness_g * (dij - d0_p) / dij;
    dij = dij - d0_p;

    cwl->forces_xyz[mi] += tmp * x;
    cwl->forces_xyz[mi+1] += tmp * y;
    cwl->forces_xyz[mi+2] += tmp * z;

    cwl->forces_xyz[mj] -= tmp * x;
    cwl->forces_xyz[mj+1] -= tmp * y;
    cwl->forces_xyz[mj+2] -= tmp * z;

    p_energy += dij * dij;
  }

  p_energy = p_energy * 0.5f * stiffness_p;

  return p_energy;

}


/*
* Compute the total energy of lipidic springs according to
* the equation of spring : e = 0.5 * k_l * (||r2-r1|| - d0_l)^2
* (euclidian norme used)
* 
* Parameters:
*             cellwall object, use coordinate array and lipidic bonds
* 
*/
double compute_energy_lbond(CellWallLipidLayer *ll){
  
  int i, mi, mj;
  double dij, x, y , z, tmp;
  double l_energy = 0.0f;
  double at_rest = ll->get_spring_d0();

  for(i=0; i<ll->get_total_lbonds()*2; i+=2){
    mi = ll->lipidic_bonds[i];
    mj = ll->lipidic_bonds[i+1];

    x = ll->coordinate_xyz[mj] - ll->coordinate_xyz[mi];
    y = ll->coordinate_xyz[mj+1] - ll->coordinate_xyz[mi+1];
    z = ll->coordinate_xyz[mj+2] - ll->coordinate_xyz[mi+2];

    dij = sqrt(x*x + y*y + z*z);

    tmp = stiffness_l * (dij - at_rest) / dij;
    dij = dij - at_rest;

    ll->forces_xyz[mi] += tmp * x;
    ll->forces_xyz[mi+1] += tmp * y;
    ll->forces_xyz[mi+2] += tmp * z;

    ll->forces_xyz[mj] -= tmp * x;
    ll->forces_xyz[mj+1] -= tmp * y;
    ll->forces_xyz[mj+2] -= tmp * z;

    l_energy += dij * dij;
  }

  l_energy = l_energy * 0.5f * stiffness_l;

  return l_energy;

}

/*
 * Compute the total energy of lipid - lipid angles according to
 * the equation  : e = 0.5 * k_ll * (Theta - Theta_0)^2
 * (euclidian norme used). Theorical/optimal angle is PI
 * 
 * Parameters:
 *             lipid layer object, use coordinate array and lipid
 *             ll_angles
 * 
 */
double compute_energy_ll_angles(CellWallLipidLayer *ll){
  
  int i, mi, mj, mk;
  double x_ij, y_ij , z_ij, norm_ij;
  double x_jk, y_jk , z_jk, norm_jk;
  double ll_energy = 0.0f, pscal, theta;

  for(i=0; i<ll->get_total_lipid_lipid_angles()*3; i+=3){

    // the two involved bonds
    mi = ll->ll_angles[i];
    mj = ll->ll_angles[i+1]; // middle mass #emoticon_lunette
    mk = ll->ll_angles[i+2];

    // mi, mj, mk are the index in the coordinate array of the "x" coordinate of the mass
    x_ij = ll->coordinate_xyz[mi] - ll->coordinate_xyz[mj];
    y_ij = ll->coordinate_xyz[mi+1] - ll->coordinate_xyz[mj+1];
    z_ij = ll->coordinate_xyz[mi+2] - ll->coordinate_xyz[mj+2];

    x_jk = ll->coordinate_xyz[mk] - ll->coordinate_xyz[mj];
    y_jk = ll->coordinate_xyz[mk+1] - ll->coordinate_xyz[mj+1];
    z_jk = ll->coordinate_xyz[mk+2] - ll->coordinate_xyz[mj+2];

    norm_ij = sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);
    norm_jk = sqrt(x_jk*x_jk + y_jk*y_jk + z_jk*z_jk);

    // cos_theta = (mimj . mjmk)/||mimj||*||mjmk||
    pscal = x_ij*x_jk + y_ij*y_jk + z_ij*z_jk;
    theta = acos( pscal/(norm_ij*norm_jk) );

    ll_energy += (theta-alpha0_gg)*(theta-alpha0_gg);
    
  }

  ll_energy = ll_energy * 0.5f * stiffness_ll;

  return ll_energy;

}