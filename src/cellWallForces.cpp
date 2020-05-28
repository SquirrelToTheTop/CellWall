
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
