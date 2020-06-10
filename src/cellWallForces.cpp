
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

/*
 * Compute the total energy due to the turgor pressure according to
 * the equation  : e = P . V
 * 
 * Volume computed based on center of a strand (coordinate computed
 * with least square method which equals to a mean) and tree 
 * tetrahedrons for each mesh element (2 triangles)
 * 
 * Parameters:
 *             lipid layer object, use coordinate array and lipid
 *             ll_angles
 * 
 */
double compute_energy_pressure(CellWallLipidLayer *ll){
  
  int i, j, ci, mi;
  int ns, nlps, meshi, nmesh_strand;
  int a, b, c, d;
  double tmp_oa[DIM], tmp_oc[DIM], tmp_ob[DIM], tmp_od[DIM];

  ns = ll->get_number_of_strands();
  nlps = ll->get_number_of_lp_strand();
  nmesh_strand = nlps;

  double centers[ns*DIM];
  double total_v, v1, v2, v3;

  total_v = 0.0f;

  // compute centers of lipid strands least square method
  ci = 0;
  mi = 0;
  for(i=0; i<ns; ++i){

    centers[ci]= 0.0f;
    centers[ci+1]= 0.0f;
    centers[ci+2]= 0.0f;
    
    for(j=0; j<nlps; ++j){
      centers[ci] += ll->coordinate_xyz[mi];
      centers[ci+1] += ll->coordinate_xyz[mi+1];
      centers[ci+2] += ll->coordinate_xyz[mi+2];
      mi += DIM;
    }

    centers[ci] /= nlps;
    centers[ci+1] /= nlps;
    centers[ci+2] /= nlps;

    // fprintf(stdout, "\n> Center #%d @ (%f,%f,%f)", i, centers[ci], centers[ci+1], centers[ci+2]);
    // fflush(stdout);
    ci += DIM;
  }

  ci = 0;
  meshi = 0;
  for(i=0; i<ns-1; ++i){

    // same center for all those mesh element (square mesh here)
    for(j=0; j<nmesh_strand; ++j){

      // first tetrahedron
      a = ll->lipidic_mesh[meshi];
      b = ll->lipidic_mesh[meshi+1];
      c = ll->lipidic_mesh[meshi+2];
      
      tmp_oc[0] = ll->coordinate_xyz[c] - centers[ci];
      tmp_oc[1] = ll->coordinate_xyz[c+1] - centers[ci+1];
      tmp_oc[2] = ll->coordinate_xyz[c+2] - centers[ci+2];

      tmp_oa[0] = ll->coordinate_xyz[a] - centers[ci];
      tmp_oa[1] = ll->coordinate_xyz[a+1] - centers[ci+1];
      tmp_oa[2] = ll->coordinate_xyz[a+2] - centers[ci+2];

      v1 = (tmp_oc[1]*tmp_oa[2] - tmp_oc[2]*tmp_oa[1]) * (ll->coordinate_xyz[b] - ll->coordinate_xyz[a]); 
      v1 += (tmp_oc[2]*tmp_oa[0] - tmp_oc[0]*tmp_oa[2]) * (ll->coordinate_xyz[b+1] - ll->coordinate_xyz[a+1]);
      v1 += (tmp_oc[0]*tmp_oa[1] - tmp_oc[1]*tmp_oa[0]) * (ll->coordinate_xyz[b+2] - ll->coordinate_xyz[a+2]);
      v1 /= 6.0f;

      meshi += 3;
      // fprintf(stdout, "\n> Volume tetrahedron 1 (%d,%d,%d,O1) = %f ", int(a/3), int(b/3), int(c/3), v1);

      // second tetrahedron
      tmp_ob[0] = ll->coordinate_xyz[b] - centers[ci];
      tmp_ob[1] = ll->coordinate_xyz[b+1] - centers[ci+1];
      tmp_ob[2] = ll->coordinate_xyz[b+2] - centers[ci+2];

      v2 = (tmp_oc[1]*tmp_ob[2] - tmp_oc[2]*tmp_ob[1]) * (centers[ci+DIM] - ll->coordinate_xyz[b]); 
      v2 += (tmp_oc[2]*tmp_ob[0] - tmp_oc[0]*tmp_ob[2]) * (centers[ci+DIM+1] - ll->coordinate_xyz[b+1]);
      v2 += (tmp_oc[0]*tmp_ob[1] - tmp_oc[1]*tmp_ob[0]) * (centers[ci+DIM+2] - ll->coordinate_xyz[b+2]);
      v2 /= 6.0f;

      // two triangle used
      // fprintf(stdout, "\n> Volume tetrahedron 2 (%d,%d,O1,O2) = %f", int(c/3), int(b/3), v2);

      // third tetrahedon
      d = ll->lipidic_mesh[meshi+1];
      
      // work with center of next strand
      ci += DIM;
      tmp_oc[0] = ll->coordinate_xyz[c] - centers[ci];
      tmp_oc[1] = ll->coordinate_xyz[c+1] - centers[ci+1];
      tmp_oc[2] = ll->coordinate_xyz[c+2] - centers[ci+2];

      tmp_od[0] = ll->coordinate_xyz[d] - centers[ci];
      tmp_od[1] = ll->coordinate_xyz[d+1] - centers[ci+1];
      tmp_od[2] = ll->coordinate_xyz[d+2] - centers[ci+2];

      v3 = (tmp_od[1]*tmp_oc[2] - tmp_od[2]*tmp_oc[1]) * (ll->coordinate_xyz[b] - ll->coordinate_xyz[d]); 
      v3 += (tmp_od[2]*tmp_oc[0] - tmp_od[0]*tmp_oc[2]) * (ll->coordinate_xyz[b+1] - ll->coordinate_xyz[d+1]);
      v3 += (tmp_od[0]*tmp_oc[1] - tmp_od[1]*tmp_oc[0]) * (ll->coordinate_xyz[b+2] - ll->coordinate_xyz[d+2]);
      v3 /= 6.0f;

      // pointor to center of current strand
      ci -= DIM;
      meshi += 3;
      // fprintf(stdout, "\n> Volume tetrahedron 3 (%d,%d,%d,O2) = %f \n", int(a/3), int(b/3), int(d/3), v3);

      total_v += v1 + v2 + v3;

    }

    ci += DIM;

  }

  // fprintf(stdout, "\n\t\t> Total volume of lipid layer (triangle mesh) :  %f nm^3", total_v);
  // fprintf(stdout, "\n\t\t> Total volume of lipid layer (classic) :  %f nm^3", 
  //         PI*ll->get_radius()*ll->get_radius()*ll->get_length());

  return total_v * inner_pressure;

}

/*
 * Compute energy of interaction between the cellwall and the lipid layer
 * 
 * Terrific in time consuming
 */
double compute_energy_lennardJones(CellWallMonolayer *cwl, CellWallLipidLayer *ll){

  int mi, mj, pg_i, lp_i;
  double xij, yij, zij, dij, tmp;
  double energy_lj, energy_tmp;

  pg_i = 0; 
  energy_lj = 0.0f;
  for(mi=0; mi<cwl->get_total_npg(); ++mi){

    lp_i = 0;
    energy_tmp = 0.0f;
    for(mj=0; mj<ll->get_total_lipids(); ++mj){
      
      xij = cwl->coordinate_xyz[pg_i] - ll->coordinate_xyz[lp_i];
      yij = cwl->coordinate_xyz[pg_i+1] - ll->coordinate_xyz[lp_i+1];
      zij = cwl->coordinate_xyz[pg_i+2] - ll->coordinate_xyz[lp_i+2];
      dij = sqrt(xij*xij + yij*yij + zij*zij);

      tmp = (cut_off/dij) * (cut_off/dij) * (cut_off/dij); // ^3
      tmp = tmp * tmp; // ^6

      energy_tmp += (tmp * tmp) - 2.0f * tmp; // tmp^12 - 2.0 * tmp^6

      lp_i += DIM;
    }

    energy_lj += energy_tmp * epsilon_lj;

    pg_i += DIM;

  }

  return energy_lj;

}