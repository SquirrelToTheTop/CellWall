#include <stdio.h>

#include "cellWallObject.h"

/*
 * Constructor
 * 
 * Parameters:
 *              nstrands : number of strands for the cell wall model
 *              npgstrand : number of masses per strand
 */
CellWallMonolayer::CellWallMonolayer(int nstrands, int npgstrand){

  fprintf(stdout, "\n> Ask for new cell wall monolayer! \n");

  if( nstrands == 0 || npgstrand == 0 ){
    fprintf(stderr, "\n> Number of strands or number of PG/strands must be != 0\n");
    fflush(stderr);
    return;
  }

  if( npgstrand%2 != 0 ){
    npgstrand++;
    fprintf(stderr, "\n> Number of PG/strands must be pair because of the bond model\n");
    fflush(stderr);
  }

  _nstrands = nstrands;
  _npgstrand = npgstrand;
  _total_npg = nstrands * npgstrand;

  _cellwall_radius = d0_g * npgstrand / (2.0f*PI);

  // total number of glycosidic bond 
  _nglyco_bonds = npgstrand * nstrands;

  _npepti_bonds = (nstrands-1) * int(npgstrand/2);

  // allocation of masses coordinate array
  coordinate_xyz = new double[DIM * _total_npg];
  _memory_consumption += (sizeof(coordinate_xyz[0])*(DIM*_total_npg)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of coordinate_xyz array");
  fflush(stdout);

  glyco_bonds = new int[_nglyco_bonds];
  _memory_consumption += (sizeof(glyco_bonds[0])*(_nglyco_bonds)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of glycosidic bonds array");
  fflush(stdout);

  pepti_bonds = new int[_npepti_bonds];
  _memory_consumption += (sizeof(pepti_bonds[0])*(_npepti_bonds)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of peptidic bonds array");
  fflush(stdout);

  fprintf(stdout, "\n\t> Memory consumption : %f MBytes \n", _memory_consumption);
  fflush(stdout);

}

/*
 * Destructor
 *
 */
CellWallMonolayer::~CellWallMonolayer(){

  if( coordinate_xyz ){
    delete [] coordinate_xyz;
    fprintf(stdout, "\n\t> Deallocation of coordinate_xyz array");
    fflush(stdout);
  }

  if( glyco_bonds ){
    delete [] glyco_bonds;
    fprintf(stdout, "\n\t> Deallocation of glyco_bonds array");
    fflush(stdout);
  }

  if( pepti_bonds ){
    delete [] pepti_bonds;
    fprintf(stdout, "\n\t> Deallocation of pepti_bonds array");
    fflush(stdout);
  }

}

/*
 * Generate all mass position in 3D coordinate system
 * 
 *        z
 *        | 
 *        |  *  *  *  *  *  *
 *        |  *  *  *  *  *  * 
 *        |__*__*__*__*__*__* y
 *       /   *  *  *  *  *  *
 *      /    *  *  *  *  *  *
 *     /
 *    x
 *
 * x: varies in clock like
 * y: strand coordinate wich increase from a lenght of peptidique spring 
 * z: varies with projection 
 */   
void CellWallMonolayer::generate_geometry(){

  fprintf(stdout, "\n\t> Generate cell wall masses position in 3D coordinate");
  fflush(stdout);

  int i, j, offset, offset_0;
  double ystrand;
  double alpha = 0.0f;
  double dalpha = (2.0f*PI) / _npgstrand;

  // make first strand, then for all strand y coordinate and z are the same
  offset=0;
  for(i=0; i<_npgstrand; ++i){
    coordinate_xyz[offset] = _cellwall_radius * cos(alpha); // x coordinate
    coordinate_xyz[offset+1] = 0.0f; // y coordinate (first strand @ (x,0,z))
    coordinate_xyz[offset+DIM-1] = _cellwall_radius * sin(alpha); // z coordinate
    alpha += dalpha;

    offset += DIM;
  }

  // set all y coordinate and propagate coordinate of first strand to others
  offset= _npgstrand*DIM; // x coordinate of first mass of second strand
  ystrand = d0_p;
  for(i=1; i<_nstrands; ++i){
    
    offset_0 = 0; // offset first strand to loop over it
    for(j=0; j<_npgstrand; ++j){
      coordinate_xyz[offset] = coordinate_xyz[offset_0]; // x coordinate
      coordinate_xyz[offset+1] = ystrand; // y coordinate
      coordinate_xyz[offset+2] = coordinate_xyz[offset_0+2]; // z coordinate

      offset += DIM;
      offset_0 += DIM;
    }
    
    // increase y coordinate by a peptidic distance
    ystrand += d0_p;

  }

  // just test the distance between two point of first strand
  // to compare to the distance at rest of a glycosidic spring in 
  // order to see if its already hardly bend because a small number of 
  // mass. Stupid for now since the cell wall radius is computed considering
  // input but later cell wall radius will be defined by input of bacteria
  // then we must be sure of the geometrical state of the spring are not to
  // far from rest
  
  double actual_d0_g = 2.0f * (_cellwall_radius * tan(dalpha*0.5f));
  double inf, sup;

  inf = d0_g * (1.0f-epsilon_g);
  sup = d0_g * (1.0f+epsilon_g);
  if( actual_d0_g > sup || actual_d0_g < inf  ){
    fprintf(stdout, "\n\t> WARNING: geometrical size of glyco-spring to far from rest ");
    fprintf(stdout, "\n\t> WARNING: borne : %f < %f  < %f", inf, actual_d0_g, sup);
    fprintf(stdout, "\n\t> WARNING: rest %f <-----> current %f ", d0_g, actual_d0_g);
    fflush(stdout);
  }

  fprintf(stdout, "\n\t\t> actual d0_g : %f", actual_d0_g);
  fprintf(stdout, "\n\t\t> theorical d0_g : %f\n", d0_g);
  fflush(stdout);

}

/*
 * Generate all mass position in 3D coordinate system
 * 
 *        z
 *        |  
 *           4  9--*  *--*  *
 *           |  |  |  |  |  | 
 *           3--8  *--*  *--* 
 *           |  |  |  |  |  |
 *           2  7--*  *--*  *    -> y
 *           |  |  |  |  |  |
 *           1--6  *--*  *--*
 *           |  |  |  |  |  | 
 *           0  5--*  *--*  * 
 *       /
 *      x
 *  
 * x: varies in clock like
 * y: strand coordinate wich increase from a lenght of peptidique spring 
 * z: varies with projection 
 */   
void CellWallMonolayer::generate_bonds(){

  fprintf(stdout, "\n\t> Generate cell wall masses bonds \n");
  fflush(stdout);




}



/* Getter */

double * CellWallMonolayer::get_coordinate_array(){
  return coordinate_xyz;
}

int CellWallMonolayer::get_total_npg(){
  return _total_npg;
}

