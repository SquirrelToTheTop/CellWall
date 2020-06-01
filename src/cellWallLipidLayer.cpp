#include <stdio.h>

#include "cellWallLipidLayer.h"

/*
 * Constructor
 * 
 * Parameters:
 *              cw_radius : radius of the cellwall
 *              cw_length : total length of the cellwall
 *              cw_nstrands : number of strands in the cell wall
 * 
 * In order to build the lipid layer, we need the cellwall radius, to create a lipid 
 * layer a few nanometer above the cell wall
 * 
 */
CellWallLipidLayer::CellWallLipidLayer(double cw_radius, double cw_length, int cw_nstrands){

  int i;

  fprintf(stdout, "\n> Ask for new lipid layer ! \n");

  if(  cw_radius < 0.0f ){
    fprintf(stderr, "\n> Cell wall Radius is negative !\n");
    fflush(stderr);
    return;
  }

  // compute layer radius according to cell wall radius
  // 90% of the cell wall radius could/should be adapted
  _layer_radius = 0.90f * cw_radius;
  _layer_length = cw_length;

   // could be 2 but 1.5 makes less lipids masses
  _nstrands = cw_nstrands * 1.5f;

  d0_l = _layer_length / (_nstrands-1);

  // number of lipid per strands
  _nlpstrand = int(_layer_radius * (2.0f*PI) / d0_l);

  // allocation of masses coordinate array
  coordinate_xyz = new double[DIM * _total_nlp];
  _memory_consumption += (sizeof(coordinate_xyz[0]) * (DIM *_total_nlp)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of coordinate_xyz array");
  fflush(stdout);

  // allocation of masses forces array
  forces_xyz = new double[DIM * _total_nlp];
  _memory_consumption += (sizeof(forces_xyz[0]) * (DIM * _total_nlp)) / (MBYTES);
  
  // set all forces to 0.0f at the begining of the simulation
  for(i=0; i<DIM * _total_nlp; ++i )
    forces_xyz[i] = 0.0f;

  fprintf(stdout, "\n\t> Allocation of forces_xyz array");
  fflush(stdout);

  // x2 because one bond is made of two masses
  // and lipidic_bonds will contains the masses "x" index in coordinate_xyz
  lipidic_bonds = new int[_nlipidic_bonds * 2];
  _memory_consumption += (sizeof(lipidic_bonds[0]) * (_nlipidic_bonds * 2)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of glycosidic bonds array");
  fflush(stdout);

  fprintf(stdout, "\n\t> Memory consumption : %f MBytes \n", _memory_consumption);
  fflush(stdout);

}

/*
 * Destructor
 *
 */
CellWallLipidLayer::~CellWallLipidLayer(){

  if( lipidic_bonds ){
    delete [] lipidic_bonds;
    fprintf(stdout, "\n\t> Deallocation of lipidic_bonds array");
    fflush(stdout);
  }

  if( forces_xyz ){
    delete [] forces_xyz;
    fprintf(stdout, "\n\t> Deallocation of forces_xyz array");
    fflush(stdout);
  }

  if( coordinate_xyz ){
    delete [] coordinate_xyz;
    fprintf(stdout, "\n\t> Deallocation of coordinate_xyz array");
    fflush(stdout);
  }

}

/*
 * Output information about the simulation parameters
 * 
 */
void CellWallLipidLayer::simulation_infos(){

  fprintf(stdout, "\n\t> Simulation parameters for lipid layer : \n");
  fprintf(stdout, "\n\t\t> Number of strands : %d", _nstrands);
  fprintf(stdout, "\n\t\t> Number of lp/strands : %d\n", _nlpstrand);
  fprintf(stdout, "\n\t\t> Lipid layer radius: %f nm", _layer_radius);
  fprintf(stdout, "\n\t\t> Lipid layer length: %f nm\n", _layer_length);
  fprintf(stdout, "\n\t\t> Total number of LP: %d", _total_nlp);
  fprintf(stdout, "\n\t\t> Total number of LP-springs: %d\n", _nlipidic_bonds);
  fflush(stdout);

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
 *     /    void generate_glyco_glyco_angles(); // generate array of glycosidic - glycosidic angles

 *    x
 *
 * x: varies in clock like
 * y: strand coordinate wich increase from a lenght of peptidique spring 
 * z: varies with projection 
 */   
void CellWallLipidLayer::generate_geometry(){

  fprintf(stdout, "\n\t> Generate lipid masses position in 3D coordinate");
  fflush(stdout);

  int i, j, offset, offset_0;
  double ystrand;
  double alpha = 0.0f;
  double dalpha = (2.0f*PI) / _nlpstrand;

  // make first strand, then for all strand y coordinate and z are the same
  offset=0;
  for(i=0; i<_nlpstrand; ++i){
    coordinate_xyz[offset] = _layer_radius * cos(alpha); // x coordinate
    coordinate_xyz[offset+1] = 0.0f; // y coordinate (first strand @ (x,0,z))
    coordinate_xyz[offset+DIM-1] = _layer_radius * sin(alpha); // z coordinate
    alpha += dalpha;

    offset += DIM;
  }

  // set all y coordinate and propagate coordinate of first strand to others
  offset= _nlpstrand*DIM; // x coordinate of first mass of second strand
  ystrand = d0_p;
  for(i=1; i<_nstrands; ++i){
    
    offset_0 = 0; // offset first strand to loop over it
    for(j=0; j<_nlpstrand; ++j){
      coordinate_xyz[offset] = coordinate_xyz[offset_0]; // x coordinate
      coordinate_xyz[offset+1] = ystrand; // y coordinate
      coordinate_xyz[offset+2] = coordinate_xyz[offset_0+2]; // z coordinate

      offset += DIM;
      offset_0 += DIM;
    }
    
    // increase y coordinate by a peptidic distance
    ystrand += d0_l;

  }

}

/*
 * Generate array of connection for glycosidic bond
 * 
 *        z
 *        |  
 *           4  9  *  *  *  *
 *           |  |  |  |  |  | 
 *           3  8  *  *  *  * 
 *           |  |  |  |  |  |
 *           2  7  *  *  *  *    -> y
 *           |  |  |  |  |  |
 *           1  6  *  *  *  *
 *           |  |  |  |  |  | 
 *           0  5  *  *  *  * 
 *       /
 *      x
 * 
 * the glyco_bond array is constructed in a way that glyco_bond[i] and glyco_bond[i+1]
 * are the couple of masses linked together. 
 *  
 */   
void CellWallLipidLayer::generate_bonds(){

  int i, j, mi, mj, ai, bi;

  fprintf(stdout, "\n\t> Generate lipidic springs ... ");
  fflush(stdout);

  // glycosidic bonds
  mi = 0;
  mj = DIM;
  bi = 0;
  ai = 0;
  for(i=0; i<_nstrands; ++i){
    
    for(j=0; j<_nlpstrand-1; ++j){
      lipidic_bonds[bi] = mi;
      lipidic_bonds[bi+1] = mj;

      // gg_angles[ai] = int(mi/3);
      // gg_angles[ai+1] = int(mj/3);

      mi += DIM;
      mj += DIM;
      bi += 2;
      ai += 2;
    }

    // link last mass of strand and first of strand (circle)
    lipidic_bonds[bi] = mj - DIM;
    lipidic_bonds[bi+1] = mi - (_nlpstrand-1)*DIM;
    bi += 2;

    // gg_angles[ai] = int((mj-DIM)/3);
    // gg_angles[ai+1] = int((mi - (_nlpstrand-1)*DIM)/3);
    // ai += 2;

    // next strands
    mi += DIM;
    mj += DIM;
  }

  if( int(bi/2) != _nlipidic_bonds ){
    fprintf(stderr, "\n\t> Error: number of lipidic springs created != theorical number of bonds !");
    fprintf(stderr, "\t> %d vs %d\n", int(bi)/2, _nlipidic_bonds );
    fflush(stdout);
  }else{
    if( int(ai)/2 != _nlipid_lipid_angles ){
      fprintf(stderr, "\n\n\t> Error: number of l-l angles created != theorical number of g-g angles !\n");
      fprintf(stderr, "\t> %d vs %d\n", int(ai)/2, _nlipid_lipid_angles );
      fflush(stdout);
    }else{
      fprintf(stdout, " SUCCESS \n");
      fflush(stdout);
    }
  }  

}
