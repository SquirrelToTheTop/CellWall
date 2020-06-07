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
  _layer_radius = 0.80f * cw_radius;
  _layer_length = cw_length;

   // could be 2 but 1.5 makes less lipids masses
  _nstrands = cw_nstrands * 1.5f;

  d0_l = _layer_length / (_nstrands-1);

  // number of lipid per strands
  _nlpstrand = int(_layer_radius * (2.0f*PI) / d0_l);

  // total number of lipid masses
  _total_nlp = _nlpstrand * _nstrands;

  // total number of springs (bonds)
  _nlipidic_bonds = _total_nlp + (_nstrands-1)*_nlpstrand;

  // total number of lipid-lipid angles (only 'Pi' angles)
  _nlipid_lipid_angles = _total_nlp + (_nstrands-2)*_nlpstrand;

  // total number of mesh element
  _nlipidic_mesh = (_total_nlp - _nlpstrand)*2;

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
  fprintf(stdout, "\n\t> Allocation of lipidic springs array");
  fflush(stdout);

  // array of index of masses for lipidic-lipidic angles. Warning, different from 
  // glyco_glyco_angles  and pepti_glyco_angles !
  // so an angle is made with 3 masses !
  ll_angles = new int[_nlipid_lipid_angles * 3];
  _memory_consumption += (sizeof(ll_angles[0]) * (_nlipid_lipid_angles * 3)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of lipid-lipid (Pi) angles array");
  fflush(stdout);

  // array of index of masses for lipidi mesh, 
  // lipidic_mesh[i], lipidic_mesh[i+1], lipidic_mesh[i+2] is one triangle
  lipidic_mesh = new int[_nlipidic_mesh * 3];
  _memory_consumption += (sizeof(lipidic_mesh[0]) * (_nlipidic_mesh * 3)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of lipidic mesh array");
  fflush(stdout);

  fprintf(stdout, "\n\t> Memory consumption : %f MBytes \n", _memory_consumption);
  fflush(stdout);

}

/*
 * Destructor
 *
 */
CellWallLipidLayer::~CellWallLipidLayer(){

  if( lipidic_mesh ){
    delete [] lipidic_mesh;
    fprintf(stdout, "\n\t> Deallocation of lipidic_mesh array");
    fflush(stdout);    
  }

  if( ll_angles ){
    delete [] ll_angles;
    fprintf(stdout, "\n\t> Deallocation of ll_angles array");
    fflush(stdout);
  }

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
  fprintf(stdout, "\n\t\t> Total number of LP-springs: %d", _nlipidic_bonds);
  fprintf(stdout, "\n\t\t> Total number of 'Pi' angles: %d", _nlipid_lipid_angles);
  fprintf(stdout, "\n\t\t> Total number of mesh element: %d\n", _nlipidic_mesh);
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
  ystrand = d0_l;
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
 * Generate array of connection for lipidic spring
 * 
 *        z
 *        |  
 *           4--9--*--*--*--*
 *           |  |  |  |  |  | 
 *           3--8--*--*--*--* 
 *           |  |  |  |  |  |
 *           2--7--*--*--*--*    -> y
 *           |  |  |  |  |  |
 *           1--6--*--*--*--*
 *           |  |  |  |  |  | 
 *           0--5--*--*--*--* 
 *       /
 *      x
 * 
 * the lipidic_bond array is constructed in a way that lipidic_bond[i] and lipidic_bond[i+1]
 * are the couple of masses linked together. 
 *  
 */   
void CellWallLipidLayer::generate_bonds(){

  int i, j, mi, mj, ai, bi;
  int offset;

  fprintf(stdout, "\n\t> Generate lipid springs & angles ... ");
  fflush(stdout);

  // lipidic bonds
  mi = 0;
  mj = DIM;
  bi = 0;

  // offset to same masse on next/previous strand
  offset = _nlpstrand*DIM;

  for(i=0; i<_nstrands-1; ++i){
    
    for(j=0; j<_nlpstrand-1; ++j){
      
      // on the strand
      lipidic_bonds[bi] = mi;
      lipidic_bonds[bi+1] = mj;
      bi += 2;

      // to the right strand
      lipidic_bonds[bi] = mi;
      lipidic_bonds[bi+1] = mi+offset;

      mi += DIM;
      mj += DIM;
      bi += 2;

    }

    // link last mass of strand and first of strand (circle)
    lipidic_bonds[bi] = mj - DIM;
    lipidic_bonds[bi+1] = mi - (_nlpstrand-1)*DIM;
    bi += 2;

    // plus to the right
    lipidic_bonds[bi] = (mj - DIM);
    lipidic_bonds[bi+1] = mi+offset;
    bi += 2;

    // next strands
    mi += DIM;
    mj += DIM;
  }

  // last strand
  for(j=0; j<_nlpstrand-1; ++j){
    
    // on the strand
    lipidic_bonds[bi] = mi;
    lipidic_bonds[bi+1] = mj;
    bi += 2;

    mi += DIM;
    mj += DIM;

  }

  // link last mass of strand and first of strand (circle)
  lipidic_bonds[bi] = mj - DIM;
  lipidic_bonds[bi+1] = mi - (_nlpstrand-1)*DIM;
  bi += 2;

  // Pi angle between tree aligned lipid masses
  ai = 0;
  mi = DIM;

  // first strand 
  // first angle on strand 
  ll_angles[ai] = offset - DIM;
  ll_angles[ai+1] = mi - DIM; // angle here
  ll_angles[ai+2] = mi;
  ai += 3;

  for(i=1; i<_nlpstrand-1; ++i){
    ll_angles[ai] = mi - DIM;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi + DIM;
    ai += 3;
    mi += DIM;
  }

  // last angle on strand 
  ll_angles[ai] = mi - DIM;
  ll_angles[ai+1] = mi; // angle here
  ll_angles[ai+2] = mi - offset + DIM;
  ai += 3;
  mi += DIM;

  // middle strand + angle with masses on left and right strand
  for(i=1; i<_nstrands-1; ++i){

    // first angle on strand 
    ll_angles[ai] = (i+1) * _nlpstrand * DIM - DIM;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi + DIM;
    ai += 3;

    ll_angles[ai] = mi - offset;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi + offset;
    ai += 3;

    mi += DIM;

    for(j=1; j<_nlpstrand-1; ++j){
      
      // angle on same strand
      ll_angles[ai] = mi - DIM;
      ll_angles[ai+1] = mi; // angle here
      ll_angles[ai+2] = mi + DIM;
      ai += 3;

      // angle with mass of left and on right stands
      ll_angles[ai] = mi - offset;
      ll_angles[ai+1] = mi; // angle here
      ll_angles[ai+2] = mi + offset;
      ai += 3;

      mi += DIM;
    }

    // skip last one for the moment
    // last angle on strand 
    ll_angles[ai] = mi - DIM;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi - offset + DIM;
    ai += 3;

    // angle with mass of left and on right stands
    ll_angles[ai] = mi - offset;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi + offset;
    ai += 3;
    mi += DIM;
  
  }

  // last strand
  // first angle on last strand 
  ll_angles[ai] = _nstrands * _nlpstrand * DIM - DIM;
  ll_angles[ai+1] = mi; // angle here
  ll_angles[ai+2] = mi + DIM;
  ai += 3;
  mi += DIM;

  for(i=1; i<_nlpstrand-1; ++i){
    ll_angles[ai] = mi - DIM;
    ll_angles[ai+1] = mi; // angle here
    ll_angles[ai+2] = mi + DIM;
    ai += 3;
    mi += DIM;
  }

  // last angle on last strand 
  ll_angles[ai] = mi - DIM;
  ll_angles[ai+1] = mi; // angle here
  ll_angles[ai+2] = mi - offset + DIM;
  ai += 3;
  mi += DIM;

  if( int(bi/2) != _nlipidic_bonds ){
    fprintf(stderr, "\n\t> Error: number of lipidic springs created != theorical number of bonds !");
    fprintf(stderr, "\t> %d vs %d\n", int(bi)/2, _nlipidic_bonds );
    fflush(stdout);
  }else{
    if( int(ai)/3 != _nlipid_lipid_angles ){
      fprintf(stderr, "\n\n\t> Error: number of l-l angles created != theorical number of l-l angles !\n");
      fprintf(stderr, "\t> %d vs %d\n", int(ai)/3, _nlipid_lipid_angles );
      fflush(stdout);
    }else{
      fprintf(stdout, " SUCCESS ");
      fflush(stdout);
    }
  }

}

/*
 * Generate array of mesh element for the lipid layer
 * 
 *        z
 *        |  
 *           4---9---*---*---*---*
 *           | \ | \ | \ | \ | \ | 
 *           3---8---*---*---*---* 
 *           | \ | \ | \ | \ | \ |
 *           2---7---*---*---*---*    -> y
 *           | \ | \ | \ | \ | \ |
 *           1---6---*---*---*---*
 *           | \ | \ | \ | \ | \ | 
 *           0---5---*---*---*---* 
 *       /
 *      x
 * 
 * the lipidic_bond array is constructed in a way that lipidic_mesh[i], lipidic_mesh[i+1]
 * and lipidic_mesh[i+2] are the couple of masses that form a triangle of the mesh. 
 *  
 */   
void CellWallLipidLayer::generate_mesh(){

  int i, j, mi, bi;
  int offset;

  fprintf(stdout, "\n\t> Generate mesh of lipid layer ... ");
  fflush(stdout);

  // lipidic bonds
  mi = 0;
  bi = 0;

  // offset to same masse on next/previous strand
  offset = _nlpstrand*DIM;

  for(i=0; i<_nstrands-1; ++i){
    
    for(j=0; j<_nlpstrand-1; ++j){
      
      /* first triangle
      * mi+DIM
      *   |  \
      *   |   \
      *  mi -- mi+offset
      */
      lipidic_mesh[bi] = mi;
      lipidic_mesh[bi+1] = mi + offset;
      lipidic_mesh[bi+2] = mi + DIM;
      bi += 3;

      /* second triangle 
      * mi+DIM -- mi + DIM + offset
      *     \     |
      *      \    |
      * mi     -- mi + offset
      */
      lipidic_mesh[bi] = mi + DIM;
      lipidic_mesh[bi+1] = mi + DIM + offset;
      lipidic_mesh[bi+2] = mi + offset;

      mi += DIM;
      bi += 3;

    }

    /* at the end of the strand connect last mass with first
    * first triangle
    * mi-offset+DIM
    *   |  \
    *   |   \
    *  mi -- mi+offset
    */
    lipidic_mesh[bi] = mi;
    lipidic_mesh[bi+1] = mi + offset;
    lipidic_mesh[bi+2] = mi - offset + DIM;
    bi += 3;

    /* second triangle 
    * mi-offset+DIM -- mi + DIM
    *     \                |
    *      \               |
    * mi            -- mi + offset
    */
    lipidic_mesh[bi] = mi - offset + DIM;
    lipidic_mesh[bi+1] = mi + DIM ;
    lipidic_mesh[bi+2] = mi + offset;

    mi += DIM;
    bi += 3;

  }

  if( int(bi)/3 != _nlipidic_mesh ){
    fprintf(stderr, "\n\n\t> Error: number of mesh element created != theorical number of mesh element !\n");
    fprintf(stderr, "\t> %d vs %d\n", int(bi)/3, _nlipidic_mesh );
    fflush(stdout);
  }else{
    fprintf(stdout, " SUCCESS \n");
    fflush(stdout);
  }

}


/* getter */

double * CellWallLipidLayer::get_coordinate_array(){
  return coordinate_xyz;
}

int * CellWallLipidLayer::get_lipidic_bonds_array(){
  return lipidic_bonds;
}

int CellWallLipidLayer::get_total_lipids(){
  return _total_nlp;
}

int CellWallLipidLayer::get_total_lbonds(){
  return _nlipidic_bonds;
}

int CellWallLipidLayer::get_total_lipid_lipid_angles(){
  return _nlipid_lipid_angles;
}

int CellWallLipidLayer::get_total_number_of_mesh(){
  return _nlipidic_mesh;
}

double CellWallLipidLayer::get_spring_d0(){
  return d0_l;
}