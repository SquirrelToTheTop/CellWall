#include <stdio.h>

#include "cellWallMonolayer.h"

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

  // compute parameters
  cellwall_parameters();

  // allocation of masses coordinate array
  coordinate_xyz = new double[DIM * _total_npg];
  _memory_consumption += (sizeof(coordinate_xyz[0])*(DIM*_total_npg)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of coordinate_xyz array");
  fflush(stdout);

  // x2 because one bond is made of two masses
  // and glyco_bonds will contains the masses "x" index in coordinate_xyz
  glyco_bonds = new int[_nglyco_bonds * 2];
  _memory_consumption += (sizeof(glyco_bonds[0])*(_nglyco_bonds*2)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of glycosidic springs array");
  fflush(stdout);

  // same idea as glyco_bonds
  pepti_bonds = new int[_npepti_bonds * 2];
  _memory_consumption += (sizeof(pepti_bonds[0])*(_npepti_bonds*2)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of peptidic bonds array");
  fflush(stdout);

  // array of index of bonds for glycosidic-glycosidic angles
  gg_angles = new int[_nglyco_glyco_angles * 3];
  _memory_consumption += (sizeof(gg_angles[0])*(_nglyco_glyco_angles*3)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of glyco-glyco angles array");
  fflush(stdout);

  // array of index of bonds for glycosidic-glycosidic angles
  gp_angles = new int[_nglyco_pepti_angles * 3];
  _memory_consumption += (sizeof(gp_angles[0])*(_nglyco_pepti_angles*3)) / (MBYTES);
  fprintf(stdout, "\n\t> Allocation of glyco-pepti angles array");
  fflush(stdout);

  // array of index of masses for mesh, 
  // mesh[i], mesh[i+1], mesh[i+2] is one triangle
  mesh = new int[_nmesh_elements * 3];
  _memory_consumption += (sizeof(mesh[0]) * (_nmesh_elements * 3)) / (MBYTES);

  fprintf(stdout, "\n\t> Memory consumption : %f MBytes \n", _memory_consumption);
  fflush(stdout);

}

/*
 * Constructor
 * 
 * Parameters:
 *              nstrands  : number of strands for the cell wall model
 *              npgstrand : number of masses per strand
 *              prank     : process ID
 *              psize     : total number of process
 */
CellWallMonolayer::CellWallMonolayer(int nstrands, int npgstrand, int prank, int psize){

  _mpi_rank = prank;
  _mpi_size = psize;

  fprintf(stdout, "\n> New CW monolayer for process %d ! \n", _mpi_rank);

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

  // compute parameters
  cellwall_parameters();

  // allocation of masses coordinate array
  coordinate_xyz = new double[DIM * _total_npg];
  _memory_consumption += (sizeof(coordinate_xyz[0])*(DIM*_total_npg)) / (MBYTES);

  // x2 because one bond is made of two masses
  // and glyco_bonds will contains the masses "x" index in coordinate_xyz
  glyco_bonds = new int[_nglyco_bonds * 2];
  _memory_consumption += (sizeof(glyco_bonds[0])*(_nglyco_bonds*2)) / (MBYTES);

  // same idea as glyco_bonds
  pepti_bonds = new int[_npepti_bonds * 2];
  _memory_consumption += (sizeof(pepti_bonds[0])*(_npepti_bonds*2)) / (MBYTES);

  // array of index of masses for glycosidic-glycosidic angles
  gg_angles = new int[_nglyco_glyco_angles * 3];
  _memory_consumption += (sizeof(gg_angles[0])*(_nglyco_glyco_angles*3)) / (MBYTES);

  // array of index of bonds for glycosidic-glycosidic angles
  gp_angles = new int[_nglyco_pepti_angles * 3];
  _memory_consumption += (sizeof(gp_angles[0])*(_nglyco_pepti_angles*3)) / (MBYTES);

  // array of index of masses for mesh, 
  // mesh[i], mesh[i+1], mesh[i+2] is one triangle
  mesh = new int[_nmesh_elements * 3];
  _memory_consumption += (sizeof(mesh[0]) * (_nmesh_elements * 3)) / (MBYTES);

  if( _mpi_size > 1 ){
    fprintf(stdout, "\n\t> Allocation of CW array (P%d) !", _mpi_rank);
    fprintf(stdout, "\n\t> Memory consumption (P%d): %f MBytes \n", _mpi_rank, _memory_consumption);
  }else{
    fprintf(stdout, "\n\t> Allocation of CW array !");
    fprintf(stdout, "\n\t> Memory consumption: %f MBytes \n", _memory_consumption);
  }
  fflush(stdout);

}

/*
 * Destructor for array deallocation
 *
 */
CellWallMonolayer::~CellWallMonolayer(){

  if( mesh )
    delete [] mesh;

  if( gp_angles )
    delete [] gp_angles;

  if( gg_angles )
    delete [] gg_angles;

  if( pepti_bonds )
    delete [] pepti_bonds;

  if( glyco_bonds )
    delete [] glyco_bonds;

  // if( forces_xyz )
  //   delete [] forces_xyz;

  if( coordinate_xyz )
    delete [] coordinate_xyz;

  if( _mpi_size > 1 ){
    fprintf(stdout, "\n\t> Deallocation of CW array (P%d) !\n", _mpi_rank);
  }else{
    fprintf(stdout, "\n\t> Deallocation of CW array !\n");
  }
  fflush(stdout);

}

/*
 * Reset du tableau des forces
 * 
 */
// void CellWallMonolayer::clean_forces(){
//   int i;
  
//   for(i=0; i<_total_npg*DIM; ++i)
//     forces_xyz[i] = 0.0f;

// }

/*
 * Calcul des parametres de la CellWall
 * 
 * _total_npg : nombre total de masses dans le CW
 * _cellwall_radius : rayon initial du CW
 * _cellwall_length : longueur initial du CW
 * 
 * _cellwall_cap_radius : rayon initial pour le CAP (CAP not implemented yet)
 * _cellwall_cap_nstrands : nombre de strands pour le CAP (CAP not implemented yet)
 * 
 * _nglyco_bonds : nombre de ressort de type glycosidic
 * _nglyco_glyco_angles : nombre d'angle type 'PI' entre deux liasisons glyco (meme strand)
 * _npepti_bonds : nombre de ressort de type peptidic
 * _nglyco_pepti_angles : nombre d'angle de type 'PI/2' entre une liaison glyco et une pepti
 * 
 */
void CellWallMonolayer::cellwall_parameters(){

  if( _mpi_size > 1 ){

    if( _mpi_rank > 0 ){
      _nghost_strands = 1;
      _nghost_npg = _npgstrand;
    }

    _nghost_nglyco_bonds = _nghost_strands;
    _nghost_gg_angles = _nghost_nglyco_bonds;
    _nghost_npepti_bonds = _nghost_strands * int(_npgstrand/2);
  }

  _total_npg = _nstrands * _npgstrand + _nghost_npg;

  _cellwall_radius = d0_g * double(_npgstrand) / (2.0f*PI);
  
  _cellwall_cap_radius = _cellwall_radius;
  _cellwall_cap_nstrands = int(_cellwall_cap_radius/d0_p);
  
  // _cellwall_length = (nstrands-1) * d0_p + 2.0f * _cellwall_cap_radius;
  _cellwall_length = double(_nstrands-1) * d0_p;

  // total number of glycosidic springs and ghost bonds not considered for now
  _nglyco_bonds = _npgstrand * _nstrands;

  // there is as much glyco-glyco angles as glycosidic bond and ghost anlges are not considered
  _nglyco_glyco_angles = _nglyco_bonds;

  _npepti_bonds = (_nstrands-1) * int(_npgstrand/2) + _nghost_npepti_bonds;

  // 4 angles for each peptidic spring
  _nglyco_pepti_angles = 4 * _npepti_bonds;

  _nmesh_elements = (_total_npg - _npgstrand) * 2;

}

/*
 * Output information about the simulation parameters
 * 
 */
void CellWallMonolayer::simulation_infos(){

  fprintf(stdout, "\n\t> Simulation parameters for the cellwall : \n");
  fprintf(stdout, "\n\t\t> Number of strands : %d", _nstrands);
  fprintf(stdout, "\n\t\t> Number of pg/strands : %d", _npgstrand);

  if( _mpi_size > 1 ){
    fprintf(stdout, "\n\t\t> Number of ghost strands (P%d): %d ", _mpi_rank, _nghost_strands);
    fprintf(stdout, "\n\t\t> Number of ghost PG      (P%d): %d ", _mpi_rank, _nghost_npg);
  }

  fprintf(stdout, "\n\t\t> Total number of PG: %d\n", _total_npg);
  fprintf(stdout, "\n\t\t> CellWall radius: %f nm", _cellwall_radius);
  fprintf(stdout, "\n\t\t> CellWall length: %f nm\n", _cellwall_length);

  fprintf(stdout, "\n\t\t> Total number of G-springs: %d", _nglyco_bonds);
  fprintf(stdout, "\n\t\t> Total number of G-G angles: %d", _nglyco_glyco_angles);
  fprintf(stdout, "\n\t\t> Total number of P-springs: %d", _npepti_bonds);
  fprintf(stdout, "\n\t\t> Total number of G-P angles: %d\n", _nglyco_pepti_angles);

  fprintf(stdout, "\n\t\t> Total number of mesh elements: %d\n", _nmesh_elements);

  if( _mpi_size > 1 ){
    fprintf(stdout, "\n\t\t> Number of ghost G-springs  (P%d): %d", _mpi_rank, _nghost_nglyco_bonds);
    fprintf(stdout, "\n\t\t> Number of ghost G-G angles (P%d): %d", _mpi_rank, _nghost_gg_angles);
    fprintf(stdout, "\n\t\t> Number of ghost P-springs  (P%d): %d\n", _mpi_rank, _nghost_npepti_bonds);
  }

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
void CellWallMonolayer::generate_geometry(){

  fprintf(stdout, "\n\t> Generate cell wall masses position in 3D coordinate");
  fflush(stdout);

  int i, j, offset, offset_loop, offset_0;
  double ystrand, dy, y_end;
  double alpha = 0.0f;
  double dalpha = (2.0f*PI) / double(_npgstrand);

  dy = d0_p*1.01f; // this makes a spring at rest

  // shift pour les ghosts
  offset = 0;
  offset_0 = 0;
  if( _mpi_size > 1 ){
    offset = _nghost_npg*DIM;
    offset_0 = offset;

    MPI_Barrier(MPI_COMM_WORLD);

    if( _mpi_rank != 0 ){
      MPI_Recv(&_y0, 1, MPI_DOUBLE, _mpi_rank-1, 52, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // _y0 += d0_p;
      _y0 += dy;
    }

    if( _mpi_rank != _mpi_size -1 ){
      y_end = _y0 + (_nstrands-1)*dy;
      MPI_Send(&y_end, 1, MPI_DOUBLE, _mpi_rank+1 , 52, MPI_COMM_WORLD);
    }
    
  }

  // make first strand (after ghost), then for all strand y coordinate and z are the same
  for(i=0; i<_npgstrand; ++i){
    coordinate_xyz[offset] = _cellwall_radius * cos(alpha); // x coordinate
    coordinate_xyz[offset+1] = _y0; // y coordinate (first strand @ (x,0,z))
    coordinate_xyz[offset+DIM-1] = _cellwall_radius * sin(alpha); // z coordinate
    alpha += dalpha;
    offset += DIM;
  }

  // set all y coordinate and propagate coordinate of first strand to others
  // offset += _npgstrand*DIM; // x coordinate of first mass of second strand
  ystrand = _y0 +  dy;
  for(i=1; i<_nstrands; ++i){
    
    offset_loop = offset_0; // offset first strand to loop over it
    for(j=0; j<_npgstrand; ++j){
      coordinate_xyz[offset] = coordinate_xyz[offset_loop]; // x coordinate
      coordinate_xyz[offset+1] = ystrand; // y coordinate
      coordinate_xyz[offset+2] = coordinate_xyz[offset_loop+2]; // z coordinate

      offset += DIM;
      offset_loop += DIM;
    }
    
    // increase y coordinate by a peptidic distance
    ystrand += dy;

  }

  // Interprocess communication
  if( _mpi_size > 1 ){

    MPI_Barrier(MPI_COMM_WORLD);

    if( _mpi_rank != _mpi_size -1 ){
      MPI_Send(&coordinate_xyz[(_total_npg-_npgstrand)*DIM], _npgstrand*DIM, MPI_DOUBLE, 
               _mpi_rank+1, 52, MPI_COMM_WORLD);
    }

    if( _mpi_rank != 0 ){
      MPI_Recv(&coordinate_xyz[0], _npgstrand*DIM, MPI_DOUBLE, _mpi_rank-1, 52, 
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
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
void CellWallMonolayer::generate_glycosidic_bonds(){

  int i, j, mi, mj, ai, bi;
  int offset;

  fprintf(stdout, "\t> Generate cell wall glycosidic springs ... ");
  fflush(stdout);

  // we skip ghost pg and do not take into account glyco spring
  // since its taken into account on _mpi_rank-1 process
  if( _mpi_size > 1 && _mpi_rank > 0 ){
    mi = _npgstrand*DIM;
    mj = mi+DIM;
  }else{
    mi = 0;
    mj = DIM;
  }

  bi = 0;
  for(i=0; i<_nstrands; ++i){
    
    for(j=0; j<_npgstrand-1; ++j){
      glyco_bonds[bi] = mi;
      glyco_bonds[bi+1] = mj;

      mi += DIM;
      mj += DIM;
      bi += 2;

    }

    // link last mass of strand and first of strand (circle)
    glyco_bonds[bi] = mj - DIM;
    glyco_bonds[bi+1] = mi - (_npgstrand-1)*DIM;

    bi += 2;

    // next strands
    mi += DIM;
    mj += DIM;
  }

  // Pi angle between tree aligned pg masses
  ai = 0;
  if( _mpi_size > 1 && _mpi_rank > 0 ){
    mi = _nghost_strands*_npgstrand*DIM;
  }else{
    mi = 0;
  }

  offset = _npgstrand*DIM;

  for(i=0; i<_nstrands; ++i){

    //first angles on first mass of strand
    gg_angles[ai] =  mi + offset - DIM;
    gg_angles[ai+1] = mi; // attached here, angle here
    gg_angles[ai+2] = mi+DIM;
    ai += 3;
    mi += DIM;

    for(j=1; j<_npgstrand-1; ++j){
      gg_angles[ai] =  mi - DIM;
      gg_angles[ai+1] = mi; // attached here, angle here
      gg_angles[ai+2] = mi+DIM;
      ai += 3;
      mi += DIM;
    }

    gg_angles[ai] =  mi - DIM;
    gg_angles[ai+1] = mi; // attached here, angle here
    gg_angles[ai+2] = mi - offset + DIM;
    ai += 3;
    mi += DIM;

  }

  if( int(bi/2) != _nglyco_bonds ){
    fprintf(stderr, "\n\t> Error: number of glycosidc bonds created != theorical number of bonds !");
    fflush(stdout);
  }else{
    if( int(ai)/3 != _nglyco_glyco_angles ){
      fprintf(stderr, "\n\n\t> Error: number of g-g angles created != theorical number of g-g angles !\n");
      fprintf(stderr, "\t> %d vs %d\n", int(ai)/3, _nglyco_glyco_angles );
      fflush(stdout);
    }else{
      fprintf(stdout, " SUCCESS \n");
      fflush(stdout);
    }
  }

}

/*
 * Generate array of connection for peptidic bond
 * 
 *        z
 *        |   
 *           3  7--*  *--*  * 
 *           |  |  |  |  |  |
 *           2--6  *--*  *--*    -> y
 *           |  |  |  |  |  |
 *           1  5--9  *--*  *
 *           |  |  |  |  |  | 
 *           0--4  8--*  *--* 
 *       /
 *      x
 *  
 */   
void CellWallMonolayer::generate_peptidic_bonds(){

  int i, j, mi, mj, ai, bi;
  int offset, laststrand_pg;

  fprintf(stdout, "\t> Generate cell wall peptidic springs ... ");
  fflush(stdout);

  // start @ DIM and +DIM, for the condition if( i%2 == 1 )
  mi=DIM;
  mj=(_npgstrand*DIM)+DIM;

  offset = _npgstrand*DIM;
  laststrand_pg = offset-DIM;

  bi=0;
  ai=0;
  // take into account ghost pg
  for(i=0; i<(_nstrands+_nghost_strands)-1; ++i){

    // this allow to switch between odd and even strand
    // to prevent linking the same mass on the strand before and after
    // see in description: example masse 0 -- 4 and 4 -- 8 with the switch
    // the first link will be 5 -- 9 
    if( i%2 == 0 ){
      mi -= DIM;
      mj -= DIM;
    }else{
      mi += DIM;
      mj += DIM;
    }

    for(j=mi; j<((i+1)*_npgstrand)*DIM; j+=(2*DIM)){

      // mi belong to current strand
      // mj belong to the next strand (left to right, y=0 to y=n)

      pepti_bonds[bi] = mi;
      pepti_bonds[bi+1] = mj;

      if( (mi%(_npgstrand*DIM)) == 0 ){

        gp_angles[ai] = mi+offset-DIM;
        gp_angles[ai+1] = mi; // angle here
        gp_angles[ai+2] = mj;
        ai += 3;

        gp_angles[ai] = mi+DIM;
        gp_angles[ai+1] = mi; // angle here
        gp_angles[ai+2] = mj;
        ai += 3;

        gp_angles[ai] = mj+offset-DIM;
        gp_angles[ai+1] = mj; // angle here
        gp_angles[ai+2] = mi;
        ai += 3;

        gp_angles[ai] = mj+DIM;
        gp_angles[ai+1] = mj; // angle here
        gp_angles[ai+2] = mi;
        ai += 3;

      }else{

        if( mi == laststrand_pg ){

          gp_angles[ai] = mi-DIM;
          gp_angles[ai+1] = mi; // angle here
          gp_angles[ai+2] = mj;
          ai += 3;

          gp_angles[ai] = mi-offset+DIM;
          gp_angles[ai+1] = mi; // angle here
          gp_angles[ai+2] = mj;
          ai += 3;

          gp_angles[ai] = mj-DIM;
          gp_angles[ai+1] = mj; // angle here
          gp_angles[ai+2] = mi;
          ai += 3;

          gp_angles[ai] = mj-offset+DIM;
          gp_angles[ai+1] = mj; // angle here
          gp_angles[ai+2] = mi;
          ai += 3;

        }else{

          gp_angles[ai] = mi-DIM;
          gp_angles[ai+1] = mi; // angle here
          gp_angles[ai+2] = mj;
          ai += 3;

          gp_angles[ai] = mi+DIM;
          gp_angles[ai+1] = mi; // angle here
          gp_angles[ai+2] = mj;
          ai += 3;

          gp_angles[ai] = mj-DIM;
          gp_angles[ai+1] = mj; // angle here
          gp_angles[ai+2] = mi;
          ai += 3;

          gp_angles[ai] = mj+DIM;
          gp_angles[ai+1] = mj; // angle here
          gp_angles[ai+2] = mi;
          ai += 3;

        }

      }

      mi += DIM+DIM;
      mj += DIM+DIM;

      bi += 2;

    }

    laststrand_pg += offset;

  }

  if( int(bi/2) != _npepti_bonds ){
    fprintf(stderr, "\n\n\t> Error: N peptidic bonds created (%d) != theorical N %d !", int(bi/2), _npepti_bonds);
    fflush(stdout);
  }else{
    if( int(ai)/3 != _nglyco_pepti_angles ){
      fprintf(stderr, "\n\n\t> Error: number of g-p angles created != theorical number of g-p angles !\n");
      fprintf(stderr, "\t> %d vs %d\n", int(ai)/3, _nglyco_pepti_angles );
      fflush(stdout);
    }else{
      fprintf(stdout, " SUCCESS \n");
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
 *           0---5---*---*---*---* Iteration # %d ramdonly selected PG : %d", iter, rand_id);
			fprintf(stdout, "\n\t\t> Force x : %f", cwl->forces_xyz[rand_id]);
			fprintf(stdout, "\n\t\t> Force y : %f", cwl->forces_xyz[rand_id+1]);
			fprintf(stdout, "\n\t\t> Force z : %f", cwl->forces_xyz[rand_id+2]);

			// move the selected PG in the force direction 
			cwl->coordinate_xyz[rand_id] +=  1.0f * cwl->forces_xyz[rand_id];
			cwl->coordinate_xyz[rand_id+1] += 1.0f * cwl->forces_xyz[rand_id+1];
			cwl->coordinate_xyz[ran
 *       /
 *      x
 * 
 *  
 */   
void CellWallMonolayer::generate_mesh(){

  int i, j, mi, bi;
  int offset;

  fprintf(stdout, "\n\t> Generate mesh element ... ");
  fflush(stdout);

  mi = 0;
  bi = 0;

  // offset to same masse on next/previous strand
  offset = _npgstrand*DIM;

  for(i=0; i<_nstrands-1; ++i){
    
    for(j=0; j<_npgstrand-1; ++j){
      
      /* first triangle
      * mi+DIM
      *   |  \
      *   |   \
      *  mi -- mi+offset
      */
      mesh[bi] = mi;
      mesh[bi+1] = mi + offset;
      mesh[bi+2] = mi + DIM;
      bi += 3;

      /* second triangle 
      * mi+DIM -- mi + DIM + offset
      *     \     |
      *      \    |
      * mi     -- mi + offset
      */
      mesh[bi] = mi + DIM;
      mesh[bi+1] = mi + DIM + offset;
      mesh[bi+2] = mi + offset;

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
    mesh[bi] = mi;
    mesh[bi+1] = mi + offset;
    mesh[bi+2] = mi - offset + DIM;
    bi += 3;

    /* second triangle 
    * mi-offset+DIM -- mi + DIM
    *     \                |
    *      \               |
    * mi            -- mi + offset
    */
    mesh[bi] = mi - offset + DIM;
    mesh[bi+1] = mi + DIM ;
    mesh[bi+2] = mi + offset;

    mi += DIM;
    bi += 3;

  }

  if( int(bi)/3 != _nmesh_elements ){
    fprintf(stderr, "\n\n\t> Error: number of mesh element created != theorical number of mesh element !\n");
    fprintf(stderr, "\t> %d vs %d\n", int(bi)/3, _nmesh_elements );
    fflush(stdout);
  }else{
    fprintf(stdout, " SUCCESS \n");
    fflush(stdout);
  }

}

/* Getter */

double * CellWallMonolayer::get_coordinate_array(){
  return coordinate_xyz;
}

int * CellWallMonolayer::get_glycosidic_bonds_array(){
  return glyco_bonds;
}

int * CellWallMonolayer::get_peptidic_bonds_array(){
  return pepti_bonds;
}

int CellWallMonolayer::get_total_npg(){
  return _total_npg;
}

int CellWallMonolayer::get_number_of_ghost_pg(){
  return _nghost_npg;
}

int CellWallMonolayer::get_number_of_ghost_strand(){
  return _nghost_strands;
}

int CellWallMonolayer::get_number_of_mesh_elements(){
  return _nmesh_elements;
}

int CellWallMonolayer::get_total_glycobonds(){
  return _nglyco_bonds;
}

int CellWallMonolayer::get_total_gg_angles(){
  return _nglyco_glyco_angles;
}

int CellWallMonolayer::get_total_gp_angles(){
  return _nglyco_pepti_angles;
}

int CellWallMonolayer::get_total_peptibonds(){
  return _npepti_bonds;
}

double CellWallMonolayer::get_cw_radius(){
  return _cellwall_radius;
}

int CellWallMonolayer::get_number_of_strands(){
  return _nstrands;
}

int CellWallMonolayer::get_number_of_pg_strand(){
  return _npgstrand;
}

double CellWallMonolayer::get_radius(){
  return _cellwall_radius;
}

double CellWallMonolayer::get_length(){
  return _cellwall_length;
}
