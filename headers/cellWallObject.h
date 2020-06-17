#ifndef CELLWALLMONOLAYER_CLASS
#define CELLWALLMONOLAYER_CLASS

#include <mpi.h>

#include "cellWallParams.h"
#include "math.h"

#include "cellWallUtils.h"

/*
 * Class to define the entire cell wall structure
 * 
 */
class CellWallMonolayer{

  public:
    CellWallMonolayer(int, int);  // constructor
    CellWallMonolayer(int, int, int, int); // constructor for MPI
    ~CellWallMonolayer(); // killor

    void clean_forces(); // clean force array

    void generate_geometry(); // generate position of masses in 3D
    void generate_glycosidic_bonds(); // generate glycosidic bonds between masses
    void generate_peptidic_bonds(); // generate peptidic bonds between masses

    // getter
    int get_total_npg(); // get total of masses
    int get_number_of_ghost_pg(); // get total of ghost masses
    int get_total_glycobonds(); // get total number of glycosidic bonds
    int get_total_gg_angles(); // get total number of glycosidic - glycosidic angles
    int get_total_peptibonds(); // get total number of peptidic bonds
    double get_cw_radius(); // get computed radius
    int get_number_of_strands(); // get number of strands
    int get_number_of_pg_strand(); // get number of pg per strands
    double get_radius(); // get the radius of the cell wall
    double get_length(); // get the length of the cell wall
    
    double * get_coordinate_array(); // get pointor to coordinate array
    int * get_glycosidic_bonds_array(); // get pointor to glycosidic array
    int * get_peptidic_bonds_array(); // get pointor to peptidic array

    // display import information about the cellwall model
    void simulation_infos();

  private:
    void cellwall_parameters();
  
  public:

    // array of coordinate (x1,y1,z1, x2,y2,z2, ...)
    double * coordinate_xyz = nullptr;

    // array of forces acting on each masses (fx1,fy1,fz1, fx2,fy2,fz2, ...)
    double * forces_xyz = nullptr;

    // arrays for glycosidic & peptidics bonds index of arrays gives the bonds number
    // and the arrays are sorted like (i1,j1, i2,j2, i3,j3, in,jn) where i1,j1 is a couple
    // of masses index for the first bond, i2,j2 the second bond, etc. 
    // i and j index refer to the "x" coordinate index in the coordinate_xyz thus
    // if i1 is the mass @ (x1,y1,z1) and j2 the mass @ (x2,y2,z2) then
    // i1=0, j1=3
    int * glyco_bonds = nullptr; 
    int * pepti_bonds = nullptr;

    // Array of bonds index which are connected to the same mass, thus forming
    // an angle. So: gg_angles[i] and gg_angles[i+1] refer to the i-th bond in glyco_bonds
    // and the (i+1)-th bond in glyco_bonds array
    int * gg_angles = nullptr;

    // memory info in MBytes
    double _memory_consumption = 0.0f;

  private:
    int _nstrands = 0; // number of strands
    int _nghost_strands = 0; // number of ghost strands

    int _npgstrand = 0; // number of peptidoglycans per strand

    int _total_npg = 0; // total number of peptidoglycans
    int _nghost_npg = 0; // number of ghost peptidogycans

    int _nglyco_bonds = 0; // number of glycosidic springs (bond between masses on same strand)
    int _nghost_nglyco_bonds = 0; // number of ghost glycosidic springs

    int _nglyco_glyco_angles = 0; // number of glycosidic-glycosidic angles
    int _nghost_gg_angles = 0; // number of ghost glycosidic-glycosidic angles

    int _npepti_bonds = 0; // number of peptidic bonds (bond between masses on different strand)
    int _nghost_npepti_bonds = 0; // number of ghost peptidic bonds 

    int _mpi_rank = 0; // rang MPI
    int _mpi_size = 1; // total number of process

    double _cellwall_radius = 0.0f;
    double _cellwall_length = 0.0f;

    double _y0 = 0.0f; // y coordinate for first strand

    int _cellwall_cap_nstrands = 0; // number of strand in one cap
    double _cellwall_cap_radius = 0.0f; // radius of caps
};

#endif