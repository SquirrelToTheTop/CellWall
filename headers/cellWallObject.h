#ifndef CELLWALLMONOLAYER_CLASS
#define CELLWALLMONOLAYER_CLASS

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
    ~CellWallMonolayer(); // killor

    void generate_geometry(); // generate position of masses in 3D
    void generate_glycosidic_bonds(); // generate glycosidic bonds between masses

    // getter
    int get_total_npg(); // get total of masses
    int get_total_glycobonds(); // get total number of glycosidic bonds
    int get_total_peptibonds(); // get total number of peptidic bonds
    
    double * get_coordinate_array(); // get pointor to coordinate array
    int * get_glycosidic_bonds_array(); // get pointor to glycosidic array

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

    // memory info in MBytes
    double _memory_consumption = 0.0f;

  private:
    int _nstrands = 0; // number of strands
    int _npgstrand = 0; // number of peptidoglycans per strand
    int _total_npg = 0; // total number of peptidoglycans
    int _nglyco_bonds = 0; // number of peptidic bonds (bond between masses on same strand)
    int _npepti_bonds = 0; // number of peptidic bonds (bond between masses on different strand)

    double _cellwall_radius = 0.0f;

};

#endif