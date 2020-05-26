#ifndef CELLWALLMONOLAYER_CLASS
#define CELLWALLMONOLAYER_CLASS

#include "cellWallParams.h"
#include "math.h"

/*
 * Class to define the entire cell wall structure
 * 
 */
class CellWallMonolayer{

  public:
    CellWallMonolayer(int, int);  // constructor
    ~CellWallMonolayer(); // killor

    void generate_geometry(); // generate position of masses in 3D
    void generate_bonds(); // generate bonds between masses

    // getter
    double * get_coordinate_array(); // get pointor to coordinate array
    int get_total_npg();

  public:
    double * coordinate_xyz = nullptr; // array of coordinate (x1,y1,z1, x2,y2,z2, ...)
    int * glyco_bonds = nullptr;
    int * pepti_bonds = nullptr;

  private:
    int _nstrands = 0; // number of strands
    int _npgstrand = 0; // number of peptidoglycans per strand
    int _total_npg = 0; // total number of peptidoglycans

    double _cellwall_radius = 0.0f;

};

#endif