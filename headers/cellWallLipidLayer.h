#ifndef CELLWALLLIPIDLAYER_CLASS
#define CELLWALLLIPIDLAYER_CLASS

#include "cellWallParams.h"
#include "math.h"

#include "cellWallUtils.h"

/*
 * Class to define the entire lipidic layer structure
 * 
 */
class CellWallLipidLayer{

  public:
    CellWallLipidLayer(double, double, int);  // constructor
    ~CellWallLipidLayer(); // killor

    void generate_geometry(); // generate position of masses in 3D
    void generate_bonds(); // generate lipidic bonds between masses

    double * get_coordinate_array(); // get pointor to coordinate array
    int * get_lipidic_bonds_array(); // get pointor to lipidic array

    // display import information about the cellwall model
    void simulation_infos();
  
  public:

    // array of coordinate (x1,y1,z1, x2,y2,z2, ...)
    double * coordinate_xyz = nullptr;

    // array of forces acting on each masses (fx1,fy1,fz1, fx2,fy2,fz2, ...)
    double * forces_xyz = nullptr;

    // arrays for lipidic bonds index of arrays gives the bonds number
    // and the arrays are sorted like (i1,j1, i2,j2, i3,j3, in,jn) where i1,j1 is a couple
    // of masses index for the first bond, i2,j2 the second bond, etc. 
    // i and j index refer to the "x" coordinate index in the coordinate_xyz thus
    // if i1 is the mass @ (x1,y1,z1) and j2 the mass @ (x2,y2,z2) then
    // i1=0, j1=3
    int * lipidic_bonds = nullptr; 


    // memory info in MBytes
    double _memory_consumption = 0.0f;

  private:
    int _nstrands = 0; // number of strands
    int _nlpstrand = 0; // number of lipid mass per strand
    int _total_nlp = 0; // total number of lipid masses
    int _nlipidic_bonds = 0; // total number of lipid springs
    int _nlipid_lipid_angles = 0; // total number of lipid-lipid angles

    double _layer_radius = 0.0f;
    double _layer_length = 0.0f;
    double _distance_from_cw = 0.0f;

    double d0_l = 0.0f; // need to be computed some how

};

#endif