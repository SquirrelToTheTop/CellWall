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
    void generate_mesh(); // generate mesh of the layer (used for turgor pressure computation)

    double * get_coordinate_array(); // get pointor to coordinate array
    int * get_lipidic_bonds_array(); // get pointor to lipidic array
    int get_total_lipids(); // get total number of lipids
    int get_number_of_strands(); // get number of strands
    int get_number_of_lp_strand(); // get number of lp/strand
    int get_total_lbonds(); // get total number of bonds
    int get_total_lipid_lipid_angles(); // get total number of lipid lipid angles
    int get_total_number_of_mesh(); // get total number of mesh element
    double get_spring_d0(); // return computed d0_l for lipidic spring

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

    // Not the same structure as glyco_angles, pepti_angles. Here, we do not
    // use the bonds index but the masses index (in *DIM mode). Firstly, because 
    // we won't remove any springs from the lipid layer. Secondly, we need to 
    // compute some crazy stuff and put the masses index is easier. So let makes 
    // our life easier =)
    int * ll_angles = nullptr;

    // array of triangle mesh index. Computed once, then used to compute the
    // volume of the lipid layer and the turgor pressure forces acting on lipid
    // and on the cellwall indirectly. One triangle is made of 3 masses
    // so lipidic_mesh[i], lipidic_mesh[i+1], lipidic_mesh[i+2] is one triangle
    int * lipidic_mesh = nullptr;

    // memory info in MBytes
    double _memory_consumption = 0.0f;

  private:
    int _nstrands = 0; // number of strands
    int _nlpstrand = 0; // number of lipid mass per strand
    int _total_nlp = 0; // total number of lipid masses
    int _nlipidic_bonds = 0; // total number of lipid springs
    int _nlipid_lipid_angles = 0; // total number of lipid-lipid angles
    int _nlipidic_mesh = 0; // number of mesh element (triangle)

    double _layer_radius = 0.0f;
    double _layer_length = 0.0f;
    double _distance_from_cw = 0.0f;

    double d0_l = 0.0f; // need to be computed some how

};

#endif