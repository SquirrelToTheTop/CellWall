#ifndef CELLWALL_IO_H
#define CELLWALL_IO_H

#define LENGHT_FILENAME 256

// thanks f****** Cplusshit 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

#include <string>

#include "cellWallObject.h"
#include "cellWallLipidLayer.h"

/*
 * Class to define the IO system
 * 
 */
class CellWallIOSystem{

  public:
    CellWallIOSystem(std::string);  // constructor
    ~CellWallIOSystem(); // killor

    void write_coordinate_ascii(CellWallMonolayer *);
    void write_coordinate_ascii_PLY(CellWallMonolayer *);

    void write_PDB(CellWallMonolayer *, int);
    void write_PDB(CellWallLipidLayer *, int);

  private:

    const std::string pdb_ext=".pdb";
    const std::string ply_ext=".ply";
    const std::string xyz_ext=".xyz";

    std::ofstream _output_file;
    std::string _filename;

    int _nwrite = 0;

    void open_file(std::string filename);
    void close_file(std::string filename);

};

#endif
