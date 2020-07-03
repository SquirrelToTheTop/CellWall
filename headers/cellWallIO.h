#ifndef CELLWALL_IO_H
#define CELLWALL_IO_H

#define LENGHT_FILENAME 256

// thanks f****** Cplusshit 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>

#include <string>

#include <mpi.h>

#include "cellWallMonolayer.h"
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

    int add_data_file(std::string);

  private:

    const std::string pdb_ext=".pdb";
    const std::string ply_ext=".ply";
    const std::string xyz_ext=".xyz";

    std::ofstream _output_file;
    std::string _filename;

    int _nwrite = 0;
    int _buffer_size_ko = 4; // en Koctets

    int _mpi_rank = 0; // rang MPI
    int _mpi_size = 1; // total number of process

    bool activate_io = true; // fast activation or desactivation of IO

    void open_file(std::string filename);
    void close_file(std::string filename);

};

#endif
