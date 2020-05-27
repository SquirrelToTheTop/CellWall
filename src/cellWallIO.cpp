#include "cellWallIO.h"

/*
 * Constructor
 */
CellWallIOSystem::CellWallIOSystem(std::string filename){

  _filename = filename;

  fprintf(stdout, "\n\t> Initialize IO system !");
  fflush(stdout);

}

/*
 * Destructor
 *
 */
CellWallIOSystem::~CellWallIOSystem(){

  
}

/*
 * Open file stream
 * 
 * Parameters:
 *             filename: string containing filename with extension
 *
 */
void CellWallIOSystem::open_file(std::string filename_ext){

  // check if not already open
  if( ! _output_file.is_open() ){

    _output_file.open(filename_ext, std::ios::out);
    if( ! _output_file.is_open() ){
      fprintf(stderr, "\n\t> Error cannot open file : %s ", filename_ext.c_str());
      fflush(stderr);
    }

    _output_file << std::setprecision(std::numeric_limits<double>::digits10);
  
  }else{
    fprintf(stderr, "\n\t> Error in open_file : %s, ofstream already open !", filename_ext.c_str());
    fflush(stderr);
  }

}

/*
 * Close file stream
 * 
 * Parameters:
 *             filename_ext: string containing filename with extension
 *
 */
void CellWallIOSystem::close_file(std::string filename_ext){
  
  if( _output_file.is_open() ){
    _output_file.close();

    fprintf(stdout, "\n\t> Closing output file : %s \n", filename_ext.c_str());
    fflush(stdout);
  }

}

/*
 * Write ASCII file with X Y Z coordinate of masses
 * 
 * file format : xyz
 * 
 * Parameters:
 *             *cwm : pointer to the CellWallObject
 *
 */
void CellWallIOSystem::write_coordinate_ascii(CellWallMonolayer *cwm){

  int i;

  double *xyz = cwm->get_coordinate_array();
  
  open_file(_filename+xyz_ext);
  
  for(i=0; i< cwm->get_total_npg()*DIM; i+=DIM){
    _output_file << std::setw(20) << xyz[i] << "\t" << xyz[i+1] << "\t" << xyz[i+2] << std::endl;
  }

  close_file(_filename + xyz_ext);

  fprintf(stdout, "\n\t> Write output file : %s (#%d)", _filename.c_str(), _nwrite);
  fflush(stdout);

}

/*
 * Write ASCII file with X Y Z coordinate of masses
 * 
 * file format: .ply
 * 
 * Parameters:
 *             *cwm : pointer to the CellWallObject
 *
 */
void CellWallIOSystem::write_coordinate_ascii_PLY(CellWallMonolayer *cwm){

  int i;

  double *xyz = cwm->get_coordinate_array();

  open_file(_filename + ply_ext);

  _output_file << "ply" << std::endl;
  _output_file << "format ascii 1.0" << std::endl;
  _output_file << "element vertex " << cwm->get_total_npg() << std::endl;
  _output_file << "property float x" << std::endl;
  _output_file << "property float y" << std::endl;
  _output_file << "property float z" << std::endl;
  _output_file << "end_header" << std::endl;

  for(i=0; i< cwm->get_total_npg()*DIM; i+=DIM){
    _output_file << std::setw(20) << xyz[i] << "\t" << xyz[i+1] << "\t" << xyz[i+2] << std::endl;
  }

  close_file(_filename + ply_ext);

  fprintf(stdout, "\n\t> Write output file : %s (#%d)", _filename.c_str(), _nwrite);
  fflush(stdout);

}

/*
 * Write ASCII file with X Y Z coordinate of masses and connection
 * between masses
 * 
 * file format: .pdb
 * 
 * Parameters:
 *             *cwm : pointer to the CellWallObject
 *
 */
void CellWallIOSystem::write_PDB(CellWallMonolayer *cwm){

  int i;

  double *xyz = cwm->get_coordinate_array();

  std::string tex = "ATOM";
  std::string res = "GLY";
  std::string typ = "LAM";
  std::string segid = "CP";
  double w1=1.0f, w2=0.0f;
  int ipg =1;
  int idm=1;

  open_file(_filename + pdb_ext);

  _output_file << std::fixed;

  char buff[13];
  std::string buffAsStdStr;

  for(i=0; i< cwm->get_total_npg()*DIM; i+=DIM){
 
    _output_file << std::left << std::setw(6) << tex;
    
    sprintf(&buff[0], "%5d", idm);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;
    
    _output_file << "  " << std::setw(3) << typ << " " << std::setw(3) << res << " ";

    sprintf(&buff[0], "%5d", ipg);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr << "   ";

    // because C++ sucks so damn much 
    sprintf(&buff[0], "%8.3f", xyz[i]);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;

    sprintf(&buff[0], "%8.3f", xyz[i+1]);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;

    sprintf(&buff[0], "%8.3f", xyz[i+2]);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;

    sprintf(&buff[0], "%6.2f", w1);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;

    sprintf(&buff[0], "%6.2f", w2);
    buffAsStdStr = buff;
    _output_file << buffAsStdStr;

    _output_file << " " << std::fixed << std::setw(4) << segid << std::endl;

    idm++;
  }

  close_file(_filename + pdb_ext);

  fprintf(stdout, "\n\t> Write output file : %s (#%d)", (_filename+pdb_ext).c_str(), _nwrite);
  fflush(stdout);

}