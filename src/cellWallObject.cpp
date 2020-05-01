#include <stdio.h>

#include "cellWallObject.h"

/*
 * Constructor
 */
CellWallMonolayer::CellWallMonolayer(int nstrands, int npgstrand){

  fprintf(stdout, "\n> Ask for new cell wall monolayer! \n");

  if( nstrands == 0 || npgstrand == 0 ){
    fprintf(stderr, "\n> Number of strands or number of PG/strands must be != 0\n");
    return;
  }

  total_npg = nstrands * npgstrand;
  
}

/*
 * Destructor
 *
 */
CellWallMonolayer::~CellWallMonolayer(){

}


