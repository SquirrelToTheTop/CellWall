#include <stdio.h>
#include <stdlib.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"
#include "cellWallIO.h"

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(50, 100);

	cwl->generate_geometry();

	cwl->generate_bonds();

	// iosystem->write_coordinate_ascii_PLY(cwl);
	iosystem->write_PDB(cwl);

	delete cwl;
	delete iosystem;

	return EXIT_SUCCESS;
}
