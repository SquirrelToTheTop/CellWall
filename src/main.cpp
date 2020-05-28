#include <stdio.h>
#include <stdlib.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"
#include "cellWallIO.h"

#include "cellWallDebug.h"

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(30, 100);

	cwl->generate_geometry();

	cwl->generate_glycosidic_bonds();

	display_glyco_bonds(cwl);

	// iosystem->write_coordinate_ascii_PLY(cwl);
	iosystem->write_PDB(cwl);

	delete cwl;
	delete iosystem;

	return EXIT_SUCCESS;
}
