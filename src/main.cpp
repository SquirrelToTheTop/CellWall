#include <stdio.h>
#include <stdlib.h>

#include "cellWallConfig.h"
#include "cellWallObject.h"

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	CellWallMonolayer *cwl = new CellWallMonolayer();

	return EXIT_SUCCESS;
}
