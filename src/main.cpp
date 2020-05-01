#include <stdio.h>
#include <stdlib.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(5, 10);

	return EXIT_SUCCESS;
}
