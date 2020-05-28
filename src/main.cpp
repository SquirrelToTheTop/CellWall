#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"
#include "cellWallForces.h"
#include "cellWallIO.h"

#include "cellWallDebug.h"

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(10,20);

	cwl->generate_geometry();

	cwl->generate_glycosidic_bonds();

	cwl->generate_peptidic_bonds();

#ifdef DEBUG
	display_glyco_bonds(cwl);
	display_pepti_bonds(cwl);
#endif

	double energy_glyco;	
	clock_t start, end;
	double cpu_time_used;

	start = clock();
	energy_glyco = compute_energy_gbond(cwl);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	fprintf(stdout,"\n\t> Total energy of glycosidic bonds : %f nJ - elapsed : %f ", energy_glyco, cpu_time_used);
	fflush(stdout);

	// iosystem->write_coordinate_ascii_PLY(cwl);
	iosystem->write_PDB(cwl);

	delete cwl;
	delete iosystem;

	return EXIT_SUCCESS;
}
