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

	CellWallMonolayer *cwl = new CellWallMonolayer(20,50);

	cwl->generate_geometry();

	cwl->generate_glycosidic_bonds();

	cwl->generate_peptidic_bonds();

#ifdef DEBUG
	display_glyco_bonds(cwl);
	display_pepti_bonds(cwl);
	display_glyco_glyco_angles(cwl);
#endif

	double energy_glyco, energy_pepti, energy_glyco_glyco;
	clock_t start, end;
	double cpu_time_used;

	start = clock();
	energy_glyco = compute_energy_gbond(cwl);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	fprintf(stdout,"\n\t> Total energy of glycosidic bonds : %f nJ - elapsed : %f ", energy_glyco, cpu_time_used);
	fflush(stdout);

	start = clock();
	energy_glyco_glyco = compute_energy_gg_angles(cwl);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	fprintf(stdout,"\n\t> Total energy of glyco-glyco angles : %f nJ - elapsed : %f ", energy_glyco_glyco, cpu_time_used);
	fflush(stdout);

	// Should be null at begining because by construction the distance between two
	// strand is equal to d0_p which is the lenght of the spring at rest
	start = clock();
	energy_pepti = compute_energy_pbond(cwl);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	fprintf(stdout,"\n\t> Total energy of peptidic bonds : %f nJ - elapsed : %f ", energy_pepti, cpu_time_used);
	fflush(stdout);

	// iosystem->write_coordinate_ascii_PLY(cwl);
	iosystem->write_PDB(cwl);

	delete cwl;
	delete iosystem;

	return EXIT_SUCCESS;
}
