#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"
#include "cellWallLipidLayer.h"
#include "cellWallForces.h"
#include "cellWallIO.h"

#ifdef DEBUG
#include "cellWallDebug.h"
#endif

int main(int argc, char *argv[]){

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(4,4);
	CellWallLipidLayer *llayer = new CellWallLipidLayer(cwl->get_radius(), cwl->get_length(), 
	                                                    cwl->get_number_of_strands());

	cwl->simulation_infos();
	cwl->generate_geometry();
	cwl->generate_glycosidic_bonds();
	cwl->generate_peptidic_bonds();

	llayer->simulation_infos();
	llayer->generate_geometry();
	llayer->generate_bonds();

#ifdef DEBUG
	display_glyco_bonds(cwl);
	display_pepti_bonds(cwl);
	display_glyco_glyco_angles(cwl);
	display_lipid_lipid_angles(llayer);
#endif

	double energy_glyco, energy_pepti, energy_glyco_glyco, energy_lipid;
	clock_t start, end;
	double cpu_time_used;

	int iter, nitermax = 2;

	for(iter=0; iter<nitermax; ++iter){

		fprintf(stdout,"\n> Iteration # %d ", iter);
		fflush(stdout);

		start = clock();
		energy_glyco = compute_energy_gbond(cwl);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

		fprintf(stdout,"\n\t> Total energy of glycosidic springs : %f nJ - elapsed : %f ", energy_glyco, cpu_time_used);
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

		fprintf(stdout,"\n\t> Total energy of peptidic springs : %f nJ - elapsed : %f ", energy_pepti, cpu_time_used);
		fflush(stdout);

		start = clock();
		energy_lipid = compute_energy_lbond(llayer);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

		fprintf(stdout,"\n\t> Total energy of lipidic springs : %f nJ - elapsed : %f \n", energy_lipid, cpu_time_used);
		fflush(stdout);


		// mini integration 
		// double dt=0.00001f;
		// for(i=0; i<cwl->get_total_npg()*DIM; i+=DIM){
		// 	cwl->coordinate_xyz[i] -= dt * cwl->forces_xyz[i];
		// 	cwl->coordinate_xyz[i+1] -= dt * cwl->forces_xyz[i+1];
		// 	cwl->coordinate_xyz[i+2] -= dt * cwl->forces_xyz[i+2];
		// }

		// iosystem->write_coordinate_ascii_PLY(cwl);
		iosystem->write_PDB(cwl, iter);
		iosystem->write_PDB(llayer, iter);

	}

	delete llayer;
	delete cwl;
	delete iosystem;

	return EXIT_SUCCESS;
}
