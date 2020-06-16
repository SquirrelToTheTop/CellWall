#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "cellWallUtils.h"
#include "cellWallConfig.h"
#include "cellWallObject.h"
#include "cellWallLipidLayer.h"
#include "cellWallForces.h"
#include "cellWallIO.h"
#include "cellWallOptimization.h"

#ifdef DEBUG
#include "cellWallDebug.h"
#endif

int main(int argc, char *argv[]){

	int prank, nrank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &prank);

	printf("Hello World from process %d of %d\n", prank, nrank);

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);


	if( prank == 0 ){
		welcome_message();

		CellWallMonolayer *cwl = new CellWallMonolayer(10,50);
		CellWallLipidLayer *llayer = new CellWallLipidLayer(cwl->get_radius(), cwl->get_length(), 
																												cwl->get_number_of_strands());

		cwl->simulation_infos();
		cwl->generate_geometry();
		cwl->generate_glycosidic_bonds();
		cwl->generate_peptidic_bonds();

		llayer->simulation_infos();
		llayer->generate_geometry();
		llayer->generate_bonds();
		llayer->generate_mesh();

#ifdef DEBUG
		display_glyco_bonds(cwl);
		display_pepti_bonds(cwl);
		display_glyco_glyco_angles(cwl);
		display_lipid_lipid_angles(llayer);
		display_lipid_mesh(llayer);
#endif

		// analyze_cpu_time_energy(cwl, llayer, 300);

		// optimize_simulated_annealing(cwl, llayer, 10000);
		optimize_simulated_annealing_force(cwl, llayer, 1);

		delete llayer;
		delete cwl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
