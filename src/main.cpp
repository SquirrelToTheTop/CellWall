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
	int nstrands, npgstrand;

	CellWallMonolayer *cwl;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &prank);

	nstrands = 5 * nrank;
	npgstrand = 16;

	if( prank == 0 ){
		welcome_message();
		fprintf(stdout, "\n> Number of MPI process involved : %d", nrank);
		fprintf(stdout,"\n\t> Number of strands : %d", nstrands);
		fprintf(stdout,"\n\t> Number of PG strands : %d\n", npgstrand);
	}

	if( nrank < 2 ){
		cwl = new CellWallMonolayer(nstrands, npgstrand);
	}else{

		int nstrand_per_proc = nstrands/nrank;
		if( nstrands%nrank != 0 ){
			fprintf(stdout, "\n\t> Load balance != per process");
			fflush(stdout);
		}

		cwl = new CellWallMonolayer(nstrand_per_proc, npgstrand, prank, nrank);
	}
	
	// CellWallLipidLayer *llayer = new CellWallLipidLayer(cwl->get_radius(), cwl->get_length(), 
	// 																										cwl->get_number_of_strands());

	cwl->simulation_infos();
	cwl->generate_geometry();
	cwl->generate_glycosidic_bonds();
	cwl->generate_peptidic_bonds();

	CellWallIOSystem *cio = new CellWallIOSystem("out");
	cio->write_PDB(cwl, prank);
	delete cio;

	// llayer->simulation_infos();
	// llayer->generate_geometry();
	// llayer->generate_bonds();
	// llayer->generate_mesh();

#ifdef DEBUG
		display_glyco_bonds(cwl);
		display_pepti_bonds(cwl);
		display_glyco_glyco_angles(cwl);
		display_lipid_lipid_angles(llayer);
		display_lipid_mesh(llayer);
#endif

	// analyze_cpu_time_energy(cwl, llayer, 300);

	// optimize_simulated_annealing(cwl, llayer, 10000);
	// optimize_simulated_annealing_force(cwl, llayer, 1);

	// delete llayer;
	delete cwl;
	

	MPI_Finalize();

	return EXIT_SUCCESS;
}
