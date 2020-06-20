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

// #ifdef DEBUG
#include "cellWallDebug.h"
// #endif

int main(int argc, char *argv[]){

	int prank, nrank;
	int nstrands, npgstrand;

	CellWallMonolayer *cwl;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	MPI_Comm_rank(MPI_COMM_WORLD, &prank);

	nstrands = 512;
	npgstrand = 256;

	if( nrank < 2 ){
		cwl = new CellWallMonolayer(nstrands, npgstrand);
	}else{

		int nstrand_per_proc = nstrands/nrank;
		if( nstrands%nrank != 0 ){
			fprintf(stdout, "\n\t> Load balance != per process");
			fflush(stdout);
		}

		if( prank == 0 ){
			welcome_message();
			fprintf(stdout, "\n> Number of MPI process involved : %d", nrank);
			fprintf(stdout,"\n\t> Number of strands : %d", nstrands);
			fprintf(stdout,"\n\t> Number of PG strands : %d\n", npgstrand);
		}

		cwl = new CellWallMonolayer(nstrand_per_proc, npgstrand, prank, nrank);
	}
	
	// CellWallLipidLayer *llayer = new CellWallLipidLayer(cwl->get_radius(), cwl->get_length(), 
	// 																										cwl->get_number_of_strands());

	cwl->simulation_infos();
	cwl->generate_geometry();
	cwl->generate_glycosidic_bonds();
	cwl->generate_peptidic_bonds();

	// if( prank == 0){
	// 	display_glyco_bonds(cwl);
	// 	display_glyco_glyco_angles(cwl);
	// }

	// benchmarks_energy_cw(cwl, prank, nrank);

	double energy_glyco = compute_energy_gbond(cwl);
	fprintf(stdout, "\n\t> (P%d) Energy G spring : %f ", prank, energy_glyco);
	fflush(stdout);

	double total_energy_glyco = 0.0f;
	MPI_Reduce(&energy_glyco, &total_energy_glyco, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if( prank == 0 ){
		fprintf(stdout, "\n\t> (Master P) Total energy of G spring : %f \n", total_energy_glyco);
		fflush(stdout);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	double energy_pepti = compute_energy_pbond(cwl);
	fprintf(stdout, "\n\t> (P%d) Energy P spring : %f ", prank, energy_pepti);
	fflush(stdout);

	double total_energy_pepti = 0.0f;
	MPI_Reduce(&energy_pepti, &total_energy_pepti, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if( prank == 0 ){
		fprintf(stdout, "\n\t> (Master P) Total energy of P spring : %f \n", total_energy_pepti);
		fflush(stdout);
	}

	double energy_gg = compute_energy_gg_angles(cwl);
	fprintf(stdout, "\n\t> (P%d) Energy G-G angles : %f ", prank, energy_gg);
	fflush(stdout);

	double total_energy_gg = 0.0f;
	MPI_Reduce(&energy_gg, &total_energy_gg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if( prank == 0 ){
		fprintf(stdout, "\n\t> (Master P) Total energy of G-G angles : %f \n", total_energy_gg);
		fflush(stdout);
	}

	CellWallIOSystem *cio;
	cio = new CellWallIOSystem("initial");
	// cio->write_PDB(cwl, prank);
	delete cio;

	// MC_simulated_annealing(cwl, prank, nrank);

	// for(int i=0; i<100; ++i)
	// 	conjugate_gradient(cwl, prank, nrank);

	// cio = new CellWallIOSystem("post_CG");
	// cio->write_PDB(cwl, prank);
	// delete cio;

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
