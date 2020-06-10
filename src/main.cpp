#include <stdio.h>
#include <stdlib.h>

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

	printf("\n> Major Version %d \n", cellWall_VERSION_MAJOR);
	printf("\n> Minor Version %d \n", cellWall_VERSION_MINOR);

	welcome_message();

	CellWallMonolayer *cwl = new CellWallMonolayer(50,50);
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

	analyze_cpu_time_energy(cwl, llayer, 1000);

	optimize_simulated_annealing(cwl, llayer, 10);

	delete llayer;
	delete cwl;

	return EXIT_SUCCESS;
}
