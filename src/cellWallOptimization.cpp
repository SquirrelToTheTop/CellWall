#include <stdio.h>
#include <string.h>
#include <time.h>

#include "cellWallForces.h"
#include "cellWallObject.h"
#include "cellWallLipidLayer.h"
#include "cellWallIO.h"
#include "cellWallOptimization.h"

/*
 *
 * 
 * 
 */
void analyze_cpu_time_energy(CellWallMonolayer *cwl, CellWallLipidLayer *ll, int niter){

  double energy_glyco, energy_pepti, energy_lipid;
	double energy_glyco_glyco, energy_lipid_lipid, energy_pressure;
  double total_energy;

	clock_t start, end;
	double cpu_time_used[6];

  std::string * ener_msg = new std::string[6] { "Spring G", "Angles G-G", "Spring P", "Spring L", "Angles L-L", "Pressure"};

	int iter;

	for(iter=0; iter<niter; ++iter){

    if( iter%100 == 0 ){
		  fprintf(stdout,"\n\t> Iteration # %d ", iter);
	  	fflush(stdout);
    }

		start = clock();
		energy_glyco = compute_energy_gbond(cwl);
		end = clock();
		cpu_time_used[0] += ((double) (end - start)) / CLOCKS_PER_SEC;

		start = clock();
		energy_glyco_glyco = compute_energy_gg_angles(cwl);
		end = clock();
		cpu_time_used[1] += ((double) (end - start)) / CLOCKS_PER_SEC;

		// Should be null at begining because by construction the distance between two
		// strand is equal to d0_p which is the lenght of the spring at rest
		start = clock();
		energy_pepti = compute_energy_pbond(cwl);
		end = clock();
		cpu_time_used[2] += ((double) (end - start)) / CLOCKS_PER_SEC;

		start = clock();
		energy_lipid = compute_energy_lbond(ll);
		end = clock();
		cpu_time_used[3] += ((double) (end - start)) / CLOCKS_PER_SEC;

		start = clock();
		energy_lipid_lipid = compute_energy_ll_angles(ll);
		end = clock();
		cpu_time_used[4] += ((double) (end - start)) / CLOCKS_PER_SEC;

		start = clock();
		energy_pressure = compute_energy_pressure(ll);
		end = clock();
		cpu_time_used[5] += ((double) (end - start)) / CLOCKS_PER_SEC;

    total_energy = energy_glyco + energy_glyco_glyco + energy_lipid + energy_lipid_lipid;
    total_energy += energy_pepti + energy_pressure;

	}

  for(iter=0; iter<6; iter++){
    fprintf(stdout, "\n\t> Average elasped time for energy %s : %f s", ener_msg[iter].c_str(), cpu_time_used[iter]/niter);
  }
  fprintf(stdout,"\n");

}

/*
 *
 * 
 * 
 */
void optimize_simulated_annealing(CellWallMonolayer *cwl, CellWallLipidLayer *ll, int niter){

  CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

  double energy_glyco, energy_pepti, energy_lipid;
	double energy_glyco_glyco, energy_lipid_lipid, energy_pressure;
  double total_energy;

	int iter;

	for(iter=0; iter<niter; ++iter){

		fprintf(stdout,"\n> Iteration # %d ", iter);
		fflush(stdout);

		energy_glyco = compute_energy_gbond(cwl);
		energy_glyco_glyco = compute_energy_gg_angles(cwl);
		energy_pepti = compute_energy_pbond(cwl);
		energy_lipid = compute_energy_lbond(ll);
		energy_lipid_lipid = compute_energy_ll_angles(ll);
		energy_pressure = compute_energy_pressure(ll);

    total_energy = energy_glyco + energy_glyco_glyco + energy_lipid + energy_lipid_lipid;
    total_energy += energy_pepti + energy_pressure;

    fprintf(stdout,"\n\t> Total energy : %f nJ \n", total_energy);
		fflush(stdout);

    if( iter%output_iter == 0 ){
      iosystem->write_PDB(cwl, iter);
      iosystem->write_PDB(ll, iter);
    }

	}

  delete iosystem;

}