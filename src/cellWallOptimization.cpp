#include <stdio.h>
#include <string.h>
#include <time.h>

// c++ hope this isn't shit
#include <random>

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

  double energy_glyco, energy_pepti, energy_lipid, energy_lj;
	double energy_glyco_glyco, energy_lipid_lipid, energy_pressure;
  double total_energy;

	clock_t start, end;
	double cpu_time_used[7];

  std::string * ener_msg = new std::string[7] { "Spring G", "Angles G-G", "Spring P", "Spring L", "Angles L-L", "Pressure", "LJ"};

	int iter;

	for(int i=0; i<7; ++i)
		cpu_time_used[i] = 0.0f;

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

		start = clock();
		energy_lj = compute_energy_lennardJones(cwl, ll);
		end = clock();
		cpu_time_used[6] += ((double) (end - start)) / CLOCKS_PER_SEC;

    total_energy = energy_glyco + energy_glyco_glyco + energy_lipid + energy_lipid_lipid;
    total_energy += energy_pepti + energy_pressure + energy_lj;

	}

  for(iter=0; iter<7; iter++){
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

  double energy_glyco, energy_pepti, energy_lipid, energy_lj;
	double energy_glyco_glyco, energy_lipid_lipid, energy_pressure;
  double total_energy, total_prev_energy;

	int iter, i, niter_reset=0;

	double temperature = fictive_temperature;
	double max_dx, dx, proba_keep;

	// save n-1
	double *cwl_coordinate = new double[cwl->get_total_npg()*DIM];
	double *ll_coordinate = new double[ll->get_total_lipids()*DIM];

	// max distance for mass deplacent at once
	// min_dx = 0.01f * d0_p;
	max_dx = 0.05f * d0_p;

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_real_distribution<> dis(0.0f, 1.0f);//uniform distribution between 0 and 1

	total_prev_energy = 0.0f;
	total_prev_energy += compute_energy_gbond(cwl);
	total_prev_energy += compute_energy_gg_angles(cwl);
	total_prev_energy += compute_energy_pbond(cwl);
	total_prev_energy += compute_energy_lbond(ll);
	total_prev_energy += compute_energy_ll_angles(ll);
	total_prev_energy -= compute_energy_pressure(ll);
	total_prev_energy += compute_energy_lennardJones(cwl, ll);
	 
	for(iter=0; iter<niter; ++iter){

		for(i=0; i<cwl->get_total_npg()*DIM; ++i){
			cwl_coordinate[i] = cwl->coordinate_xyz[i];

			dx = dis(rng)*max_dx;
			// fprintf(stdout,"\n> moving dx : %f", dx);
			if( dis(rng) > 0.05f )
				cwl->coordinate_xyz[i] += dx;
			else
				cwl->coordinate_xyz[i] -= dx;
		}

		for(i=0; i<ll->get_total_lipids()*DIM; ++i){
			ll_coordinate[i] = ll->coordinate_xyz[i];

			dx = dis(rng)*max_dx;
			// fprintf(stdout,"\n> moving dx : %f", dx);
			if( dis(rng) > 0.05f )
				ll->coordinate_xyz[i] += dx;
			else
				ll->coordinate_xyz[i] -= dx;
		}

		energy_glyco = compute_energy_gbond(cwl);
		energy_glyco_glyco = compute_energy_gg_angles(cwl);
		energy_pepti = compute_energy_pbond(cwl);
		energy_lipid = compute_energy_lbond(ll);
		energy_lipid_lipid = compute_energy_ll_angles(ll);
		energy_pressure = compute_energy_pressure(ll);
		energy_lj = compute_energy_lennardJones(cwl, ll);

    // total_energy = energy_glyco + energy_glyco_glyco + energy_pepti;
    total_energy =   energy_lipid + energy_lipid_lipid - energy_pressure; //+ energy_lj;

		fprintf(stdout, "\n\t> G %f P %f L %f GG %f LL %f V %f LJ %f\n", energy_glyco, energy_pepti, energy_lipid, 
		energy_glyco_glyco, energy_lipid_lipid, energy_pressure, energy_lj);

		if( total_energy > total_prev_energy ){
			// compute probability of keep it anyway
			proba_keep = exp(-(total_energy - total_prev_energy) / (fictive_kb*temperature) );

			fprintf(stdout,"\n\t\t Energy k+1 = %f", total_energy);
			fprintf(stdout,"\n\t\t Energy k = %f", total_prev_energy);
			fprintf(stdout,"\n\t\t Delta E = %f", total_energy - total_prev_energy);
			fprintf(stdout,"\n\t\t Proba %f \n", proba_keep);
			fflush(stdout);

			// on reset sinon on accepte quand meme la montée en énergie
			if( dis(rng) > proba_keep ){
				for(i=0; i<cwl->get_total_npg()*DIM; ++i){
					cwl->coordinate_xyz[i] = cwl_coordinate[i];
				}

				for(i=0; i<ll->get_total_lipids()*DIM; ++i){
					ll->coordinate_xyz[i] = ll_coordinate[i];
				}

				niter_reset++;
				if( niter_reset == 200 ){
					niter_reset = 0;
					// reduction de temp
					temperature *= 0.98;
					fprintf(stdout,"\n\t> Reduction of fictive temperarture : %f K\n", temperature);
					fflush(stdout);	
				}
				// fprintf(stdout,"\n> Reset masses positions !\n");

			}else{
				
				total_prev_energy = total_energy;

				// fprintf(stdout,"\n> Energy increased !\n");
				fprintf(stdout,"\n\tAccepted increased energy %d %f nJ \n", iter, total_energy);
				fflush(stdout);
				iosystem->write_PDB(cwl, iter);
        iosystem->write_PDB(ll, iter);

			}

		}else{
			fprintf(stdout,"\n\tLower state of energy %d %f nJ \n", iter, total_energy);
			fflush(stdout);

			total_prev_energy = total_energy;

			iosystem->write_PDB(cwl, iter);
      iosystem->write_PDB(ll, iter);

		}

	}

  delete iosystem;

}