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
#include "cellWallUtils.h"

/*
 * Benchmarks of energy computation for CW: glycosidic springs, peptidic springs
 * and glycosidic-glycosidic angles
 * 
 * Parameters:
 * 						 cwl : pointor to CW object
 * 						 mpi_rank : ID of current MPI process
 * 						 mpi_size : number of process
 */
void benchmarks_energy_cw(CellWallMonolayer *cwl, int mpi_rank, int mpi_size){

  double energy_glyco, energy_pepti;
	double energy_glyco_glyco;

	clock_t start, end;
	double cpu_time_used[3];

  std::string * ener_msg = new std::string[7] { "Spring G", "Angles G-G", "Spring P"};

	int iter;

	for(int i=0; i<3; ++i)
		cpu_time_used[i] = 0.0f;

	for(iter=0; iter<1000; ++iter){

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

	}

  for(iter=0; iter<3; iter++){
    fprintf(stdout, "\n\t> (P%d) Average elasped time for energy %s : %f s", mpi_rank, ener_msg[iter].c_str(), cpu_time_used[iter]/1000.0f);
  }
  fprintf(stdout,"\n");

}

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
 * Simulated annealing methode with a fictive temperature and boltzmann constante
 * 
 * Not great.
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

/*
 * Simulated annealing methode with a fictive temperature and boltzmann constante
 * 
 * Not great.
 * 
 */
void optimize_simulated_annealing_force(CellWallMonolayer *cwl, CellWallLipidLayer *ll, int niter){

  CellWallIOSystem *iosystem = new CellWallIOSystem("out_test");

  double energy_glyco, energy_pepti, energy_lipid, energy_lj;
	double energy_glyco_glyco, energy_lipid_lipid, energy_pressure;
  double total_energy, total_prev_energy;

	int iter, i, niter_reset=0, id_pg, id_lp;

	double temperature = fictive_temperature;
	double norm_f, proba_keep, tmp;

	// save n-1
	double *cwl_coordinate = new double[cwl->get_total_npg()*DIM];
	double *ll_coordinate = new double[ll->get_total_lipids()*DIM];
	double *cwl_computed_forces = new double[cwl->get_total_npg()*DIM];
	double *ll_computed_forces = new double[ll->get_total_lipids()*DIM];

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_real_distribution<> dis(0.0f, 1.0f);//uniform distribution between 0.0 and 0.5

	cwl->clean_forces();
	ll->clean_forces();

	total_prev_energy = 0.0f;
	energy_glyco = compute_energy_gbond(cwl);
	energy_glyco_glyco = compute_energy_gg_angles(cwl);
	energy_pepti = compute_energy_pbond(cwl);
	energy_lipid = compute_energy_lbond(ll);
	energy_lipid_lipid = compute_energy_ll_angles(ll);
	energy_pressure = compute_energy_pressure(ll);
	compute_force_pressure(ll);
	energy_lj = compute_energy_lennardJones(cwl, ll);

	total_prev_energy = energy_glyco + energy_glyco_glyco + energy_pepti;

	iter = 0;
	iosystem->write_PDB(cwl, iter);
	iosystem->write_PDB(ll, iter);

	fprintf(stdout,"\n\t\t> Total initial energy of CW: %f ", total_prev_energy);
	fprintf(stdout, "\n\t\t> G %f | P %f | GG %f |\n", energy_glyco, energy_pepti, energy_glyco_glyco);

	fprintf(stdout,"\n\t\t> Total energy of lipid Layer: %f", energy_lipid + energy_lipid_lipid - energy_pressure+energy_lj);
	fprintf(stdout, "\n\t\t> L %f | LL %f | V %f | LJ %f\n", energy_lipid, energy_lipid_lipid, -energy_pressure, energy_lj);

	total_prev_energy += energy_lipid + energy_lipid_lipid - energy_pressure + energy_lj;
	 
	for(iter=1; iter<niter; ++iter){

		id_pg = 0;
		for(i=0; i<cwl->get_total_npg(); ++i){
			cwl_coordinate[id_pg] = cwl->coordinate_xyz[id_pg];
			cwl_coordinate[id_pg+1] = cwl->coordinate_xyz[id_pg+1];
			cwl_coordinate[id_pg+2] = cwl->coordinate_xyz[id_pg+2];

			cwl_computed_forces[id_pg] = cwl->forces_xyz[id_pg];
			cwl_computed_forces[id_pg+1] = cwl->forces_xyz[id_pg+1];
			cwl_computed_forces[id_pg+2] = cwl->forces_xyz[id_pg+2];

			norm_f = cwl->forces_xyz[id_pg]*cwl->forces_xyz[id_pg];
			norm_f += cwl->forces_xyz[id_pg+1]*cwl->forces_xyz[id_pg+1];
			norm_f += cwl->forces_xyz[id_pg+2]*cwl->forces_xyz[id_pg+2];
			norm_f = sqrt(norm_f);

			norm_f = 0.002f; // take only 10%

			// fprintf(stdout, "\n> displacement of about : %f", dis(rng) * norm_f * cwl->forces_xyz[id_pg]);

			cwl->coordinate_xyz[id_pg] += norm_f * cwl->forces_xyz[id_pg];
			cwl->coordinate_xyz[id_pg+1] += norm_f * cwl->forces_xyz[id_pg+1];
			cwl->coordinate_xyz[id_pg+2] += norm_f * cwl->forces_xyz[id_pg+2];

			// fprintf(stdout, "\n> displacement : %f", dis(rng) * norm_f * cwl->forces_xyz[id_pg]);

			id_pg += DIM;

		}

		id_lp = 0;
		for(i=0; i<ll->get_total_lipids(); ++i){
			ll_coordinate[id_lp] = ll->coordinate_xyz[id_lp];
			ll_coordinate[id_lp+1] = ll->coordinate_xyz[id_lp+1];
			ll_coordinate[id_lp+2] = ll->coordinate_xyz[id_lp+2];

			ll_computed_forces[id_lp] = ll->forces_xyz[id_lp];
			ll_computed_forces[id_lp+1] = ll->forces_xyz[id_lp+1];
			ll_computed_forces[id_lp+2] = ll->forces_xyz[id_lp+2];

			norm_f = ll->forces_xyz[id_lp]*ll->forces_xyz[id_lp];
			norm_f += ll->forces_xyz[id_lp+1]*ll->forces_xyz[id_lp+1];
			norm_f += ll->forces_xyz[id_lp+2]*ll->forces_xyz[id_lp+2];
			norm_f = sqrt(norm_f);

			norm_f = 0.002f; // take only 10%f

			ll->coordinate_xyz[id_lp] += norm_f * ll->forces_xyz[id_lp];
			ll->coordinate_xyz[id_lp+1] += norm_f * ll->forces_xyz[id_lp+1];
			ll->coordinate_xyz[id_lp+2] += norm_f * ll->forces_xyz[id_lp+2];

			id_lp += DIM;

		}

		// compute new energy after displacement
		cwl->clean_forces();
		ll->clean_forces();

		energy_glyco = compute_energy_gbond(cwl);
		energy_glyco_glyco = compute_energy_gg_angles(cwl);
		energy_pepti = compute_energy_pbond(cwl);
		energy_lipid = compute_energy_lbond(ll);
		energy_lipid_lipid = compute_energy_ll_angles(ll);
		energy_pressure = compute_energy_pressure(ll);
		compute_force_pressure(ll);
		// energy_lj = compute_energy_lennardJones(cwl, ll);

		total_energy = energy_glyco + energy_glyco_glyco + energy_pepti;

		fprintf(stdout, "\n\t> Step # %d ", iter);
		fprintf(stdout,"\n\t\t> Total energy of CW: %f", total_energy);
		fprintf(stdout, "\n\t\t> G %f | P %f | GG %f |\n", energy_glyco, energy_pepti, energy_glyco_glyco);

		fprintf(stdout,"\n\t\t> Total energy of lipid Layer: %f", energy_lipid + energy_lipid_lipid - energy_pressure+energy_lj);
		fprintf(stdout, "\n\t\t> L %f | LL %f | V %f | LJ %f\n", energy_lipid, energy_lipid_lipid, -energy_pressure, energy_lj);

		total_energy += energy_lipid + energy_lipid_lipid - energy_pressure + energy_lj;

		if( total_energy > total_prev_energy ){

			// fprintf(stdout, "\n\t> State of higher energy ! ");

			// compute probability of keep it anyway
			proba_keep = exp(-(total_energy - total_prev_energy) / (fictive_kb*temperature) );

			// fprintf(stdout,"\n\t\t Energy k+1 = %f", total_energy);
			// fprintf(stdout,"\n\t\t Energy k = %f", total_prev_energy);
			// fprintf(stdout,"\n\t\t Delta E = %f", total_energy - total_prev_energy);
			// fprintf(stdout,"\n\t\t Proba to keep %f", proba_keep);
			tmp = dis(rng);
			// fprintf(stdout,"\n\t\t Random number %f \n", tmp);
			// fflush(stdout);

			// Si le rand est inferieur a la proba de keep
			// alors on accepte sinon on reset les positions et les forces pour un new try
			if( tmp < proba_keep ){
				// sinon on accepte quand meme
				total_prev_energy = total_energy;

				fprintf(stdout,"\n\t\t> ACCEPTED ! \n");
				fflush(stdout);

				iosystem->write_PDB(cwl, iter);
        iosystem->write_PDB(ll, iter);

			}else{

				for(i=0; i<cwl->get_total_npg()*DIM; ++i){
					cwl->coordinate_xyz[i] = cwl_coordinate[i];
					cwl->forces_xyz[i] = cwl_computed_forces[i];
				}

				for(i=0; i<ll->get_total_lipids()*DIM; ++i){
					ll->coordinate_xyz[i] = ll_coordinate[i];
					ll->forces_xyz[i] = ll_computed_forces[i];
				}

				fprintf(stdout,"\n\t\t> REFUSED ! \n");
				fflush(stdout);

				niter_reset++;
				if( niter_reset == 200 ){
					niter_reset = 0;
					// reduction de temp
					temperature *= 0.98;
					fprintf(stdout,"\n\t> Reduction of fictive temperarture : %f K\n", temperature);
					fflush(stdout);	
				}

			}

		}else{
			fprintf(stdout,"\n\tState of lower energy %d \n", iter);
			fprintf(stdout,"\n\t\t Energy k+1 = %f", total_energy);
			fprintf(stdout,"\n\t\t Energy k = %f", total_prev_energy);
			fprintf(stdout,"\n\t\t Delta E = %f", total_energy - total_prev_energy);
			fflush(stdout);

			total_prev_energy = total_energy;

			iosystem->write_PDB(cwl, iter);
      iosystem->write_PDB(ll, iter);

		}

	}

  delete iosystem;

}

/*
 * Attempt to use a conjugate gradient
 * 
 */
void conjugate_gradient(CellWallMonolayer *cwl, int mpi_rank, int mpi_size){

	int npg_dim=cwl->get_total_npg()*DIM;
	
	int i, k;
	double cw_energy;;
	double cw_gbond, cw_pepti;

	double beta_k = 0.001f;
	double alpha_k = 0.001f;
	double d_k[npg_dim];

	if( mpi_rank == 0 ){
		fprintf(stdout, "\n\t> Conjugate gradient algorithm: ");
	}

	cw_gbond = compute_energy_gbond(cwl);
	fprintf(stdout,"\n\t\t> (P%d) CW G-bond energy : %f", mpi_rank, cw_gbond);
	fflush(stdout);
	MPI_Reduce(&cw_gbond, &cw_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if( mpi_rank == 0 ){
		fprintf(stdout,"\n\t\t> CW G-spring energy : %f", cw_energy);
		fflush(stdout);
	}

	cw_pepti = compute_energy_pbond(cwl);
	fprintf(stdout,"\n\t\t> (P%d) CW P-spring energy : %f", mpi_rank, cw_pepti);
	fflush(stdout);
	MPI_Reduce(&cw_pepti, &cw_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if( mpi_rank == 0 ){
		fprintf(stdout,"\n\t\t> CW P-bond energy : %f", cw_energy);
		fflush(stdout);
	}

	for(i=0; i<npg_dim; ++i){
		d_k[i] = -cwl->forces_xyz[i];
	}

	for(k=0; k<10; ++k){

		// fprintf(stdout,"\n\t\t> Iteration # %d", k);
		// fflush(stdout);
		for(i=0; i<npg_dim; ++i){
			d_k[i] = -cwl->forces_xyz[i] + beta_k * d_k[i];
		}

	}

	for(i=0; i<npg_dim; ++i){
		cwl->coordinate_xyz[i] = cwl->coordinate_xyz[i] + alpha_k * d_k[i];
	}

	cw_gbond = compute_energy_gbond(cwl);
	fprintf(stdout,"\n\t\t> (P%d) CW G-bond energy : %f", mpi_rank, cw_gbond);
	fflush(stdout);
	MPI_Reduce(&cw_gbond, &cw_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if( mpi_rank == 0 ){
		fprintf(stdout,"\n\t\t> CW G-spring energy : %f", cw_energy);
		fflush(stdout);
	}

	cw_pepti = compute_energy_pbond(cwl);
	fprintf(stdout,"\n\t\t> (P%d) CW P-spring energy : %f", mpi_rank, cw_pepti);
	fflush(stdout);
	MPI_Reduce(&cw_pepti, &cw_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if( mpi_rank == 0 ){
		fprintf(stdout,"\n\t\t> CW P-bond energy : %f", cw_energy);
		fflush(stdout);
	}

	// this need to be checker and modified for MPI
	// cw_gg_angles = compute_energy_gg_angles(cwl);

	


}

/*
 * This must be corrected for MPI
 * 
 */
void MC_simulated_annealing(CellWallMonolayer *cwl, int mpi_rank, int mpi_size){
	
	int rand_id, i;
	double cw_energy_n, cw_energy_np;
	double cw_gbond_n, cw_gbond_np, cw_pepti_n, cw_pepti_np;

	double proba;

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_real_distribution<> dis(0.0f, 1.0f);//uniform distribution between 0.0 and 0.5
	std::uniform_int_distribution<> rand_pg(0,cwl->get_total_npg()-1);

	double *cwl_coordinate = new double[cwl->get_total_npg()*DIM];
	double *cwl_forces = new double[cwl->get_total_npg()*DIM];

	CellWallIOSystem *iosystem = new CellWallIOSystem("MC_test");

	fprintf(stdout, "\n\t> Monte-Carlo Simulated Annealing minimisation");

	cwl->clean_forces();
	cw_gbond_n = compute_energy_gbond(cwl);
	cw_pepti_n = compute_energy_pbond(cwl);
	cw_energy_n = cw_gbond_n + cw_pepti_n;

	fprintf(stdout,"\n\t\t> (P%d) CW G-bond energy 0 : %f", mpi_rank, cw_gbond_n);
	fprintf(stdout,"\n\t\t> (P%d) CW P-spring energy 0 : %f", mpi_rank, cw_pepti_n);
	fprintf(stdout,"\n\t\t> (P%d) CW total energy 0 : %f", mpi_rank, cw_energy_n);
	fprintf(stdout,"\n\t> Norm vector forces : %f\n", norm_vect(&(cwl->forces_xyz[0]),cwl->get_total_npg()*DIM ));

	for(i=0; i<cwl->get_total_npg()*DIM; ++i){
		cwl_coordinate[i] = cwl->coordinate_xyz[i];
		cwl_forces[i] = cwl->forces_xyz[i];
	}

	for(int iter=0; iter<2000000; ++iter){

		rand_id = rand_pg(rng) * DIM;

		// fprintf(stdout, "\n\t> Iteration # %d ramdonly selected PG : %d", iter, rand_id);
		// fprintf(stdout, "\n\t\t> Force x : %f", cwl->forces_xyz[rand_id]);
		// fprintf(stdout, "\n\t\t> Force y : %f", cwl->forces_xyz[rand_id+1]);
		// fprintf(stdout, "\n\t\t> Force z : %f", cwl->forces_xyz[rand_id+2]);

		// move the selected PG in the force direction 
		cwl->coordinate_xyz[rand_id] +=  0.01f * cwl->forces_xyz[rand_id];
		cwl->coordinate_xyz[rand_id+1] += 0.01f * cwl->forces_xyz[rand_id+1];
		cwl->coordinate_xyz[rand_id+2] += 0.01f * cwl->forces_xyz[rand_id+2];

		cwl->clean_forces();
		cw_gbond_np = compute_energy_gbond(cwl);
		cw_pepti_np = compute_energy_pbond(cwl);
		cw_energy_np = cw_gbond_np + cw_pepti_np;

		// check new total energy
		if( cw_energy_np < cw_energy_n ){
			cw_energy_n = cw_energy_np;

			fprintf(stdout, "\n\t\t> Accepted new state # %f", cw_energy_np);
			fprintf(stdout,"\n\t\t> (P%d) CW G-bond energy: %f", mpi_rank, cw_gbond_np);
			fprintf(stdout,"\n\t\t> (P%d) CW P-spring energy: %f", mpi_rank, cw_pepti_np);
			fprintf(stdout,"\n\t\t> (P%d) CW total energy: %f", mpi_rank, cw_energy_np);
			fprintf(stdout,"\n\t\t> Norm vector forces : %f\n", norm_vect(&(cwl->forces_xyz[0]),cwl->get_total_npg()*DIM ));
			iosystem->write_PDB(cwl, iter);

		}else{

			proba = exp(-(cw_energy_np - cw_gbond_n) / (fictive_kb*fictive_temperature) );

			if( dis(rng) < proba ){
				cw_energy_n = cw_energy_np;

				fprintf(stdout, "\n\t\t> Accepted new state # %f - P(n+1) = %f", cw_energy_np, proba);
				fprintf(stdout,"\n\t\t> (P%d) CW G-bond energy: %f", mpi_rank, cw_gbond_np);
				fprintf(stdout,"\n\t\t> (P%d) CW P-spring energy: %f", mpi_rank, cw_pepti_np);
				fprintf(stdout,"\n\t\t> (P%d) CW total energy: %f", mpi_rank, cw_energy_np);
				fprintf(stdout,"\n\t\t> Norm vector forces : %f\n", norm_vect(&(cwl->forces_xyz[0]),cwl->get_total_npg()*DIM ));
				iosystem->write_PDB(cwl, iter);
			}else{

				// fprintf(stdout, "\n\t\t> Refused new state %f -> P(n+1) = %f \n", cw_energy_np, proba);
				cwl->coordinate_xyz[rand_id] = cwl_coordinate[rand_id];
				cwl->coordinate_xyz[rand_id+1] = cwl_coordinate[rand_id+1];
				cwl->coordinate_xyz[rand_id+2] = cwl_coordinate[rand_id+2];

				for(i=0; i<cwl->get_total_npg()*DIM; ++i){
					cwl->forces_xyz[i] = cwl_forces[i];
				}
			}

		}


	}

}



