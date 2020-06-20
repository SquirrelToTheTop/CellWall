#ifndef CELLWALL_FORCES_H
#define CELLWALL_FORCES_H

#include <math.h>
#include "cellWallObject.h"
#include "cellWallLipidLayer.h"

/* Compute energy of glycosidic springs */
double compute_energy_gbond(CellWallMonolayer *cwl);

/* Compute energy of peptidic springs */
double compute_energy_pbond(CellWallMonolayer *cwl);

/* Compute energy of lipidic springs */
double compute_energy_lbond(CellWallLipidLayer *ll);

/* Compute energy of bending for glycosidic - glycosidic */
double compute_energy_gg_angles(CellWallMonolayer *cwl);

/* Compute energy of bending for glycosidic - peptidic */
double compute_energy_gp_angles(CellWallMonolayer *cwl);

/* Compute energy of bending for lipid - lipid */
double compute_energy_ll_angles(CellWallLipidLayer *ll);

/* Compute energy of pressure in/on the lipid layer */
double compute_energy_pressure(CellWallLipidLayer *ll);

/* Compute energy of Lennard-Jones interaction between cellwall and lipid layer */
double compute_energy_lennardJones(CellWallMonolayer *cwl, CellWallLipidLayer *ll);

/* Compute force of pressure on lipid mass */
void compute_force_pressure(CellWallLipidLayer *ll);

#endif