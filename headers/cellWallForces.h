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

#endif