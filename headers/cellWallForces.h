#ifndef CELLWALL_FORCES_H
#define CELLWALL_FORCES_H

#include <math.h>
#include "cellWallObject.h"

/* Compute energy of glycosidic springs */
double compute_energy_gbond(CellWallMonolayer *cwl);

/* Compute energy of peptidic springs */
double compute_energy_pbond(CellWallMonolayer *cwl);

#endif