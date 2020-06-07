#ifndef CELLWALL_DEBUG_H
#define CELLWALL_DEBUG_H

#include <stdio.h>
#include "cellWallObject.h"
#include "cellWallLipidLayer.h"

// display glycosidic bonds for debug 
void display_glyco_bonds(CellWallMonolayer *cwl);

// display peptidic bonds for debug 
void display_pepti_bonds(CellWallMonolayer *cwl);

// display glycosidic - glycosidic angles and masses number for debug
void display_glyco_glyco_angles(CellWallMonolayer *cwl);

// display lipid - lipid angles and masses number for debug
void display_lipid_lipid_angles(CellWallLipidLayer *ll);

// display lipid - lipid mesh and masses number for debug
void display_lipid_mesh(CellWallLipidLayer *ll);

#endif