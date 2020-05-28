#include "cellWallDebug.h"

/*
 * Display glycosidic bonds
 * 
 */
void display_glyco_bonds(CellWallMonolayer *cwl){
  int i;

  fprintf(stdout, "\n\t> Glycosidic springs : %d ", cwl->get_total_glycobonds());
  
  for(i=0; i<cwl->get_total_glycobonds()*2; i+=2){
    fprintf(stdout, "\n\t\t> bond #%5d : %5d <-----> %5d", i/2, cwl->glyco_bonds[i], cwl->glyco_bonds[i+1]);
  }
  fprintf(stdout, "\n");

}