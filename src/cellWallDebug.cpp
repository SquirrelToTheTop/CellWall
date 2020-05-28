#include "cellWallDebug.h"

/*
 * Display glycosidic bonds in a debug and sexy way
 * 
 * Parameters:
 *               pointor to cellwall object
 */
void display_glyco_bonds(CellWallMonolayer *cwl){
  int i;

  fprintf(stdout, "\n\t> Glycosidic springs : %d ", cwl->get_total_glycobonds());
  
  for(i=0; i<cwl->get_total_glycobonds()*2; i+=2){
    fprintf(stdout, "\n\t\t> bond #%5d : %5d <-----> %5d", i/2, cwl->glyco_bonds[i], cwl->glyco_bonds[i+1]);
  }
  fprintf(stdout, "\n");

}

/*
 * Display peptidic bonds in a debug and sexy way
 * 
 * Parameters:
 *               pointor to cellwall object
 */
void display_pepti_bonds(CellWallMonolayer *cwl){
  int i;

  fprintf(stdout, "\n\t> Peptidic springs : %d ", cwl->get_total_peptibonds());
  
  for(i=0; i<cwl->get_total_peptibonds()*2; i+=2){
    fprintf(stdout, "\n\t\t> bond #%5d : %5d <-----> %5d", i/2, cwl->pepti_bonds[i], cwl->pepti_bonds[i+1]);
  }
  fprintf(stdout, "\n");

}