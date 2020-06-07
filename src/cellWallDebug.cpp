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

/*
 * Display glycosidic - glycosidic angles in a debug and sexy way
 * 
 * Parameters:
 *               pointor to cellwall object
 */
void display_glyco_glyco_angles(CellWallMonolayer *cwl){
  int i;
  int mi, mj;

  fprintf(stdout, "\n\t> Glyco-Glyco angles : %d ", cwl->get_total_gg_angles());
  
  for(i=0; i<cwl->get_total_gg_angles()*2; i+=2){
    mi = cwl->gg_angles[i];
    mj = cwl->gg_angles[i+1];
    fprintf(stdout, "\n\t\t> angle #%5d between bonds %5d ---^--- %5d", i/2, mi, mj);

    // because to get the masses index
    mi *= 2;
    mj *= 2;
    fprintf(stdout, "\n\t\t\t Involved masses %5d <-> %5d ^ ", int(cwl->glyco_bonds[mi]/3), int(cwl->glyco_bonds[mi+1]/3));
    fprintf(stdout, "%5d <-> %5d \n", int(cwl->glyco_bonds[mj]/3), int(cwl->glyco_bonds[mj+1]/3));
  }
  fprintf(stdout, "\n");

}

/*
 * Display lipid - lipid angles in a debug and sexy way
 * 
 * Parameters:
 *               pointor to lipid layer
 */
void display_lipid_lipid_angles(CellWallLipidLayer *ll){
  int i;
  int mi, mj, mk;

  fprintf(stdout, "\n\t> Lipid-Lipid angles : %d ", ll->get_total_lipid_lipid_angles());
  
  for(i=0; i<ll->get_total_lipid_lipid_angles()*3; i+=3){
    mj = int(ll->ll_angles[i])/3;
    mi = int(ll->ll_angles[i+1])/3;
    mk = int(ll->ll_angles[i+2])/3;
    
    fprintf(stdout, "\n\t\t> angle #%5d between lipid masses : %d -- %d -- %d\n", int(i/3), mj, mi, mk);

  }
  fprintf(stdout, "\n");

}

/*
 * Display lipid mesh debug and sexy way
 * 
 * Parameters:
 *               pointor to lipid layer
 */
void display_lipid_mesh(CellWallLipidLayer *ll){
  int i;
  int mi, mj, mk;

  fprintf(stdout, "\n\t> Lipid mesh : %d ", ll->get_total_number_of_mesh());
  
  for(i=0; i<ll->get_total_number_of_mesh()*3; i+=3){
    mj = int(ll->lipidic_mesh[i])/3;
    mi = int(ll->lipidic_mesh[i+1])/3;
    mk = int(ll->lipidic_mesh[i+2])/3;
    
    fprintf(stdout, "\n\t\t> mesh #%5d between lipid masses : %d -- %d -- %d\n", int(i/3), mj, mi, mk);

  }
  fprintf(stdout, "\n");

}