#include <stdio.h>

#include "cellWallConfig.h"
#include "cellWallUtils.h"

// Sexy and entertaining welcome message
void welcome_message(){

  fprintf(stdout, "\n    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    \n");
  fprintf(stdout, "   @----------------------------------@   \n");
  fprintf(stdout, "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@  \n");
  fprintf(stdout, " @---@                              @---@  \n");
  fprintf(stdout, "@---@  CW Bacteria Simulator %d.%d    @---@ \n", 
          cellWall_VERSION_MAJOR, cellWall_VERSION_MINOR);
  fprintf(stdout, "@---@       Gram + / Gram -          @---@\n");                           
  fprintf(stdout, " @---@                              @---@  \n");
  fprintf(stdout,  "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@   \n");
  fprintf(stdout,  "   @----------------------------------@    \n");
  fprintf(stdout,  "    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     \n");

}
