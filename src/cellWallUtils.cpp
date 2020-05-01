#include <stdio.h>
#include "cellWallUtils.h"

void welcome_message(){

  fprintf(stdout, "\n    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    \n");
  fprintf(stdout, "   @--------------------------------@   \n");
  fprintf(stdout, "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@  \n");
  fprintf(stdout, " @---@                            @---@  \n");
  fprintf(stdout, "@---@    CW Bacteria Simulator     @---@ \n");
  fprintf(stdout, "@---@       Gram + / Gram -        @---@\n");                           
  fprintf(stdout, " @---@                            @---@  \n");
  fprintf(stdout,  "  @---@@@@@@@@@@@@@@@@@@@@@@@@@@@@---@   \n");
  fprintf(stdout,  "   @--------------------------------@    \n");
  fprintf(stdout,  "    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     \n");

}