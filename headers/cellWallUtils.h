#ifndef CELLWALL_UTILS_H
#define CELLWALL_UTILS_H

#include <math.h>

const double MBYTES = 1024.0f*1024.0f;

// display a lovely introducing message (ascii art)
void welcome_message();

double norm_vect(double *, int);
#endif
