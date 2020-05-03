#ifndef CELLWALL_PARAMS_H
#define CELLWALL_PARAMS_H

const short int DIM=3; // 3Dimensions by default, 2D not possible 

const double PI=3.1415926535897932;

const double stiffness_g = 5000.0; // pN/nm glycosidic
const double stiffness_p = 100.0; // pN/nm peptidic
const double stiffness_l = 1.0; // pN/nm lipidic

const double stiffness_gg = 0.1; // pN/nm angle glyco-glyco
const double stiffness_pg = 0.01; // pN/nm angle pepti-glyco
const double stiffness_ll = 0.05; // pN/nm angle lipi-lipi

const double inner_pressure = 0.5; // pN/nm^2

const double d0_g = 2.0; // nm, distance at rest for glycosidic spring
const double d0_p = 1.0; // nm, distance at rest for peptidic spring
const double alpha0_gg = PI; // rad, at rest for glyco-glyco angle
const double alpha0_pg = 0.5 * PI; // rad, at rest for pepti-glyco angle
const double alpha0_ll = PI; // rad, at rest for lipi-lipi angle

const double epsilon_lj = 0.01; // pJ, epsilon for Lennard-Jones potential
const double cut_off = 4.0; // nm cut off for LJ potential

#endif