#ifndef CELLWALL_PARAMS_H
#define CELLWALL_PARAMS_H

const short int DIM=3; // 3Dimensions by default, 2D not possible 

const double PI=3.1415926535897932;

const double stiffness_g = 1000.0f; // pN/nm glycosidic
const double stiffness_p = 100.0f; // pN/nm peptidic
const double stiffness_l = 1.0f; // pN/nm lipidic

const double stiffness_gg = 0.1f; // pN/nm angle glyco-glyco
const double stiffness_pg = 0.01f; // pN/nm angle pepti-glyco
const double stiffness_ll = 0.1f; // pN/nm angle lipi-lipi

const double inner_pressure = 5.0f; // pN/nm^2

const double d0_g = 2.0f; // nm, distance at rest for glycosidic spring
const double d0_p = 1.0f; // nm, distance at rest for peptidic spring
const double alpha0_gg = PI; // rad, at rest for glyco-glyco angle
const double alpha0_pg = 0.5f * PI; // rad, at rest for pepti-glyco angle
const double alpha0_ll = PI; // rad, at rest for lipi-lipi angle

const double epsilon_g = 0.1f; // percentage of difference between rest and current
                               // state of spring when the structure is goemetrically 
                               // build

const double epsilon_lj = 0.001f; // pJ, epsilon for Lennard-Jones potential
const double cut_off = 2.5f; // nm cut off for LJ potential

const int output_iter = 200; // output data each # iter

const double fictive_temperature = 273.5f; // Kelvin for annealing
const double fictive_kb = 0.31415f; // fictive boltzman constant for annealing

#endif