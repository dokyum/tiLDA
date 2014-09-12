#ifndef UTILS_H
#define UTILS_H

#include "datastructure.h"
#include "common.h"

double log_sum(double log_a, double log_b);
double trigamma(double x);
double digamma(double x);
double tetragamma(double x);
void setup_output_directory(t_setting* setting);

#endif
