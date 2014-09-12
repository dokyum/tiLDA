#ifndef OPT_ALPHA_H
#define OPT_ALPHA_H

#include "common.h"
#include "utils.h"

double opt_alpha(double& alpha, const int ntopics, const int nchildren, const double tau, const double* kappa, const double* digamma_sum_over_children, const int node_index);

#endif
