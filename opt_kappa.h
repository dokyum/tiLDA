#ifndef OPT_KAPPA_H
#define OPT_KAPPA_H

#include "common.h"
#include "utils.h"
#include "model.h"

double opt_kappa(double* kappa, int ntopics, int nchildren,
		double* dirichlet_prior,
		double alpha, double tau,
		double* digamma_sum_over_children,
		int node_index);

#endif
