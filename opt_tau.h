#ifndef OPT_TAU_H
#define OPT_TAU_H

#include "common.h"
#include "utils.h"

double opt_tau(double& tau, const double* kappa,
		const int ntopics, const int nchildren,
		const double* dirichlet_prior,
		const double alpha,
		const int node_index);

#endif
