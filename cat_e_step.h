/*
 * cat_e_step.h
 *
 *  Created on: Feb 2, 2012
 *      Author: dok027
 */

#ifndef CAT_E_STEP_H_
#define CAT_E_STEP_H_

#include "common.h"
#include "utils.h"
#include "opt_kappa.h"
#include "opt_tau.h"
#include "opt_alpha.h"

void cat_e_step(double& dep_likelihood, double& indep_likelihood,
		double& tau, double* kappa, double &alpha,
		const int K, const int num_children, const t_setting* setting,
		const double children_sum, double* dirichlet_prior, double alpha_pi, double* digamma_sum_over_children, double* digamma_sum_over_children_for_kappa,
		int node_index);

#endif /* CAT_E_STEP_H_ */
