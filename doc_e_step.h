/*
 * doc_e_step.h
 *
 *  Created on: Feb 1, 2012
 *      Author: dok027
 */

#ifndef DOC_E_STEP_H_
#define DOC_E_STEP_H_

#include "common.h"
#include "datastructure.h"
#include "utils.h"

double doc_e_step(const t_document* doc, const double* dirichlet_prior, double* nu,
		double** digamma_lambda, double* digamma_lambda_sum, const t_setting* setting,
		const int doc_id, double** rho, double* old_rho);
#endif /* DOC_E_STEP_H_ */
