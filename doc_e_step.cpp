/*
 * doc_e_step.cpp
 *
 *  Created on: Feb 1, 2012
 *      Author: dok027
 */

#include "doc_e_step.h"

#define DOC_DECREASE_ALLOWANCE 1e-6

extern double** 		digamma_nu;
extern double* 			digamma_nu_sum;
extern double			oneoverk;

double doc_e_step(const t_document* doc, const double* dirichlet_prior, double* nu,
		double** digamma_lambda, double* digamma_lambda_sum, const t_setting* setting,
		const int doc_id, double** rho, double* old_rho)
{
	const int& numtopics = setting->num_topics;
	const int& doclength = doc->length;
	const double doctotaloverk = (doc->total) / (double) (numtopics);

	// Initialize rho, nu
	for (int i = 0; i < numtopics; ++i) {
		for (int l = 0; l < doclength; ++l) {
			rho[l][i] = oneoverk;
		}
		nu[i] = dirichlet_prior[i] + doctotaloverk;
		digamma_nu[doc_id][i] = digamma(nu[i]);
	}

	int doc_loop = 0;
	double doc_likelihood = 0;
	double doc_likelihood_old = 0;
	double doc_converged = 1;
	double nu_sum = 0;
	double indep_part_likelihood = 0;
	double dep_part_likelihood = 0;

	while ((doc_loop < 2)
			|| ((doc_converged > setting->doc_converged) && (doc_loop < setting->doc_max_iter))) {
		doc_loop += 1;
		for (int l = 0; l < doclength; ++l) {
			double rhosum = 0;
			for (int i = 0; i < numtopics; ++i) {
				old_rho[i] = rho[l][i];
				rho[l][i] = digamma_nu[doc_id][i] + digamma_lambda[i][doc->words[l]] - digamma_lambda_sum[i];
				assert(rho[l][i] != 0);
				assert(!std::isnan(rho[l][i]));
				if (i > 0) {
					rhosum = log_sum(rhosum, rho[l][i]);
				} else {
					rhosum = rho[l][i];
				}
				assert(!std::isnan(rhosum));
			}
			for (int i = 0; i < numtopics; ++i) {
				rho[l][i] = exp(rho[l][i] - rhosum);
				nu[i] = nu[i] + (doc->counts[l]) * (rho[l][i] - old_rho[i]);
				// !!! a lot of extra digamma's here because of how we're computing it
				// !!! but its more automatically updated too.
				digamma_nu[doc_id][i] = digamma(nu[i]);
				assert(!std::isnan(digamma_nu[doc_id][i]));
			}
		}

		nu_sum = 0;
		for (int i = 0; i < numtopics; ++i) {
			nu_sum += nu[i];
		}
		digamma_nu_sum[doc_id] = digamma(nu_sum);

		indep_part_likelihood = -lgamma(nu_sum);
		dep_part_likelihood = 0;
		for (int i = 0; i < numtopics; ++i) {
			double delta = (digamma_nu[doc_id][i] - digamma_nu_sum[doc_id]);
			indep_part_likelihood += lgamma(nu[i]) - delta * nu[i];
			dep_part_likelihood += delta * dirichlet_prior[i];
			for (int l = 0; l < doclength; ++l) {
				if (rho[l][i] > 0) {
					indep_part_likelihood += rho[l][i] * (doc->counts[l])
							* (delta + digamma_lambda[i][doc->words[l]] - digamma_lambda_sum[i] - log(rho[l][i]));
				}
			}
		}
		assert(!std::isnan(indep_part_likelihood));
		assert(!std::isnan(dep_part_likelihood));
		doc_likelihood = indep_part_likelihood + dep_part_likelihood;

		doc_converged = (doc_likelihood_old - doc_likelihood) / doc_likelihood_old;
//		if (0 != doc_likelihood_old && doc_likelihood < doc_likelihood_old) {
//			printf("Warning: doc_likelihood is decreasing. doc_id: %d \t step: %d \t old: %.8f \t new: %.8f \t ratio: %.8f\n",
//					doc_id, doc_loop, doc_likelihood_old, doc_likelihood, doc_converged);
//		}

		assert((doc_loop == 1) || (doc_likelihood >= doc_likelihood_old)
				|| (((doc_likelihood_old - doc_likelihood) / fabs(doc_likelihood_old) < DOC_DECREASE_ALLOWANCE)
						&& (doc_loop >= 4)));


		doc_likelihood_old = doc_likelihood;
	}

	if (doc_loop >= setting->doc_max_iter) {
		printf("doc loop max reached %d\n", doc_id);
		exit(-1);
	}

	return indep_part_likelihood;
}
