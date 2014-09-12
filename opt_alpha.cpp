#include "opt_alpha.h"

#define ALPHA_NEWTON_THRESH 1e-6
#define ALPHA_MAX_ITER 5000
#define DECREASE_ALLOWANCE 1e-8

double opt_alpha(double& alpha, const int ntopics, const int nchildren, const double tau, const double* kappa, const double* digamma_sum_over_children, const int node_index)
{
	int		iter = 0;
	double	d1 = 0;
	double	d2 = 0;

	double	likelihood = 0;
	double	old_likelihood = 0;

	double	init_alpha = 100;
	double	log_alpha = 0;

	double precompute = 0;

	precompute = nchildren * (digamma(tau) - (ntopics - 1) / tau);
	for (int i = 0; i < ntopics; ++i) {
		precompute += nchildren * kappa[i] * (log(kappa[i]) - digamma(tau * kappa[i]));
		precompute += kappa[i] * digamma_sum_over_children[i];
	}

	alpha = init_alpha;
	log_alpha = log(alpha);

	do {
		iter++;

//		likelihood = alpha * precompute + nchildren * lgamma(alpha);
//		d1 = precompute + nchildren * digamma(alpha);
//		d2 = nchildren * trigamma(alpha);
//
//		for (int i = 0; i < ntopics; ++i) {
//			const double alphakappai = alpha * kappa[i];
//
//			likelihood -= nchildren * lgamma(alphakappai);
//
//			d1 -= nchildren * kappa[i] * digamma(alphakappai);
//			d2 -= nchildren * kappa[i] * kappa[i] * trigamma(alphakappai);
//		}

		likelihood = 0;
		d1 = 0;
		d2 = 0;

		for (int i = 0; i < ntopics; ++i) {
			const double alphakappai = alpha * kappa[i];

			likelihood -= lgamma(alphakappai);

			d1 -= kappa[i] * digamma(alphakappai);
			d2 -= kappa[i] * kappa[i] * trigamma(alphakappai);
		}

		likelihood = nchildren * (likelihood +  lgamma(alpha))+ alpha * precompute;
		d1 = nchildren * (d1 + digamma(alpha)) + precompute;
		d2 = nchildren * (d2 + trigamma(alpha));

		assert(!std::isnan(likelihood));
		assert(!std::isnan(d1));
		assert(!std::isnan(d2));

#ifdef _DEBUG
		printf("alpha maximization %d : alpha %5.5f \t L %5.5f \t d1 %5.5f\n", node_index, alpha, likelihood, d1);
#endif
//		assert( (old_likelihood == 0) || (likelihood >= old_likelihood)
//			|| ( ((old_likelihood - likelihood) / fabs(old_likelihood) < DECREASE_ALLOWANCE)
//				&& (iter >= 5) )
//		);

		if (0 != old_likelihood && likelihood < old_likelihood) {
			printf("Warning: alpha_likelihood is decreasing. node_index: %d \t step: %d \t old: %.8f \t new: %.8f n", node_index, iter, old_likelihood, likelihood);
		}

		old_likelihood = likelihood;

		if (fabs(d1) < ALPHA_NEWTON_THRESH) {
			break;
		}

		log_alpha = log_alpha - d1 / (d2 * alpha + d1);
		alpha = exp(log_alpha);
		if (std::isnan(alpha) || alpha < 1e-10) {
			init_alpha = init_alpha * 10;
			printf("warning: alpha is nan; new init = %5.5f\n", init_alpha);
			alpha = init_alpha;
			log_alpha = log(alpha);
			old_likelihood = 0;
		}

	} while (iter < ALPHA_MAX_ITER);

	if (iter >= ALPHA_MAX_ITER) {
		printf("alpha iter max reached\n");
		exit(-1);
	}

	return likelihood;
}
