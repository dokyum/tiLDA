#include "opt_tau.h"

#define TAU_NEWTON_THRESH 1e-6
#define TAU_MAX_ITER 5000

double opt_tau(double& tau, const double* kappa,
		const int ntopics, const int nchildren,
		const double* dirichlet_prior,
		const double alpha,
		const int node_index)
{
	int		iter = 0;
	double	d1 = 0;
	double	d2 = 0;

	double	init_tau = 100;
	double	log_tau = 0;

	double	likelihood = 0;
	double	old_likelihood = 0;

	tau = init_tau;
	log_tau = log(tau);
	do {
		iter++;

		d1 = 0;
		d2 = 0;
		likelihood = 0;

		double const common2 = nchildren * alpha * (ntopics - 1) / tau;
		double const trigammatau = trigamma(tau);
		double const digammatau = digamma(tau);
		for (int i = 0; i < ntopics; ++i) {
			double const taukappai = tau * kappa[i];
			double const trigammataukappai = trigamma(taukappai);
			double const common = dirichlet_prior[i] - taukappai + nchildren * (1 - alpha * kappa[i]);

			d1 += (trigammataukappai * kappa[i] - trigammatau) * common;
			d2 += kappa[i] * kappa[i] * ( tetragamma(taukappai) * common - trigammataukappai);
			d2 -= tetragamma(tau) * common;
			likelihood += (digamma(taukappai) - digammatau) * common
				+ lgamma(taukappai);
		}
		d1 += common2 / tau;
		d2 += trigammatau - 2 * common2 / tau / tau;

		likelihood -= lgamma(tau);
		likelihood -= common2;

		assert(!std::isnan(d1));
		assert(!std::isnan(d2));
		assert(!std::isnan(likelihood));

		assert( (old_likelihood == 0) || (likelihood >= old_likelihood) );
		if (0 != old_likelihood && likelihood < old_likelihood) {
			printf("Warning: tau_likelihood is decreasing. node_index: %d \t step: %d \t old: %.8f \t new: %.8f n",
					node_index, iter, old_likelihood, likelihood);
		}
		old_likelihood = likelihood;


#ifdef _DEBUG
		printf("tau maximization %d : tau %5.10f \t L %5.10f \t d1 %5.10f\n", node_index, tau, likelihood, d1);
#endif

		if (fabs(d1) < TAU_NEWTON_THRESH) {
			break;
		}

		log_tau = log_tau - d1 / (d2 * tau + d1);
		tau = exp(log_tau);
		if (std::isnan(tau) || tau < 1e-10) {
			init_tau = init_tau * 10;
			printf("warning: tau is nan; new init = %5.5f\n", init_tau);
			tau = init_tau;
			log_tau = log(tau);
			old_likelihood = 0;
		}
	} while (iter < TAU_MAX_ITER);

	if (iter >= TAU_MAX_ITER) {
		printf("tau iter max reached\n");
		exit(-1);
	}

	return likelihood;
}
