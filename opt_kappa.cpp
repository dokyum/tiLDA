#include "opt_kappa.h"

#define KAPPA_NEWTON_THRESH 1e-6
#define KAPPA_MAX_ITER 5000
#define LIKELIHOOD_DECREASE_ALLOWANCE 1e-5

extern double oneoverk;

double opt_kappa(double* kappa, int ntopics, int nchildren,
		double* dirichlet_prior,
		double alpha, double tau,
		double* digamma_sum_over_children,
		int node_index)
{
	double* g = NULL;
	double* h = NULL;
	double* delta_kappa = NULL;
	double* new_kappa = NULL;

	g = zero_init_double_array(ntopics);
	h = zero_init_double_array(ntopics);
	delta_kappa = zero_init_double_array(ntopics);
	new_kappa = zero_init_double_array(ntopics);

	for (int i = 0; i < ntopics; ++i) {
		kappa[i] = oneoverk;
	}

#ifdef _DEBUG
//	printf("kappa opt start %d : nchildren %d \t alpha %5.15f \t tau %5.15f \n",
//				node_index, nchildren, alpha, tau);
//	for (int i = 0; i < ntopics; ++i) {
//		printf("kappa opt start %d %d : dirichlet_prior %5.15f \t kappa %5.15f \n",
//				node_index, i, dirichlet_prior[i], kappa[i]);
//	}
#endif


	double	invhsum = 0;
	double	goverhsum = 0;
	double	coefficient = 0;
	double	old_likelihood = 0;
//	double	likelihood = 0;
	double	sqr_newton_decrement = 0;
	double	step_size;
	double	indep_new_likelihood = 0;
	double	dep_new_likelihood = 0;
	double	new_likelihood;
	double	expected_increase;
#ifdef _DEBUG
	double	initial_likelihood;
#endif

	int		iter = 0;

	for (int i = 0; i < ntopics; ++i) {
		double const taukappai = tau * kappa[i];
		double const alphakappai = alpha * kappa[i];
		double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
		double const digammataukappai = digamma(taukappai);
		double const logkappai = log(kappa[i]);
		dep_new_likelihood += digammataukappai * common;
		indep_new_likelihood -= nchildren * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
		indep_new_likelihood += alphakappai * digamma_sum_over_children[i];
		dep_new_likelihood += lgamma(taukappai);
	}
	new_likelihood = indep_new_likelihood + dep_new_likelihood;

#ifdef _DEBUG
	initial_likelihood = new_likelihood;
#endif

	do {
		iter++;

		invhsum = 0;
		goverhsum = 0;
		coefficient = 0;

		for (int i = 0; i < ntopics; ++i) {
			double const taukappai = tau * kappa[i];
			double const alphakappai = alpha * kappa[i];
			double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
			double const digammataukappai = digamma(taukappai);
			double const trigammataukappai = trigamma(taukappai);
			double const logkappai = log(kappa[i]);

			g[i] = tau * trigammataukappai * common
				- nchildren * alpha * (digamma(alphakappai) - logkappai + digammataukappai - 1)
				- nchildren / kappa[i]
				+ alpha * digamma_sum_over_children[i];

			h[i] = tau * tau * tetragamma(taukappai) * common
				- tau * trigammataukappai * (tau + 2 * alpha * nchildren)
				- alpha * alpha * trigamma(alphakappai) * nchildren
				+ alpha * nchildren / kappa[i]
				+ nchildren / (kappa[i] * kappa[i]);

			invhsum += 1 / h[i];
			goverhsum += g[i] / h[i];
		}

		old_likelihood = new_likelihood;

		coefficient = goverhsum / invhsum;
		sqr_newton_decrement = 0;
		expected_increase = 0;
		step_size = 1;
		for (int i = 0; i < ntopics; ++i) {
			delta_kappa[i] = (coefficient - g[i]) / h[i];
			sqr_newton_decrement -= h[i] * delta_kappa[i] * delta_kappa[i]; // this one is maximization
			expected_increase += g[i] * delta_kappa[i];
			//sqr_newton_decrement += g[i] * delta_kappa[i];
			if (delta_kappa[i] < 0) {
				double limit = (kappa[i] - 1e-10) / -(delta_kappa[i]);
				if (step_size > limit) {
					step_size = limit;
				}
			}
		}
#ifdef _DEBUG
		printf("kappa maximization %d : L %5.15f \t dL %5.15f \t indL %5.15f \t newton %5.15f \t %5.15f\n",
				node_index, old_likelihood, dep_new_likelihood, indep_new_likelihood, sqr_newton_decrement / 2, step_size);
//		for (int i = 0; i < ntopics; ++i) {
//			printf("kappa maximization %d %d: delta_kappa %5.15f \n",
//					node_index, i, delta_kappa[i]);
//		}
#endif

		if (sqr_newton_decrement < KAPPA_NEWTON_THRESH * 2 || step_size < 1e-8 ) {
			break;
		}

		// backtracking line search
		while(1) {
//			double sum_new_kappa = 0.0;
			indep_new_likelihood = 0.0;
			dep_new_likelihood = 0.0;

			for (int i = 0; i < ntopics; ++i) {
				new_kappa[i] = kappa[i] + step_size * delta_kappa[i];

				double const taukappai = tau * new_kappa[i];
				double const alphakappai = alpha * new_kappa[i];
				double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
				double const logkappai = log(new_kappa[i]);

				dep_new_likelihood += digamma(taukappai) * common;
				indep_new_likelihood -= nchildren * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
				indep_new_likelihood += alphakappai * digamma_sum_over_children[i];
				dep_new_likelihood += lgamma(taukappai);
			}
			new_likelihood = indep_new_likelihood + dep_new_likelihood;
#ifdef _DEBUG
			printf("line search  %d : nL: %5.15f \t bound: %5.15f \t step_size: %5.15f \n",
					node_index, new_likelihood, old_likelihood + 0.4 * step_size * expected_increase, step_size);
#endif

			if (new_likelihood > old_likelihood + 0.4 * step_size * expected_increase) {
//			if (new_likelihood > old_likelihood + 0.4 * step_size * sqr_newton_decrement) {
				break;
			}
			step_size *= 0.9;
			if (step_size < 1e-8) break;
		}
		if (step_size < 1e-8) break;

		for (int i = 0; i < ntopics; ++i) {
			kappa[i] = new_kappa[i];
			assert(!std::isnan(kappa[i]));
			assert(kappa[i] > 0);
		}

	} while (iter < KAPPA_MAX_ITER);

	if (iter >= KAPPA_MAX_ITER) {
		printf("KAPPA_MAX_ITER reached\n");
		exit(-1);
	}

#ifdef _DEBUG
	printf("%d TL %5.15f \t IL %5.15f \n", node_index, new_likelihood, initial_likelihood);
	assert(new_likelihood >= initial_likelihood ||
		( ((initial_likelihood - new_likelihood) / fabs(initial_likelihood) < LIKELIHOOD_DECREASE_ALLOWANCE)
		 && (iter >= 3) )
		);
#endif

	free(new_kappa);
	free(delta_kappa);
	free(g);
	free(h);

	return new_likelihood;
}
