#include "cat_e_step.h"

void cat_e_step(double& dep_likelihood, double& indep_likelihood,
		double& tau, double* kappa, double &alpha,
		const int K, const int num_children, const t_setting* setting,
		const double children_sum, double* dirichlet_prior, double alpha_pi, double* digamma_sum_over_children, double* digamma_sum_over_children_for_kappa,
		int node_index)
{
	dep_likelihood = 0;
	indep_likelihood = 0;

	double kappa_tau_likelihood = 0;
	double kappa_tau_likelihood_old = 0;
	double kappa_tau_converged = 1;
	int kappa_tau_loop = 0;

	while ((kappa_tau_loop < 2) ||
			((kappa_tau_converged > setting->kappa_tau_converged) && (kappa_tau_loop < setting->kappa_tau_max_iter))) {
		kappa_tau_loop += 1;

		kappa_tau_likelihood = 0;
		opt_tau(tau,
				kappa,
				K, num_children,
				dirichlet_prior,
				alpha,
				node_index);

		kappa_tau_likelihood += opt_kappa(kappa, K, num_children,
				dirichlet_prior,
				alpha, tau,
				digamma_sum_over_children_for_kappa, node_index);

		for (int i = 0; i < K; ++i) {
			kappa_tau_likelihood += digamma(tau) * (-dirichlet_prior[i] + tau * kappa[i] + (alpha * kappa[i] - 1) * num_children);
		}
		kappa_tau_likelihood -= lgamma(tau);
		kappa_tau_likelihood -= alpha * num_children * (K - 1) / tau;

#ifdef _DEBUG

		printf("kappa_tau maximization %d: step %d \t nL %5.8f \t oL %5.8f \n",
				node_index, kappa_tau_loop, kappa_tau_likelihood, kappa_tau_likelihood_old);
#endif

		assert((kappa_tau_loop == 1) || (kappa_tau_likelihood >= kappa_tau_likelihood_old)
								|| (((kappa_tau_likelihood_old - kappa_tau_likelihood) / fabs(kappa_tau_likelihood_old) < 1e-2)
									&& (kappa_tau_loop >= 2)));

		if (0 != kappa_tau_likelihood_old && kappa_tau_likelihood < kappa_tau_likelihood_old) {
			printf("Warning: kappa_tau_likelihood is decreasing. node_index: %d \t step: %d \t old: %.8f \t new: %.8f n",
					node_index, kappa_tau_loop, kappa_tau_likelihood_old, kappa_tau_likelihood);
		}

		kappa_tau_converged = (kappa_tau_likelihood_old - kappa_tau_likelihood) / kappa_tau_likelihood_old;
		kappa_tau_likelihood_old = kappa_tau_likelihood;
	}

	if (kappa_tau_loop >= setting->kappa_tau_max_iter) {
		printf("kappa_tau_loop max reached\n");
		exit(-1);
	}

	if (setting->estimate_alpha) {
		indep_likelihood += opt_alpha(alpha, K, num_children, tau, kappa, digamma_sum_over_children, node_index);
	} else {
		double precompute = 0;
		double alpha_likelihood = 0;

		precompute = num_children * (digamma(tau) - (K - 1) / tau);
		for (int i = 0; i < K; ++i) {
			const double alphakappai = alpha * kappa[i];

			alpha_likelihood -= lgamma(alphakappai);
			precompute += num_children * kappa[i] * (log(kappa[i]) - digamma(tau * kappa[i]));
			precompute += kappa[i] * digamma_sum_over_children[i];
		}

		alpha_likelihood = num_children * (alpha_likelihood + lgamma(alpha)) + alpha * precompute;
		indep_likelihood += alpha_likelihood;
	}
	indep_likelihood += children_sum;

	const double digamma_tau = digamma(tau);

	indep_likelihood -= num_children * K * digamma_tau;

	for (int i = 0; i < K; ++i) {
		const double& kappai = kappa[i];
		const double taukappai = tau * kappai;
		const double digammataukappai = digamma(taukappai);
		const double common = (digammataukappai - digamma_tau);

		indep_likelihood -= num_children * (log(kappai) - digammataukappai);
		indep_likelihood += lgamma(taukappai);
		indep_likelihood -= taukappai * common;
		dep_likelihood += dirichlet_prior[i] * common;
	}

	indep_likelihood -= lgamma(tau);

	assert(!std::isnan(indep_likelihood));
	assert(!std::isnan(dep_likelihood));
}
