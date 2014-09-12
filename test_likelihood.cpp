#include "test_likelihood.h"

extern double** 		digamma_nu;
extern double* 			digamma_nu_sum;
extern double 			oneoverk;
extern double**			digamma_lambda;
extern double*			digamma_lambda_sum;

void split_document(t_document*& inference_doc, t_document*& test_doc, const t_document* original_doc)
{
	inference_doc = (t_document*) malloc(sizeof(t_document));
	test_doc = (t_document*) malloc(sizeof(t_document));

	test_doc->length = original_doc->length / 2;
	inference_doc->length = original_doc->length - test_doc->length;

	test_doc->total = 0;
	test_doc->words = (int*) malloc(sizeof(int) * (test_doc->length));
	test_doc->counts = (int*) malloc(sizeof(int) * (test_doc->length));
	test_doc->parent_index = -1;

	inference_doc->total = 0;
	inference_doc->words = (int*) malloc(sizeof(int) * (inference_doc->length));
	inference_doc->counts = (int*) malloc(sizeof(int) * (inference_doc->length));
	inference_doc->parent_index = -1;

	for (int l = 0; l < original_doc->length; ++l) {
		int index_in_split_doc = l / 2;
		if (0 == l % 2) {
			inference_doc->words[index_in_split_doc] = original_doc->words[l];
			inference_doc->counts[index_in_split_doc] = original_doc->counts[l];
			inference_doc->total += original_doc->counts[l];
		} else {
			test_doc->words[index_in_split_doc] = original_doc->words[l];
			test_doc->counts[index_in_split_doc] = original_doc->counts[l];
			test_doc->total += original_doc->counts[l];
		}
	}
}

void free_document(t_document* doc)
{
	free(doc->words);
	free(doc->counts);
	free(doc);
}

void compute_lambda_statistics(t_tilda_var_model* model, double** expected_beta)
{
    for (int i = 0; i < model->num_topics; ++i) {
    	double lambda_sum = 0.0;

        for (int v = 0; v < model->num_terms; ++v) {
        	digamma_lambda[i][v] = digamma(model->lambda[i][v]);
        	lambda_sum += model->lambda[i][v];
        }
        digamma_lambda_sum[i] = digamma(lambda_sum);

        for (int v = 0; v < model->num_terms; ++v) {
        	expected_beta[i][v] = model->lambda[i][v] / lambda_sum;
        }
    }
}

void test_likelihood(t_setting* setting, const t_corpus* corpus, const std::vector<t_cat> tree_structure)
{
	FILE* 	fileptr_lowerbound_result;
	FILE* 	fileptr_lowerbound_summary;
	FILE* 	fileptr_document_completion_result;
	FILE* 	fileptr_document_completion_summary;
	char 	filename[MAX_BUF];

	sprintf(filename, "%s_tilda", setting->model_path);
	t_tilda_model* trained_tilda_model = load_tilda_model(filename);

	sprintf(filename, "%s_var", setting->model_path);
	t_tilda_var_model* trained_var_model = load_var_model(filename, corpus);

	double**	rho = NULL;
	double* 	old_rho = NULL;
	double* 	nu = NULL;
	double*		dirichlet_prior = NULL;
	double*		expected_theta = NULL;
	double**	expected_beta = NULL;


	const int& K = trained_tilda_model->num_topics;
	setting->num_topics = K;

	oneoverk = 1 / (double) K;

	double 	document_completion_sum_ll = 0.0;
	int		document_completion_sum_num_words = 0;

	double	lowerbound_sum_likelihood = 0;
	int  	lowerbound_sum_num_words = 0;

	nu = zero_init_double_array(K);
	rho = zero_init_double_matrix(corpus->max_length, K);
	old_rho = zero_init_double_array(K);
	dirichlet_prior = zero_init_double_array(K);
	expected_theta = zero_init_double_array(K);
	expected_beta = zero_init_double_matrix(K, corpus->num_terms);

	digamma_nu = zero_init_double_matrix(corpus->num_docs, K);
	digamma_nu_sum = zero_init_double_array(corpus->num_docs);

	digamma_lambda = zero_init_double_matrix(K, corpus->num_terms);
	digamma_lambda_sum = zero_init_double_array(K);

	compute_lambda_statistics(trained_var_model, expected_beta);

	sprintf(filename, "%s_lowerbound_result", setting->output_path);
	fileptr_lowerbound_result = fopen(filename, "w");

	sprintf(filename, "%s_document_completion_result", setting->output_path);
	fileptr_document_completion_result = fopen(filename, "w");

	for (int i = 0; i < tree_structure.size(); ++i) {
		const double&	alpha_t = trained_tilda_model->alpha[i];
		const double*	kappa_t = trained_var_model->kappa[i];
		const double&	tau_t = trained_var_model->tau[i];

		for (int j = 0; j < K; ++j) {
			dirichlet_prior[j] = alpha_t * kappa_t[j];
		}

		for (int d = 0; d < tree_structure[i].docids.size(); ++d) {
			const int& docid = tree_structure[i].docids[d];

			// evaluation using variational bound
			double this_doc_lowerbound = doc_e_step(&(corpus->docs[docid]), dirichlet_prior, nu,
													digamma_lambda, digamma_lambda_sum, setting,
													docid, rho, old_rho);

			assert(!std::isnan(this_doc_lowerbound));

			this_doc_lowerbound += lgamma(alpha_t);
			this_doc_lowerbound -= (K - alpha_t) * digamma(tau_t);
			this_doc_lowerbound -= alpha_t * (K - 1) / tau_t;
			for (int j = 0; j < K; ++j) {
				this_doc_lowerbound -= lgamma(alpha_t * kappa_t[j]) +
										(1 - alpha_t * kappa_t[j]) * (log(kappa_t[j]) - digamma(tau_t * kappa_t[j]));
			}

			for (int j = 0; j < K; ++j) {
				this_doc_lowerbound += dirichlet_prior[j] * (digamma_nu[docid][j] - digamma_nu_sum[docid]);
			}

			assert(!std::isnan(this_doc_lowerbound));

			fprintf(fileptr_lowerbound_result, "docid %d\tlower_bound %5.5f\tnum_words %d\n", docid, this_doc_lowerbound, corpus->docs[docid].total);

			lowerbound_sum_likelihood += this_doc_lowerbound;
			lowerbound_sum_num_words += corpus->docs[docid].total;

			// evaluation using document completion
			t_document*	inference_doc = NULL;
			t_document*	test_doc = NULL;
			split_document(inference_doc, test_doc, &(corpus->docs[docid]));
			double half_doc_lowerbound = doc_e_step(inference_doc, dirichlet_prior, nu,
													digamma_lambda, digamma_lambda_sum, setting,
													docid, rho, old_rho);

			assert(!std::isnan(half_doc_lowerbound));

			half_doc_lowerbound += lgamma(alpha_t);
			half_doc_lowerbound -= (K - alpha_t) * digamma(tau_t);
			half_doc_lowerbound -= alpha_t * (K - 1) / tau_t;
			for (int j = 0; j < K; ++j) {
				half_doc_lowerbound -= lgamma(alpha_t * kappa_t[j]) +
										(1 - alpha_t * kappa_t[j]) * (log(kappa_t[j]) - digamma(tau_t * kappa_t[j]));
			}

			for (int j = 0; j < K; ++j) {
				half_doc_lowerbound += dirichlet_prior[j] * (digamma_nu[docid][j] - digamma_nu_sum[docid]);
			}

			assert(!std::isnan(half_doc_lowerbound));

			double document_completion_log_likelihood = 0.0;
			double nu_sum = 0.0;
			for (int j = 0; j < K; ++j) {
				nu_sum += nu[j];
			}
			for (int j = 0; j < K; ++j) {
				expected_theta[j] = nu[j] / nu_sum;
			}

			for (int n = 0; n < test_doc->length; n++) {
				double this_word_likelihood = 0.0;
				for (int j = 0; j < K; ++j) {
					this_word_likelihood += expected_theta[j] * expected_beta[j][test_doc->words[n]];
				}
				document_completion_log_likelihood += log(this_word_likelihood + 1e-100) * test_doc->counts[n];
			}

			fprintf(fileptr_document_completion_result, "docid %d\thalf_lower_bound %5.5f\tscore %5.5f\ttest_num_words %d\n",
					docid, half_doc_lowerbound, document_completion_log_likelihood, test_doc->total);

			document_completion_sum_ll += document_completion_log_likelihood;
			document_completion_sum_num_words += test_doc->total;

			free_document(inference_doc);
			free_document(test_doc);
		}
	}

	fclose(fileptr_lowerbound_result);
	fclose(fileptr_document_completion_result);

	double perplexity = exp(-lowerbound_sum_likelihood / (double) lowerbound_sum_num_words);
	sprintf(filename, "%s_lowerbound_summary", setting->output_path);
	fileptr_lowerbound_summary = fopen(filename, "w");
	fprintf(fileptr_lowerbound_summary, "sum_lowerbound %5.5f\n", lowerbound_sum_likelihood);
	fprintf(fileptr_lowerbound_summary, "sum_num_words %d\n", lowerbound_sum_num_words);
	fprintf(fileptr_lowerbound_summary, "perplexity %5.5f\n", perplexity);
	fclose(fileptr_lowerbound_summary);

	double per_word_ll = document_completion_sum_ll / (double) document_completion_sum_num_words;
	sprintf(filename, "%s_document_completion_summary", setting->output_path);
	fileptr_document_completion_summary = fopen(filename, "w");
	fprintf(fileptr_document_completion_summary, "sum_num_words %d\n", document_completion_sum_num_words);
	fprintf(fileptr_document_completion_summary, "per_word_ll %5.5f\n", per_word_ll);
	fprintf(fileptr_document_completion_summary, "perplexity %5.5f\n", exp(-per_word_ll));
	fclose(fileptr_document_completion_summary);

	free_double_matrix(digamma_lambda);
	free(digamma_lambda_sum);

	free_double_matrix(digamma_nu);
	free(digamma_nu_sum);

	free_double_matrix(expected_beta);
	free(expected_theta);
	free(dirichlet_prior);
	free(nu);
	free_double_matrix(rho);
	free(old_rho);
	free_var_model(trained_var_model);
	free_tilda_model(trained_tilda_model);

}
