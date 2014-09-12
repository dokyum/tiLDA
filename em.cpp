#include "em.h"

#define DECREASE_ALLOWANCE 1e-5
#define MAX_THREADS	32

double** 	digamma_nu = NULL;
double* 	digamma_nu_sum = NULL;
double**	dirichlet_prior = NULL;
double*		dirichlet_prior_root = NULL;
double** 	digamma_sum_over_children = NULL;
double**	digamma_sum_over_children_for_kappa = NULL;
double**	digamma_lambda = NULL;
double*		digamma_lambda_sum = NULL;
double*		lgamma_lambda_sum = NULL;

double 			oneoverk;
double			etaoverv;

// controls for threads
pthread_attr_t 	attr;
pthread_t 		thread_handle[MAX_THREADS];
bool			finalizeThread = false;

// variables that are accessed by threads.
FILE* 						likelihood_file = NULL;
t_corpus const *			thread_use_corpus = NULL;
t_setting const *			thread_use_setting = NULL;
t_tilda_model* 				thread_use_tilda_model = NULL;
t_tilda_var_model* 			thread_use_var_model = NULL;
t_tilda_suffstats* 			thread_use_suffstats = NULL;
std::vector<t_cat>			thread_use_tree_structure;
int*						thread_use_waiting_children = NULL;
double*						thread_use_returned_likelihood_from_doc = NULL;
double*						thread_use_returned_likelihood_from_cat = NULL;
double*						thread_use_cat_old_likelihood = NULL;
int*						thread_use_cat_em_step = NULL;
double						whole_likelihood_old = 0;

void push_task(int node_index, t_task_mode mode, int parent_index)
{
	t_task to_push_task;
	to_push_task.node_index = node_index;
	to_push_task.mode = mode;
	pthread_mutex_lock(&mutex_for_queue);
	task_queue.push(to_push_task);
	if (parent_index != -1) {
		thread_use_waiting_children[parent_index]++;
	}
	pthread_mutex_unlock(&mutex_for_queue);
}

void push_task_wo_lock(int node_index, t_task_mode mode, int parent_index)
{
	t_task to_push_task;
	to_push_task.node_index = node_index;
	to_push_task.mode = mode;
	task_queue.push(to_push_task);
	if (parent_index != -1) {
		thread_use_waiting_children[parent_index]++;
	}
}

void push_children(int node_index)
{
	pthread_mutex_lock(&mutex_for_queue);
	for (unsigned int i = 0; i < thread_use_tree_structure[node_index].catids.size(); ++i) {
		push_task_wo_lock(thread_use_tree_structure[node_index].catids[i], PRE, node_index);
	}
	for (unsigned int i = 0; i < thread_use_tree_structure[node_index].docids.size(); ++i) {
		push_task_wo_lock(thread_use_tree_structure[node_index].docids[i], DOC, node_index);
	}
	pthread_mutex_unlock(&mutex_for_queue);
}

void decrease_waiting_children(int parent_index)
{
	pthread_mutex_lock(&mutex_for_queue);

	if (thread_use_waiting_children[parent_index] <= 0) {
		printf("Error: Wrong count of waiting children. %d has %d\n", parent_index, thread_use_waiting_children[parent_index]);
		exit(-1);
	}

	thread_use_waiting_children[parent_index]--;
	if (thread_use_waiting_children[parent_index] == 0) {
		push_task_wo_lock(parent_index, POST, -1);
	}
	pthread_mutex_unlock(&mutex_for_queue);
}

void *my_thread_function(void *lpParam) {
	int			thread_index = *((int*)lpParam);
	double* 	old_rho = NULL;
	const int&	K = thread_use_tilda_model->num_topics;
	time_t 		rawtime;

	old_rho = zero_init_double_array(K);

	while(!finalizeThread) {
		pthread_mutex_lock(&mutex_for_queue);
		if (task_queue.empty()) {
//			printf("queue empty -> yield thread %d \n", thread_index);
			pthread_mutex_unlock(&mutex_for_queue);
			pthread_yield();
		} else {
			t_task atask = task_queue.front();
			task_queue.pop();
			pthread_mutex_unlock(&mutex_for_queue);

			if (atask.mode == PRE) {
				#ifdef _DEBUG
				printf("pre %d at thread %d \n", atask.node_index, thread_index);
				#endif

				thread_use_cat_em_step[atask.node_index] = 0;
				thread_use_cat_old_likelihood[atask.node_index] = 0.0;

				thread_use_tilda_model->alpha[atask.node_index] = thread_use_setting->initial_alpha;
				thread_use_var_model->tau[atask.node_index] = (double) K;
				for (int i = 0; i < K; ++i) {
					thread_use_var_model->kappa[atask.node_index][i] = oneoverk;
					dirichlet_prior[atask.node_index][i] = (thread_use_setting->initial_alpha) * oneoverk;
				}

                push_children(atask.node_index);
			} else if (atask.mode == POST) {
				#ifdef _DEBUG
				printf("post %d at thread %d \n", atask.node_index, thread_index);
				#endif

				int&		em_step = thread_use_cat_em_step[atask.node_index];
				double&		old_node_likelihood = thread_use_cat_old_likelihood[atask.node_index];
				double 		node_likelihood = 0;
				double 		node_indep_likelihood = 0;
				double 		node_dep_likelihood = 0;
				double 		node_converged = 0;

				const int 	num_children = thread_use_tree_structure[atask.node_index].catids.size()
						+ thread_use_tree_structure[atask.node_index].docids.size();
				double 		children_sum = 0.0;
				const int& 	parent_index = thread_use_tree_structure[atask.node_index].parent_index;

				em_step++;

				memset(digamma_sum_over_children[atask.node_index], 0, sizeof(double) * K);
				memset(digamma_sum_over_children_for_kappa[atask.node_index], 0, sizeof(double) * K);

				for (unsigned int j = 0; j < thread_use_tree_structure[atask.node_index].catids.size(); ++j) {
					const int& c = thread_use_tree_structure[atask.node_index].catids[j];
					const double digamma_tau = digamma(thread_use_var_model->tau[c]);

					children_sum += thread_use_returned_likelihood_from_cat[c];
					for (int i = 0; i < K; ++i) {
						const double digamma_taukappai = digamma(thread_use_var_model->tau[c] * thread_use_var_model->kappa[c][i]);
						digamma_sum_over_children_for_kappa[atask.node_index][i] += digamma_taukappai;
						digamma_sum_over_children[atask.node_index][i] += digamma_taukappai - digamma_tau;
					}
				}
				for (unsigned int j = 0; j < thread_use_tree_structure[atask.node_index].docids.size(); ++j) {
					const int& d = thread_use_tree_structure[atask.node_index].docids[j];
					children_sum += thread_use_returned_likelihood_from_doc[d];
					for (int i = 0; i < K; ++i) {
						digamma_sum_over_children_for_kappa[atask.node_index][i] += digamma_nu[d][i];
						digamma_sum_over_children[atask.node_index][i] += digamma_nu[d][i] - digamma_nu_sum[d];
					}
				}

				if (0 < atask.node_index) { // for other than root
					cat_e_step(node_dep_likelihood, node_indep_likelihood,
							thread_use_var_model->tau[atask.node_index],
							thread_use_var_model->kappa[atask.node_index],
							thread_use_tilda_model->alpha[atask.node_index],
							K, num_children, thread_use_setting,
							children_sum, dirichlet_prior[parent_index], thread_use_tilda_model->alpha[parent_index],
							digamma_sum_over_children[atask.node_index], digamma_sum_over_children_for_kappa[atask.node_index], atask.node_index);
				} else { // for root
					cat_e_step(node_dep_likelihood, node_indep_likelihood,
							thread_use_var_model->tau[atask.node_index],
							thread_use_var_model->kappa[atask.node_index],
							thread_use_tilda_model->alpha[atask.node_index],
							K, num_children, thread_use_setting,
							children_sum, dirichlet_prior_root, thread_use_tilda_model->gamma,
							digamma_sum_over_children[atask.node_index], digamma_sum_over_children_for_kappa[atask.node_index], atask.node_index);
				}

				node_likelihood = node_dep_likelihood + node_indep_likelihood;
				node_converged = (old_node_likelihood - node_likelihood) / old_node_likelihood;

				if (0 != old_node_likelihood && node_likelihood < old_node_likelihood) {
					printf("Warning: node_likelihood is decreasing. node_id: %d \t step: %d \t old: %.8f \t new: %.8f \t ratio: %.8f\n",
							atask.node_index, em_step, old_node_likelihood, node_likelihood, node_converged);
				}

#ifdef _DEBUG
				time(&rawtime);
				printf("node %d loop %d\t L %.10f \t oL %.10f \t conv %.10f \t %s\n",
						atask.node_index, em_step, node_likelihood, old_node_likelihood, node_converged, ctime(&rawtime));
#endif

//				assert( (em_step == 1)
//						|| (node_likelihood >= old_node_likelihood)
//						|| ((-node_converged < DECREASE_ALLOWANCE) && (em_step >= 20) ) );


				if (0 < atask.node_index) { // for other than root
					if (em_step > 1 && node_converged < thread_use_setting->cat_converged) { // converged
						time(&rawtime);
						printf("node %d converged at loop %d with %.10f %s\n", atask.node_index, em_step, node_likelihood, ctime(&rawtime));
						thread_use_returned_likelihood_from_cat[atask.node_index] = node_indep_likelihood;
						decrease_waiting_children(parent_index);
					} else if (em_step > thread_use_setting->cat_max_iter) { // error
						printf("Error: Cat loop max reached %d\n", atask.node_index);
						exit(-1);
					} else { // loop. retry children
						old_node_likelihood = node_likelihood;

						// update dirichlet_prior
						for (int i = 0; i < K; ++i) {
							dirichlet_prior[atask.node_index][i] = thread_use_tilda_model->alpha[atask.node_index] * thread_use_var_model->kappa[atask.node_index][i];
						}

		                push_children(atask.node_index);
					}

				} else { // for root
					// At the moment, we do not estimate gamma and eta.
					double whole_likelihood = 0;
					double whole_converged = 0;

					if (!task_queue.empty()) {
						printf("Error: There are remaining tasks at the end of routine for the root.\n");
						exit(-1);
					}

					collect_lambda_ss(thread_use_suffstats, thread_use_var_model, thread_use_corpus);
					whole_likelihood += node_likelihood + opt_lambda(thread_use_suffstats, thread_use_var_model);
					whole_converged = (whole_likelihood_old - whole_likelihood) / whole_likelihood_old;

					if (0 != whole_likelihood_old && whole_likelihood < whole_likelihood_old) {
						printf("Warning: whole_likelihood is decreasing. step: %d \t old: %.8f \t new: %.8f \t ratio: %.8f\n",
								em_step, whole_likelihood_old, whole_likelihood, whole_converged);
					}

					// output model and likelihood
					fprintf(likelihood_file, "END STEP %d \t L %10.10f \t ratio %5.5e\n", em_step, whole_likelihood, whole_converged);
					fflush(likelihood_file);

					if ((em_step % thread_use_setting->model_save_freq) == 0) {
						char 	filename[MAX_BUF];
						sprintf(filename, "%s/%06d_tilda", thread_use_setting->output_path, em_step);
						save_tilda_model(thread_use_tilda_model, filename);
						sprintf(filename, "%s/%06d_var", thread_use_setting->output_path, em_step);
						save_var_model(thread_use_var_model, filename, thread_use_corpus);
					}

					if (em_step > 1 && whole_converged < thread_use_setting->em_converged) { // converged
						printf("set final %s\n", ctime(&rawtime));
						finalizeThread = true;
						break;
					} else if (em_step > thread_use_setting->em_max_iter) { // error
						printf("em loop max reached %d\n", atask.node_index);
						exit(-1);
					} else { // loop. retry children
						old_node_likelihood = node_likelihood;
						whole_likelihood_old = whole_likelihood;

						// update dirichlet_prior
						for (int i = 0; i < K; ++i) {
							dirichlet_prior[atask.node_index][i] = thread_use_tilda_model->alpha[atask.node_index] * thread_use_var_model->kappa[atask.node_index][i];
						}

		                push_children(atask.node_index);
					}

				}
			} else if (atask.mode == DOC) {
				//printf("doc %d at thread %d \n", atask.node_index, thread_index);

				const int& parent_index = thread_use_corpus->docs[atask.node_index].parent_index;

				// procedure for a doc. converges here.
				thread_use_returned_likelihood_from_doc[atask.node_index] =
						doc_e_step(&(thread_use_corpus->docs[atask.node_index]), dirichlet_prior[parent_index],
								thread_use_var_model->nu[atask.node_index],
								digamma_lambda, digamma_lambda_sum, thread_use_setting,
								atask.node_index, thread_use_var_model->rho[atask.node_index], old_rho);
#ifdef _DEBUG
				//printf("doc %d L %.10f \n", atask.node_index, thread_use_returned_likelihood_from_doc[atask.node_index]);
#endif

				decrease_waiting_children(parent_index);
			} else {
				assert(0);
			}
		}
	}
	free(old_rho);

	printf("thread %d fin\n", thread_index);

	pthread_exit(NULL);
}

void random_initialize_lambda(t_tilda_var_model* model)
{
    int& num_topics = model->num_topics;
    int& num_terms = model->num_terms;
    int& num_docs = model->num_docs;
    double denominator = (double) (num_topics * num_terms);

	Generator	g;
	Gamma_Dist	gamma_dist(1.0, 1.0);
	g.seed(static_cast<unsigned int>(std::time(0)));

    for (int i = 0; i < num_topics; ++i) {
    	double lambda_sum = 0.0;
        for (int v = 0; v < num_terms; ++v) {
        	model->lambda[i][v] = gamma_dist(g) * num_docs * 100 / denominator;
        	digamma_lambda[i][v] = digamma(model->lambda[i][v]);
        	lambda_sum += model->lambda[i][v];
        }
        digamma_lambda_sum[i] = digamma(lambda_sum);
        lgamma_lambda_sum[i] = lgamma(lambda_sum);
    }
}

double opt_lambda(t_tilda_suffstats* ss, t_tilda_var_model* model)
{
	double lambda_likelihood = 0;

    for (int i = 0; i < model->num_topics; ++i) {
    	double lambda_sum = 0.0;

    	lambda_likelihood -= lgamma_lambda_sum[i];

        for (int v = 0; v < model->num_terms; ++v) {
        	// compute likelihood before updating lambda
        	lambda_likelihood += lgamma(model->lambda[i][v]);
        	lambda_likelihood += (etaoverv - model->lambda[i][v]) * (digamma_lambda[i][v] - digamma_lambda_sum[i]);

        	model->lambda[i][v] = etaoverv + ss->class_word[i][v];
        	digamma_lambda[i][v] = digamma(model->lambda[i][v]);
        	lambda_sum += model->lambda[i][v];
        }
        digamma_lambda_sum[i] = digamma(lambda_sum);
        lgamma_lambda_sum[i] = lgamma(lambda_sum);
    }

    return lambda_likelihood;
}

void run_em(t_setting* setting, t_corpus* corpus, const std::vector<t_cat> tree_structure)
{
	t_tilda_model* 		tilda_model = NULL;
	t_tilda_var_model* 	var_model = NULL;
	t_tilda_suffstats* 	suffstats = NULL;

	int 	thread_param[MAX_THREADS];
	char 	filename[MAX_BUF];
	int 	rc;

	tilda_model = new_tilda_model(corpus->num_terms, setting->num_topics, tree_structure.size());
	var_model = new_var_model(tilda_model->num_cats, corpus->num_docs, setting->num_topics, corpus);
	suffstats = new_tilda_suffstats(tilda_model);
	digamma_lambda = zero_init_double_matrix(setting->num_topics, corpus->num_terms);
	digamma_lambda_sum = zero_init_double_array(setting->num_topics);
	lgamma_lambda_sum = zero_init_double_array(setting->num_topics);

	oneoverk = 1 / (double) (setting->num_topics);

	if (!setting->warm_start) {
		tilda_model->gamma = setting->initial_gamma;
		tilda_model->eta = setting->initial_eta;
		etaoverv = tilda_model->eta / (double) (corpus->num_terms);

		if (setting->corpus_init) {
			corpus_initialize_ss(suffstats, tilda_model, corpus, setting->num_docs_for_init);
			opt_lambda(suffstats, var_model);
		} else {
			random_initialize_lambda(var_model);
		}
	} else {
		tilda_model->gamma = setting->initial_gamma;
		tilda_model->eta = setting->initial_eta;
		etaoverv = tilda_model->eta / (double) (corpus->num_terms);

		warm_start_var_model(var_model, setting->warm_start_path);
	}

	printf("Setting #topics: %d \t output: %s \n ", setting->num_topics, setting->output_path);

	digamma_nu = zero_init_double_matrix(corpus->num_docs, setting->num_topics);
	digamma_nu_sum = zero_init_double_array(corpus->num_docs);
	dirichlet_prior = zero_init_double_matrix(tree_structure.size(), setting->num_topics);
	dirichlet_prior_root = zero_init_double_array(setting->num_topics);
	digamma_sum_over_children = zero_init_double_matrix(tree_structure.size(), setting->num_topics);
	digamma_sum_over_children_for_kappa = zero_init_double_matrix(tree_structure.size(), setting->num_topics);


	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	sprintf(filename, "%s/000000_tilda", setting->output_path);
	save_tilda_model(tilda_model, filename);

	sprintf(filename, "%s/likelihood.dat", setting->output_path);
	likelihood_file = fopen(filename, "w");


	pthread_mutex_init(&mutex_for_queue, NULL);

	for (int k = 0; k < setting->num_topics; ++k) {
		dirichlet_prior_root[k] = tilda_model->gamma / (double) (setting->num_topics);
	}
	push_task(0, PRE, -1);

	thread_use_corpus = corpus;
	thread_use_setting = setting;
	thread_use_tilda_model = tilda_model;
	thread_use_var_model = var_model;
	thread_use_suffstats = suffstats;
	thread_use_tree_structure = tree_structure;
	thread_use_waiting_children = (int*) malloc(sizeof(int) * tree_structure.size());
	memset(thread_use_waiting_children, 0, sizeof(int) * tree_structure.size());
	thread_use_returned_likelihood_from_doc = zero_init_double_array(corpus->num_docs);
	thread_use_returned_likelihood_from_cat = zero_init_double_array(tree_structure.size());
	thread_use_cat_old_likelihood = zero_init_double_array(tree_structure.size());
	thread_use_cat_em_step = (int*) malloc(sizeof(int) * tree_structure.size());
	memset(thread_use_cat_em_step, 0, sizeof(int) * tree_structure.size());

	for (int i = 0; i < setting->num_threads; ++i) {
		thread_param[i] = i;

		rc = pthread_create(&thread_handle[i], &attr, my_thread_function, (void*) &(thread_param[i]));
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	for (int i = 0; i < setting->num_threads; ++i) {
		void *status;
		pthread_join(thread_handle[i], &status);
	}

	// output the final model
	sprintf(filename, "%s/final_tilda", setting->output_path);
	save_tilda_model(tilda_model, filename);
	sprintf(filename, "%s/final_var", setting->output_path);
	save_var_model(var_model, filename, corpus);

	fclose(likelihood_file);

	free(thread_use_cat_em_step);
	free(thread_use_cat_old_likelihood);
	free(thread_use_returned_likelihood_from_doc);
	free(thread_use_returned_likelihood_from_cat);
	free(thread_use_waiting_children);

	free_tilda_suffstats(suffstats);
	free_var_model(var_model);
	free_tilda_model(tilda_model);

	free(digamma_nu_sum);
	free_double_matrix(digamma_nu);
	free_double_matrix(digamma_sum_over_children);
	free_double_matrix(digamma_sum_over_children_for_kappa);
	free(dirichlet_prior_root);
	free_double_matrix(dirichlet_prior);

	free(lgamma_lambda_sum);
	free(digamma_lambda_sum);
	free_double_matrix(digamma_lambda);

	pthread_mutex_destroy(&mutex_for_queue);
	pthread_attr_destroy(&attr);
}

