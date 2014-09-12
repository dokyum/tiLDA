#include "common.h"
#include "em.h"
#include "cokus.h"
#include "data.h"
#include "utils.h"
#include "test_likelihood.h"

void read_input(char* input_path, t_setting* setting)
{
    FILE*	fileptr;

    char	mode[MAX_BUF];
	int		temp;

    printf("loading %s\n", input_path);
    fileptr = fopen(input_path, "r");
	fscanf(fileptr, "MODE %s\n", mode);
	if (strcmp(mode, "EST") == 0) {
		setting->mode = EST;
	} else if (strcmp(mode, "INF") == 0) {
		setting->mode = INF;
	} else {
		printf("%s\n", mode);
		assert(0);
	}

	if (setting->mode == EST) {
		fscanf(fileptr, "initial_gamma %lf\n", &(setting->initial_gamma));
		fscanf(fileptr, "initial_alpha %lf\n", &(setting->initial_alpha));
		fscanf(fileptr, "initial_eta %lf\n", &(setting->initial_eta));

		fscanf(fileptr, "estimate_alpha %d\n", &(temp));
		if (temp > 0) {
			setting->estimate_alpha = true;
		} else {
			setting->estimate_alpha = false;
		}

		fscanf(fileptr, "num_topics %d\n", &(setting->num_topics));

		fscanf(fileptr, "corpus_path %s\n", setting->corpus_path);
		fscanf(fileptr, "tree_structure_path %s\n", setting->tree_structure_path);
		fscanf(fileptr, "node_to_docids_path %s\n", setting->node_to_docids_path);
		fscanf(fileptr, "output_path %s\n", setting->output_path);

		fscanf(fileptr, "num_threads %d\n", &(setting->num_threads));

		fscanf(fileptr, "em_max_iter %d\n", &(setting->em_max_iter));
		fscanf(fileptr, "cat_max_iter %d\n", &(setting->cat_max_iter));
		fscanf(fileptr, "kappa_tau_max_iter %d\n", &(setting->kappa_tau_max_iter));
		fscanf(fileptr, "doc_max_iter %d\n", &(setting->doc_max_iter));

		fscanf(fileptr, "em_converged %lf\n", &(setting->em_converged));
		fscanf(fileptr, "cat_converged %lf\n",  &(setting->cat_converged));
		fscanf(fileptr, "kappa_tau_converged %lf\n",  &(setting->kappa_tau_converged));
		fscanf(fileptr, "doc_converged %lf\n", &(setting->doc_converged));

		fscanf(fileptr, "model_save_freq %d\n", &(setting->model_save_freq));
		fscanf(fileptr, "random_seed %d\n", &(setting->random_seed));

		fscanf(fileptr, "corpus_init %d\n", &(temp));
		if (temp > 0) {
			setting->corpus_init = true;
		} else {
			setting->corpus_init = false;
		}
		fscanf(fileptr, "num_docs_for_init %d\n", &(setting->num_docs_for_init));

		fscanf(fileptr, "warm_start %d\n", &(temp));
		if (temp > 0) {
			setting->warm_start = true;
		} else {
			setting->warm_start = false;
		}
		fscanf(fileptr, "warm_start_path %s\n", setting->warm_start_path);

	} else if (setting->mode == INF) {
		fscanf(fileptr, "model_path %s\n", setting->model_path);
		fscanf(fileptr, "corpus_path %s\n", setting->corpus_path);
		fscanf(fileptr, "tree_structure_path %s\n", setting->tree_structure_path);
		fscanf(fileptr, "node_to_docids_path %s\n", setting->node_to_docids_path);
		fscanf(fileptr, "output_path %s\n", setting->output_path);
		fscanf(fileptr, "doc_max_iter %d\n", &(setting->doc_max_iter));
		fscanf(fileptr, "doc_converged %lf\n", &(setting->doc_converged));
		setting->random_seed = -1;
	} else {
		assert(0);
	}
    fclose(fileptr);
}

int main(int argc, char* argv[])
{
	t_corpus*			corpus = NULL;
	t_setting*			setting = NULL;
	std::vector<t_cat>	tree_structure;

	if (argc != 2) {
		printf("usage: only one argv\n");
		return -1;
	}

	setting = (t_setting*) malloc(sizeof(t_setting));
	read_input(argv[1], setting);

	if (setting->random_seed < 0) {
	    time_t t1;
	    (void) time(&t1);
        setting->random_seed = t1;
	}
	seedMT(setting->random_seed);


	if (setting->mode == EST) {
		corpus = read_data(setting->corpus_path);
		tree_structure = read_tree_structure(setting->tree_structure_path);
		read_nodeid_to_docids(tree_structure, setting->node_to_docids_path, corpus);

		for (unsigned int i = 0; i < tree_structure.size(); ++i) {
			assert(tree_structure[i].catids.size() > 0 || tree_structure[i].docids.size() > 0);
		}

        setup_output_directory(setting);
		run_em(setting, corpus, tree_structure);
	} else if (setting->mode == INF) {
		corpus = read_data(setting->corpus_path);
		tree_structure = read_tree_structure(setting->tree_structure_path);
		read_nodeid_to_docids(tree_structure, setting->node_to_docids_path, corpus);

		test_likelihood(setting, corpus, tree_structure);

	} else {
		assert(0);
	}

	if (corpus)
		free_corpus(corpus);
	free(setting);

	pthread_exit(NULL);
}
