/*
 * datastructure.h
 *
 *  Created on: Jan 27, 2012
 *      Author: dok027
 */

#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include "common.h"

enum Mode {
	EST = 1,
	INF = 2
};

typedef struct
{
	Mode	mode;
	double	initial_gamma;
	double	initial_alpha;
	double	initial_eta;
	bool	estimate_alpha;

	int		num_topics;

	char	corpus_path[MAX_BUF];
	char	tree_structure_path[MAX_BUF];
	char	node_to_docids_path[MAX_BUF];
	char	output_path[MAX_BUF];

	int		num_threads;

	int		em_max_iter;
	int		cat_max_iter;
	int		kappa_tau_max_iter;
	int		doc_max_iter;

	double	em_converged;
    double	cat_converged;
	double	kappa_tau_converged;
    double	doc_converged;

    int		model_save_freq;
    int		random_seed;

    bool	corpus_init;
    int     num_docs_for_init;

    bool	warm_start;
	char	warm_start_path[MAX_BUF];

    // Options for evaluation
	char	model_path[MAX_BUF];
} t_setting;

typedef struct
{
    int* words;
    int* counts;
    int length;
    int total;
    int parent_index;
} t_document;

typedef struct
{
    t_document* docs;
    int 	num_terms;
    int 	num_docs;
    int		max_length;
} t_corpus;

typedef struct
{
	std::vector<int> 	docids;
	std::vector<int>	catids;
	int				parent_index;
} t_cat;

typedef struct
{
    double		gamma;
    double		eta;
	double*		alpha;
    int 		num_topics;
    int 		num_terms;
	int 		num_cats;
} t_tilda_model;

typedef struct
{
    int 		num_topics;
	int			num_terms;
	int 		num_cats;
	int			num_docs;
	double**	kappa;
	double*		tau;
	double**	nu;
	double***	rho;
	double**	lambda;
} t_tilda_var_model;

typedef struct
{
	double**	class_word;
} t_tilda_suffstats;

#endif /* DATASTRUCTURE_H_ */
