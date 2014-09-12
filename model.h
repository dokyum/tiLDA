#ifndef MODEL_H
#define MODEL_H

#include <memory.h>
#include "common.h"
#include "datastructure.h"
#include "cokus.h"

#define myrand() (double) (((unsigned long) randomMT()) / 4294967296.)

t_tilda_model* new_tilda_model(int num_terms, int num_topics, int num_cats);
void free_tilda_model(t_tilda_model*);
void save_tilda_model(t_tilda_model*, char*);
t_tilda_model* load_tilda_model(char*);

t_tilda_var_model* new_var_model(int num_cats, int num_docs, int num_topics, const t_corpus* corpus);
void free_var_model(t_tilda_var_model*);
void save_var_model(const t_tilda_var_model* model, char* model_root, const t_corpus* corpus);
t_tilda_var_model* load_var_model(char*, const t_corpus* corpus);

t_tilda_suffstats* new_tilda_suffstats(t_tilda_model* model);
void free_tilda_suffstats(t_tilda_suffstats* suffstats);
void corpus_initialize_ss(t_tilda_suffstats* ss, const t_tilda_model* model, const t_corpus* c, const int num_docs_for_init);
void collect_lambda_ss(t_tilda_suffstats* ss, const t_tilda_var_model* var_model, const t_corpus* c);

double** zero_init_double_matrix(const int height, const int width);
double* zero_init_double_array(const int size);
void free_double_matrix(double** p);

void warm_start_var_model(t_tilda_var_model* model, char* warm_start_path);

#endif
