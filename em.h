#ifndef EM_H
#define EM_H

#include "common.h"
#include "datastructure.h"
#include "task_queue.h"
#include "model.h"
#include "doc_e_step.h"
#include "cat_e_step.h"
#include "utils.h"

void run_em(t_setting* setting, t_corpus* corpus, const std::vector<t_cat> tree_structure);
void *my_thread_function(void *lpParam);
double opt_lambda(t_tilda_suffstats* ss, t_tilda_var_model* model);

#endif
