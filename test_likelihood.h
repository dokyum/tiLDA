#ifndef TEST_LIKELIHOOD_H
#define TEST_LIKELIHOOD_H

#include "common.h"
#include "datastructure.h"
#include "model.h"
#include "doc_e_step.h"
#include "utils.h"

void test_likelihood(t_setting* setting, const t_corpus* corpus, const std::vector<t_cat> tree_structure);

#endif
