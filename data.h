#ifndef DATA_H
#define DATA_H

#include "common.h"
#include "datastructure.h"

t_corpus*	read_data(char* data_filename);
void		free_corpus(t_corpus* c);

std::vector<t_cat>	read_tree_structure(char* data_filename);
void			read_nodeid_to_docids(std::vector<t_cat> &tree_structure, char* data_filename, t_corpus *c);

#endif
