#include "data.h"

t_corpus* read_data(char* data_filename) {
	FILE *fileptr;
	int count, word, n, l;
	t_corpus* c;
	printf("reading corpus data from %s\n", data_filename);
	c = (t_corpus*) malloc(sizeof(t_corpus));
	fileptr = fopen(data_filename, "r");
	fscanf(fileptr, "%10d", &(c->num_docs));
	fscanf(fileptr, "%10d", &(c->num_terms));
	c->max_length = 0;
	c->docs = (t_document*) malloc(sizeof(t_document) * (c->num_docs));
	for (n = 0; n < c->num_docs; n++) {
		fscanf(fileptr, "%10d", &(c->docs[n].length));
		if (c->docs[n].length > c->max_length) {
			c->max_length = c->docs[n].length;
		}
		c->docs[n].total = 0;
		c->docs[n].words = (int*) malloc(sizeof(int) * (c->docs[n].length));
		c->docs[n].counts = (int*) malloc(sizeof(int) * (c->docs[n].length));
		c->docs[n].parent_index = -1;
		for (l = 0; l < c->docs[n].length; l++) {
			fscanf(fileptr, "%10d:%10d", &word, &count);
			c->docs[n].counts[l] = count;
			c->docs[n].words[l] = word;
			c->docs[n].total += count;

			assert(word >= 0 && word < c->num_terms);
			assert(count >= 0);
		}
	}
	fclose(fileptr);
	printf("number of docs    : %d\n", c->num_docs);
	printf("number of terms   : %d\n", c->num_terms);
	return c;
}

std::vector<t_cat> read_tree_structure(char* data_filename) {
	FILE 			*fileptr;
	std::vector<t_cat>	tree_structure;
	int				num_nodes;
	int 			nodeid;
	int				childid;

	printf("reading tree structures from %s\n", data_filename);
	fileptr = fopen(data_filename, "r");
	fscanf(fileptr, "%10d", &num_nodes);
	for (int i = 0; i < num_nodes; ++i) {
		t_cat	onecat;
		onecat.catids.clear();
		onecat.docids.clear();
		onecat.parent_index = -1;
		tree_structure.push_back(onecat);
	}

	for (int i = 0; i < num_nodes; ++i) {
		fscanf(fileptr, "%10d\t%10d", &nodeid, &childid);
		assert(nodeid == i);
		while(childid != -1) {
			tree_structure[i].catids.push_back(childid);
			assert(tree_structure[childid].parent_index == -1);
			tree_structure[childid].parent_index = i;
			fscanf(fileptr, "%10d", &childid);
		}
	}
	fclose(fileptr);
	printf("number of nodes    : %d\n", num_nodes);
	assert(tree_structure.size() == num_nodes);

	return tree_structure;
}

void read_nodeid_to_docids(std::vector<t_cat> &tree_structure, char* data_filename, t_corpus *c)
{
	FILE 			*fileptr;
	int				num_nodes;
	int				docid;

	num_nodes = (int) tree_structure.size();

	printf("reading docids from %s\n", data_filename);
	fileptr = fopen(data_filename, "r");

	for (int i = 0; i < num_nodes; ++i) {
		fscanf(fileptr, "%10d", &docid);
		while(docid != -1) {
			assert(docid < c->num_docs);
			tree_structure[i].docids.push_back(docid);
			assert(-1 == c->docs[docid].parent_index);
			c->docs[docid].parent_index = i;
			fscanf(fileptr, "%10d", &docid);
		}
	}
	fclose(fileptr);

}

void free_corpus(t_corpus* c) {
	for (int i = 0; i < c->num_docs; ++i) {
		free(c->docs[i].words);
		free(c->docs[i].counts);
	}
	free(c->docs);
	free(c);
}

