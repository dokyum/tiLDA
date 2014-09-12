#include "model.h"

t_tilda_model* new_tilda_model(int num_terms, int num_topics, int num_cats)
{
	printf("Creating new tilda model with num_terms %d\tnum_topics %d\tnum_cats %d\n", num_terms, num_topics, num_cats);
    t_tilda_model* model;

    model = (t_tilda_model*) malloc(sizeof(t_tilda_model));
    model->num_topics = num_topics;
    model->num_terms = num_terms;
    model->num_cats = num_cats;

    model->gamma = 0.0;
    model->eta = 0.0;
    model->alpha = zero_init_double_array(num_cats);

    return(model);
}

/*
 * deallocate hier corpus model
 *
 */
void free_tilda_model(t_tilda_model* model)
{
	free(model->alpha);
	free(model);
}

void save_tilda_model(t_tilda_model* model, char* model_root)
{
    char filename[MAX_BUF];
    FILE* fileptr;

    sprintf(filename, "%s.alpha", model_root);
    fileptr = fopen(filename, "w");
	for (int i = 0; i < model->num_cats; i++)
    {
		fprintf(fileptr, "%5.10f\n", model->alpha[i]);
    }
    fclose(fileptr);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
	fprintf(fileptr, "num_cats %d\n", model->num_cats);
    fprintf(fileptr, "gamma %5.10f\n", model->gamma);
    fprintf(fileptr, "eta %5.10f\n", model->eta);
    fclose(fileptr);
}

t_tilda_model* load_tilda_model(char* model_root)
{
    char    filename[MAX_BUF];
    FILE*   fileptr;
    int		num_terms, num_topics, num_cats;
    float	x;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
	fscanf(fileptr, "num_cats %d\n", &num_cats);

    t_tilda_model* model = new_tilda_model(num_terms, num_topics, num_cats);
    fscanf(fileptr, "gamma %f\n", &x);
    model->gamma = x;
    fscanf(fileptr, "eta %f\n", &x);
    model->eta = x;
    fclose(fileptr);

    sprintf(filename, "%s.alpha", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (int i = 0; i < num_cats; i++)
    {
		fscanf(fileptr, "%f", &x);
		model->alpha[i] = x;
    }
    fclose(fileptr);

    return(model);
}

t_tilda_var_model* new_var_model(int num_cats, int num_docs, int num_topics, const t_corpus* corpus)
{
	printf("Creating new var model with num_cats %d\tnum_docs %d\tnum_topics %d\n", num_cats, num_docs, num_topics);

	t_tilda_var_model* model;

    model = (t_tilda_var_model*) malloc(sizeof(t_tilda_var_model));
    model->num_cats = num_cats;
    model->num_docs = num_docs;
	model->num_topics = num_topics;
	model->num_terms = corpus->num_terms;

	model->kappa = zero_init_double_matrix(num_cats, num_topics);
	model->tau = zero_init_double_array(num_cats);
	model->nu = zero_init_double_matrix(num_docs, num_topics);
	model->rho = (double***) malloc(sizeof(double**) * num_docs);
	for (int i = 0; i < num_docs; ++i) {
		model->rho[i] = zero_init_double_matrix(corpus->docs[i].length, num_topics);
	}
	model->lambda = zero_init_double_matrix(num_topics, corpus->num_terms);
    return(model);
}

void free_var_model(t_tilda_var_model* model)
{
	for (int i = 0; i < model->num_docs; ++i) {
		free_double_matrix(model->rho[i]);
	}
	free(model->rho);
	free_double_matrix(model->kappa);
	free_double_matrix(model->nu);
	free(model->tau);
	free_double_matrix(model->lambda);
	free(model);
}

/*
 * save a variational model
 *
 */

void save_var_model(const t_tilda_var_model* model, char* model_root, const t_corpus* corpus)
{
    char filename[MAX_BUF];
    FILE* fileptr;

    sprintf(filename, "%s.kappa", model_root);
    fileptr = fopen(filename, "w");
    for (int i = 0; i < model->num_cats; i++)
    {
		for (int j = 0; j < model->num_topics; j++)
		{
			fprintf(fileptr, " %5.10f", model->kappa[i][j]);
		}
		fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.tau", model_root);
    fileptr = fopen(filename, "w");
	for (int i = 0; i < model->num_cats; i++)
    {
		fprintf(fileptr, " %5.10f\n", model->tau[i]);
    }
    fclose(fileptr);

//    sprintf(filename, "%s.nu", model_root);
//    fileptr = fopen(filename, "w");
//    for (int i = 0; i < model->num_docs; i++)
//    {
//		for (int j = 0; j < model->num_topics; j++)
//		{
//			fprintf(fileptr, " %5.10f", model->nu[i][j]);
//		}
//		fprintf(fileptr, "\n");
//    }
//    fclose(fileptr);
//
//    sprintf(filename, "%s.rho", model_root);
//    fileptr = fopen(filename, "w");
//    for (int i = 0; i < model->num_docs; i++)
//    {
//    	for (int l = 0; l < corpus->docs[i].length; l++) {
//			for (int j = 0; j < model->num_topics; j++)
//			{
//				fprintf(fileptr, " %5.10f", model->rho[i][l][j]);
//			}
//			fprintf(fileptr, "\n");
//    	}
//    }
//    fclose(fileptr);

    sprintf(filename, "%s.lambda", model_root);
    fileptr = fopen(filename, "w");
    for (int i = 0; i < model->num_topics; i++)
    {
		for (int j = 0; j < model->num_terms; j++)
		{
			fprintf(fileptr, " %5.10f", model->lambda[i][j]);
		}
		fprintf(fileptr, "\n");
    }
    fclose(fileptr);

    sprintf(filename, "%s.other", model_root);
    fileptr = fopen(filename, "w");
    fprintf(fileptr, "num_topics %d\n", model->num_topics);
    fprintf(fileptr, "num_terms %d\n", model->num_terms);
    fprintf(fileptr, "num_cats %d\n", model->num_cats);
    fprintf(fileptr, "num_docs %d\n", model->num_docs);
    fclose(fileptr);
}

t_tilda_var_model* load_var_model(char* model_root, const t_corpus* corpus)
{
    FILE* 	fileptr;
    char 	filename[MAX_BUF];
    int		num_cats, num_docs, num_topics, num_terms;
	double	x;

    sprintf(filename, "%s.other", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "num_cats %d\n", &num_cats);
    fscanf(fileptr, "num_docs %d\n", &num_docs);
    fclose(fileptr);

    assert(num_terms == corpus->num_terms);
    // ignore num_docs from the model
    num_docs = corpus->num_docs;

    t_tilda_var_model* model = new_var_model(num_cats, num_docs, num_topics, corpus);

    sprintf(filename, "%s.kappa", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (int i = 0; i < num_cats; i++)
    {
        for (int j = 0; j < num_topics; j++)
        {
            fscanf(fileptr, "%lf", &x);
            model->kappa[i][j] = x;
        }
    }
 	fclose(fileptr);

    sprintf(filename, "%s.tau", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (int i = 0; i < num_cats; i++)
    {
        fscanf(fileptr, "%lf", &x);
        model->tau[i] = x;
    }
 	fclose(fileptr);

//    sprintf(filename, "%s.nu", model_root);
//    printf("loading %s\n", filename);
//    fileptr = fopen(filename, "r");
//    for (int i = 0; i < num_docs; i++)
//    {
//        for (int j = 0; j < num_topics; j++)
//        {
//            fscanf(fileptr, "%f", &x);
//            model->nu[i][j] = x;
//        }
//    }
//    fclose(fileptr);

    sprintf(filename, "%s.lambda", model_root);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (int i = 0; i < num_topics; i++)
    {
        for (int j = 0; j < num_terms; j++)
        {
            fscanf(fileptr, "%lf", &x);
            model->lambda[i][j] = x;
        }
    }
    fclose(fileptr);

//    sprintf(filename, "%s.rho", model_root);
//    printf("loading %s\n", filename);
//    fileptr = fopen(filename, "r");
//    for (int i = 0; i < num_docs; i++)
//    {
//    	for (int l = 0; l < corpus->docs[i].length; l++) {
//			for (int j = 0; j < num_topics; j++)
//			{
//	            fscanf(fileptr, "%f", &x);
//	            model->rho[i][l][j] = x;
//			}
//    	}
//    }
//    fclose(fileptr);

    return(model);
}

t_tilda_suffstats* new_tilda_suffstats(t_tilda_model* model)
{
    int num_topics = model->num_topics;
    int num_terms = model->num_terms;

    t_tilda_suffstats* ss = (t_tilda_suffstats*) malloc(sizeof(t_tilda_suffstats));
	ss->class_word = zero_init_double_matrix(num_topics, num_terms);

    return(ss);
}

void free_tilda_suffstats(t_tilda_suffstats* suffstats)
{
	free_double_matrix(suffstats->class_word);
	free(suffstats);
}

void corpus_initialize_ss(t_tilda_suffstats* ss, const t_tilda_model* model, const t_corpus* c, const int num_docs_for_init)
{
    int num_topics = model->num_topics;
    int i, k, d, n;

    memset(ss->class_word[0], 0, sizeof(double) * (model->num_topics) * (model->num_terms));

    for (k = 0; k < num_topics; k++)
    {
        for (i = 0; i < num_docs_for_init; i++)
        {
        	do {
        		d = floor(myrand() * c->num_docs);
        	} while (c->docs[d].parent_index == -1); // only use documents that are in the hierarchy.

            printf("initialized with document %d\n", d);
            const t_document* doc = &(c->docs[d]);
            for (n = 0; n < doc->length; n++)
            {
                ss->class_word[k][doc->words[n]] += doc->counts[n];
            }
        }
        for (n = 0; n < model->num_terms; ++n) {
        	ss->class_word[k][n] += 1;
        }
    }
}

void collect_lambda_ss(t_tilda_suffstats* ss, const t_tilda_var_model* var_model, const t_corpus* c)
{
	memset(ss->class_word[0], 0, sizeof(double) * (var_model->num_topics) * (c->num_terms));

	for (int d = 0; d < c->num_docs; ++d) {
		const t_document* doc = &(c->docs[d]);
		if (-1 != doc->parent_index) { // only use documents that are in the hierarchy.
			for (int l = 0; l < doc->length; ++l) {
				for (int i = 0; i < var_model->num_topics; ++i) {
					ss->class_word[i][doc->words[l]] += doc->counts[l] * var_model->rho[d][l][i];
				}
			}
		}
	}
}

double** zero_init_double_matrix(const int height, const int width)
{
	double*		matrix;
	double**	index;
	double*		cursor;

	if (height == 0 || width == 0) {
		return NULL;
	}

	matrix = (double*) malloc(sizeof(double) * height * width);
	memset(matrix, 0, sizeof(double) * height * width);

	index = (double**) malloc(sizeof(double*) * height);
	cursor = matrix;
	for (int i = 0; i < height; ++i) {
		index[i] = cursor;
		cursor += width;
	}

	return index;
}

void free_double_matrix(double** p)
{
	if (p == NULL) return;
	free(p[0]);
	free(p);
}

double* zero_init_double_array(const int size)
{
	double*		index;

	index = (double*) malloc(sizeof(double) * size);
	memset(index, 0, sizeof(double) * size);

	return index;
}

void warm_start_var_model(t_tilda_var_model* model, char* warm_start_path)
{
    char 	filename[MAX_BUF];
    FILE* 	fileptr;
    int		num_cats, num_docs, num_topics, num_terms;
	float	x;

    sprintf(filename, "%s.other", warm_start_path);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "num_topics %d\n", &num_topics);
    fscanf(fileptr, "num_terms %d\n", &num_terms);
    fscanf(fileptr, "num_cats %d\n", &num_cats);
    fscanf(fileptr, "num_docs %d\n", &num_docs);
    fclose(fileptr);

    if (model->num_terms != num_terms || model->num_topics != num_topics) {
    	printf("Error: the specified warm start model does not match the setting.\n");
    	exit(-1);
    }

    sprintf(filename, "%s.lambda", warm_start_path);
    printf("loading %s\n", filename);
    fileptr = fopen(filename, "r");
    for (int i = 0; i < num_topics; i++)
    {
        for (int j = 0; j < num_terms; j++)
        {
            fscanf(fileptr, "%f", &x);
            model->lambda[i][j] = x;
        }
    }
    fclose(fileptr);

}
