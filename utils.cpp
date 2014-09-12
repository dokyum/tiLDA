#include "utils.h"

/*
 * given log(a) and log(b), return log(a + b)
 *
 */
double log_sum(double log_a, double log_b)
{
  double v;

  if (log_a < log_b)
  {
      v = log_b+log(1 + exp(log_a-log_b));
  }
  else
  {
      v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

 /**
   * Proc to calculate the value of the trigamma, the second
   * derivative of the loggamma function. Accepts positive matrices.
   * From Abromowitz and Stegun.  Uses formulas 6.4.11 and 6.4.12 with
   * recurrence formula 6.4.6.  Each requires workspace at least 5
   * times the size of X.
   *
   **/

double trigamma(double x)
{
    double p;
    int i;

    x=x+6;
    p=1/(x*x);
    p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
    for (i=0; i<6 ;i++)
    {
        x=x-1;
        p=1/(x*x)+p;
    }
    return(p);
}


/*
 * taylor approximation of first derivative of the log gamma function
 *
 */

double digamma(double x)
{
    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
	0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}

double tetragamma(double x)
{
	double p;
    int i;

    x=x+6;
    p=1/(x*x);
    p=(((((0.3 - 0.833333333333333 * p) * p - 0.166666666666666) * p + 0.166666666666666) * p - 0.5) * p - 1/x - 1) * p;
    for (i=0; i<6 ;i++)
    {
        x=x-1;
        p = p - 2 / (x*x*x);
    }
    return(p);
}

void setup_output_directory(t_setting* setting)
{
    char 	setting_path[MAX_BUF];
    FILE*	setting_file;
    time_t 	rawtime;
    struct 	tm* timeinfo;

    mkdir(setting->output_path, S_IRUSR|S_IWUSR|S_IXUSR);
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    sprintf(setting_path, "%s/%04d%02d%02d_%02d%02d%02d_inference_setting", setting->output_path, 1900 + timeinfo->tm_year, 1 + timeinfo->tm_mon, timeinfo->tm_mday,
    		timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
    setting_file = fopen(setting_path, "w");
    if (EST == setting->mode) {
    	fprintf(setting_file, "MODE EST\n");
    	fprintf(setting_file, "initial_gamma %lf\n", setting->initial_gamma);
    	fprintf(setting_file, "initial_alpha %lf\n", setting->initial_alpha);
    	fprintf(setting_file, "initial_eta %lf\n", setting->initial_eta);

    	if (setting->estimate_alpha) {
    		fprintf(setting_file, "estimate_alpha 1\n");
    	} else {
    		fprintf(setting_file, "estimate_alpha 0\n");
    	}

		fprintf(setting_file, "num_topics %d\n", setting->num_topics);

		fprintf(setting_file, "corpus_path %s\n", setting->corpus_path);
		fprintf(setting_file, "tree_structure_path %s\n", setting->tree_structure_path);
		fprintf(setting_file, "node_to_docids_path %s\n", setting->node_to_docids_path);
		fprintf(setting_file, "output_path %s\n", setting->output_path);

		fprintf(setting_file, "num_threads %d\n", setting->num_threads);

		fprintf(setting_file, "em_max_iter %d\n", setting->em_max_iter);
		fprintf(setting_file, "cat_max_iter %d\n", setting->cat_max_iter);
		fprintf(setting_file, "kappa_tau_max_iter %d\n", setting->kappa_tau_max_iter);
		fprintf(setting_file, "doc_max_iter %d\n", setting->doc_max_iter);

		fprintf(setting_file, "em_converged %lf\n", setting->em_converged);
		fprintf(setting_file, "cat_converged %lf\n",  setting->cat_converged);
		fprintf(setting_file, "kappa_tau_converged %lf\n",  setting->kappa_tau_converged);
		fprintf(setting_file, "doc_converged %lf\n", setting->doc_converged);

		fprintf(setting_file, "model_save_freq %d\n", setting->model_save_freq);
		fprintf(setting_file, "random_seed %d\n", setting->random_seed);

    	if (setting->corpus_init) {
    		fprintf(setting_file, "corpus_init 1\n");
    	} else {
    		fprintf(setting_file, "corpus_init 0\n");
    	}

    	fprintf(setting_file, "num_docs_for_init %d\n", setting->num_docs_for_init);

    	if (setting->warm_start) {
    		fprintf(setting_file, "warm_start 1\n");
    	} else {
    		fprintf(setting_file, "warm_start 0\n");
    	}

    	fprintf(setting_file, "warm_start_path %s\n", setting->warm_start_path);

    } else {
    	assert(0);
    }

    fclose(setting_file);
}
