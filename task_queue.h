#ifndef TASK_QUEUE_H_
#define TASK_QUEUE_H_

#include <queue>
#include "common.h"

using namespace std;

enum t_task_mode {
	PRE = 1,
	POST = 2,
	DOC
};

struct t_task {
	int				node_index;
	t_task_mode		mode;
};

extern queue<t_task>	task_queue;
extern pthread_mutex_t 	mutex_for_queue;


#endif /* TASK_QUEUE_H_ */
