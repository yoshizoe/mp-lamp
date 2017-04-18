/*
 * mpi_tag.h
 *
 *  Created on: Apr 17, 2017
 *      Author: yuu
 */

#ifndef MP_SRC_MPI_TAG_H_
#define MP_SRC_MPI_TAG_H_

// use these as mpi tag
struct Tag {
	enum TaskType {
		// control tasks
		CONTROL_TASK_BEGIN = 0,
		DTD_REQUEST,
		DTD_REPLY,

		DTD_ACCUM_REQUEST, // request reporting accum count
		DTD_ACCUM_REPLY, // reduce closed set count

		BCAST_FINISH,
		CONTROL_TASK_END,

		// basic tasks
		BASIC_TASK_BEGIN,

		LAMBDA, // for lamp

		REQUEST, // request tasks
		REJECT, // reject requests
		GIVE, // send tasks

		BASIC_TASK_END,

		// third phase tasks
		THIRD_PHASE_BEGIN,

		RESULT_REQUEST,
		RESULT_REPLY,

		THIRD_PHASE_END,
	};
};

#endif /* MP_SRC_MPI_TAG_H_ */
