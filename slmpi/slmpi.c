#include "slmpi.h"

// This model includes mpi-Routines, which are used by the 
// function mpi_fit_pars(), which is available in the Remeis
// "isisscripts" at 
// http://www.sternwarte.uni-erlangen.de/git.public/?p=isisscripts

// Written by Thomas Dauser & Fritz Schwarm, 03/06/2016
// Few minors changed by Ashkbiz Danehkar, 08/03/2017

int g_master	= 0;
int g_rank	= 0;
int g_numtasks	= 0;
MPI_Comm g_mpicomm	= MPI_COMM_WORLD;
MPI_Status g_status;

RCL_MPI_MSG_LIST g_msg_list;


/// GENERAL MPI ROUTINES ///

/// Initialize MPI and return number of processes available
int rcl_mpi_init(void) {
	int mpi_active = 0;
	MPI_Initialized(&mpi_active);
	if(!mpi_active) {
		MPI_Init(NULL, NULL);
		MPI_Comm_dup(MPI_COMM_WORLD, &g_mpicomm);
	} else  MPI_Comm_dup(MPI_COMM_WORLD, &g_mpicomm);
	MPI_Comm_rank(g_mpicomm, &g_rank);
	MPI_Comm_size(g_mpicomm, &g_numtasks);
	return g_rank;
}
/// Probe for messages to be received, update the message list and return the number of new messages received
int rcl_mpi_update(void) {
	int receiving = 1, received = 0;
	while(receiving > 0) {
		receiving = rcl_mpi_msg_recv_double(&g_msg_list, MPI_ANY_TAG);
		if(receiving > 0) received++;
	}
	rcl_mpi_msg_update(&(g_msg_list.head));
	return received;
}
/// Shutdown MPI
void rcl_mpi_finalize(void) {
	MPI_Finalize();
	return;
}
/// Print MPI error
void rcl_mpi_printerror(MPI_Status status) {
	printf("MPI ERROR: %d\n", status.MPI_ERROR);
	return;
}
/// Get the message tag stored during the last rcl_mpi_... call
int rcl_mpi_status_tag() { return g_status.MPI_TAG; }
/// Get the message source stored during the last rcl_mpi_... call
int rcl_mpi_status_src() { return g_status.MPI_SOURCE; }
/// Return rank of process calling the function
int rcl_mpi_rank() { return g_rank; }
/// Return rank of process being the master process
int rcl_mpi_rank_master() { return g_master; }
/// Return number of processes within current MPI communicator
int rcl_mpi_numtasks() { return g_numtasks; }
/// Return 1 if the calling process is the master process or 0 otherwise
int rcl_mpi_master() {
	if(g_rank == g_master) return 1;
	else return 0;
}
/// Return 1 if the calling process is the master process or 0 otherwise
int rcl_mpi_client() {
	if(g_rank == g_master) return 0;
	else return 1;
}


/// MPI_FIT_PARS()
/// Recv <count> integers
int rcl_mpi_recv_state(int *buffer, int count) {
	return rcl_mpi_recv_int(buffer, count);
}

/// Send new best fit parameter combination to master process
int rcl_mpi_return_state(int *buffer, int count, int tag) {
	return rcl_mpi_send_int(buffer, count, g_master, tag);
}
/// Bcast state from master to clients (collective blocking communication)
int rcl_mpi_bcast_state(int *state) {
	return rcl_mpi_bcast_int(state, 1);
}
/// Bcast state from master to clients (collective blocking communication)
int rcl_mpi_dist_state(int *buffer, int count, int tag) {
	int i = 0;
	for(i = 0; i < g_numtasks; i++)
		if(i != g_rank) if(MPI_Send(buffer, count, MPI_INT, i, tag, g_mpicomm) != MPI_SUCCESS)
					rcl_mpi_printerror(g_status);
	return count;
}

/// Recv <count> doubles
int rcl_mpi_recv_param(double *buffer, int count) {
	return rcl_mpi_recv_double(buffer, count);
}
/// Send new best fit parameter combination to master process
int rcl_mpi_return_param(double *param, int count, int tag) {
	return rcl_mpi_send_double(param, count, g_master, tag);
}
/// Broadcast <N> parameters (doubles) to ALL ranks BUT the one calling the function
int rcl_mpi_bcast_param(double *param, int count) {
	if(MPI_Bcast(param, count, MPI_DOUBLE, g_master, g_mpicomm) != MPI_SUCCESS) {
		printf("%d: Bcast failed...\n", g_rank);
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Broadcast <N> parameters (doubles) from the current process to ALL ranks BUT the one calling the function
// Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_param_current(double *param, int count) {
	if(MPI_Bcast(param, count, MPI_DOUBLE, g_rank, g_mpicomm) != MPI_SUCCESS) {
		printf("%d: Bcast failed...\n", g_rank);
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Bcast state from master to clients (collective blocking communication)
int rcl_mpi_dist_param(double *buffer, int count, int tag) {
	int i = 0;
	for(i = 0; i < g_numtasks; i++)
		if(i != g_rank) if(MPI_Send(buffer, count, MPI_DOUBLE, i, tag, g_mpicomm) != MPI_SUCCESS)
					rcl_mpi_printerror(g_status);
	return count;
}


/// Collect confidence levels at master
int rcl_mpi_gather_conf(double *conf_min, double *conf_max, double *buf_min, double *buf_max, int N_par) {
	if(rcl_mpi_gather_double(conf_min, 1, buf_min, 1, g_master) != MPI_SUCCESS) return -1;
	if(rcl_mpi_gather_double(conf_max, 1, buf_max, 1, g_master) != MPI_SUCCESS) return -1;
/*
	if(rcl_mpi_master()) {
		int i = 0;
		printf("MASTER: HAVE %d PARAMETERS\n", N_par);
		for(i=0; i<N_par; i++)  {
			printf("buf_min[%d] = %lf\n", i, buf_min[i]);
			printf("buf_max[%d] = %lf\n", i, buf_max[i]);
		}
	}
*/
	return 0;
}



/// BLOCKING MPI ROUTINES ///

/// Wait until there is a message in the message queue and return the number of bytes of the message once there is one
int rcl_mpi_probe() {
	int msglen = 0;
	if(MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, g_mpicomm, &g_status) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
	if(MPI_Get_count(&g_status, MPI_BYTE, &msglen) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
	return msglen;
}

/// DOUBLES ///
/// Send <count> doubles to rank <dest> with tag <tag>
int rcl_mpi_send_double(double *buffer, int count, int dest, int tag) {
	return MPI_Send(buffer, count, MPI_DOUBLE, dest, tag, g_mpicomm);
}
/// Recv <count> doubles
int rcl_mpi_recv_double(double *buffer, int count) {
	return MPI_Recv(buffer, count, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, g_mpicomm, &g_status);
}
/// Broadcast <N> doubles to ALL ranks BUT the one calling the function
int rcl_mpi_bcast_double(double *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_DOUBLE, g_master, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Broadcast <N> chars to ALL ranks BUT the one calling the function
// Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_char(char *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_CHAR, g_master, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Broadcast <N> doubles from the current process to ALL ranks BUT the one calling the function
// Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_double_current(double *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_DOUBLE, g_rank, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Broadcast <N> chars from the current process to ALL ranks BUT the one calling the function
// Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_char_current(char *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_CHAR, g_rank, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Distribute <count> doubles over all processes
int rcl_mpi_scatter_double(double *input, int count, double *output, int root) {
	return MPI_Scatter(input, count, MPI_DOUBLE, output, count, MPI_DOUBLE, root, g_mpicomm);
}
/// Gathers together double values from all processes
int rcl_mpi_gather_double(double *input, int in_count, double *output, int out_count, int root) {
	return MPI_Gather(input, in_count, MPI_DOUBLE, output, out_count, MPI_DOUBLE, root, g_mpicomm);
}

/// INTEGERS ///
/// Send <count> integers to rank <dest> with tag <tag>
int rcl_mpi_send_int(int *buffer, int count, int dest, int tag) {
	return MPI_Send(buffer, count, MPI_INT, dest, tag, g_mpicomm);
}
/// Recv <count> integers
int rcl_mpi_recv_int(int *buffer, int count) {
	return MPI_Recv(buffer, count, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, g_mpicomm, &g_status);
}
/// Broadcast <N> doubles to ALL ranks BUT the one calling the function
int rcl_mpi_bcast_int(int *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_INT, g_master, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Broadcast <N> doubles from the current process to ALL ranks BUT the one calling the function
// Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_int_current(int *buffer, int count) {
	if(MPI_Bcast(buffer, count, MPI_INT, g_rank, g_mpicomm) != MPI_SUCCESS) {
		printf("Bcast failed...\n");
		rcl_mpi_printerror(g_status);
	}
	return count;
}
/// Distribute <count> integers over all processes
int rcl_mpi_scatter_int(int *input, int count, int *output, int root) {
	return MPI_Scatter(input, count, MPI_INT, output, count, MPI_INT, root, g_mpicomm);
}
/// Gathers together integer values from all processes
int rcl_mpi_gather_int(int *input, int count, int *output, int root) {
	return MPI_Gather(input, count, MPI_INT, output, count, MPI_INT, root, g_mpicomm);
}

/// Wait until all processes have reached this point
int rcl_mpi_barrier() {
	return MPI_Barrier(g_mpicomm);
}


/// NON-BLOCKING MPI ROUTINES ///

/// Test wether there is a message in the message queue or not and return the number of bytes (or 0)
int rcl_mpi_iprobe() {
	int msglen = 0, flag = 0;
	if(MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, g_mpicomm, &flag, &g_status) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
	if(flag) {
		if(MPI_Get_count(&g_status, MPI_BYTE, &msglen) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
		return msglen;
	} else return 0;
}

/// Test wether there is a message with tag <tag> in the message queue or not and return the number of bytes (or 0)
int rcl_mpi_iprobe_tag(int tag) {
	int msglen = 0, flag = 0;
	if(MPI_Iprobe(MPI_ANY_SOURCE, tag, g_mpicomm, &flag, &g_status) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
	if(flag) {
		if(MPI_Get_count(&g_status, MPI_BYTE, &msglen) != MPI_SUCCESS) rcl_mpi_printerror(g_status);
		return msglen;
	} else return 0;
}

/// Receive a maximum number of <count> doubles and return the number of doubles actually received
int rcl_mpi_irecv_double(double *buffer, int count) {
	int msglen = rcl_mpi_iprobe() / sizeof(double);
	if(msglen > 0) {
		if(MPI_Recv(buffer, msglen, MPI_DOUBLE, g_status.MPI_SOURCE, g_status.MPI_TAG, g_mpicomm, &g_status) != MPI_SUCCESS)
			rcl_mpi_printerror(g_status);
		return msglen;
	} else return 0;
}



/// FUTURE MESSAGE LIST STUFF FOR NON-BLOCKING COMMUNICATION API
/// Copy a message from one message structure to another one
void rcl_mpi_msg_cp(RCL_MPI_MSG *msg_out, RCL_MPI_MSG *msg_in)
{
	msg_out->index = msg_in->index;
	msg_out->tag = msg_in->tag;
	msg_out->src = msg_in->src;
	msg_out->dest = msg_in->dest;
	msg_out->flag = msg_in->flag;
	msg_out->type = msg_in->type;
	msg_out->count = msg_in->count;
	msg_out->bytes = msg_in->bytes;
	msg_out->request = msg_in->request;
	if(msg_out->data != NULL) free(msg_out->data);
	msg_out->data = (void*) calloc(msg_in->bytes, 1);
	memcpy(msg_out->data, msg_in->data, msg_in->bytes);
	msg_out->next = msg_in->next;
}
/// Free message data
void rcl_mpi_msg_free(RCL_MPI_MSG *msg) {
	if(msg->data != NULL) free(msg->data);
}
/// Add a message to the beginning of the message list
void rcl_mpi_msg_add(RCL_MPI_MSG_LIST *msg_list, RCL_MPI_MSG *msg)
{
	if(msg_list->N == 0) {
		msg_list->head = (RCL_MPI_MSG*) calloc(1, sizeof(RCL_MPI_MSG));
		rcl_mpi_msg_cp(msg_list->head, msg);
		msg_list->head->next = NULL;
	} else {
		RCL_MPI_MSG *old_head = msg_list->head;
		msg_list->head = (RCL_MPI_MSG*) calloc(1, sizeof(RCL_MPI_MSG));
		rcl_mpi_msg_cp(msg_list->head, msg);
		msg_list->head->next = old_head;
	}
	msg_list->N++;
}

/// Delete the very first message from a message list
void rcl_mpi_msg_del(RCL_MPI_MSG_LIST *msg_list)
{
	RCL_MPI_MSG *new_head = msg_list->head->next;
	if(msg_list->N > 1) new_head = msg_list->head->next;
	else new_head = NULL;
	rcl_mpi_msg_free(msg_list->head);
	free(msg_list->head);
	msg_list->head = new_head;
	msg_list->N--;
}
/// Non-blocking send of <count> doubles to rank <dest> with tag <tag>
int rcl_mpi_isend_double(double *buffer, int count, int dest, int tag)
{
	RCL_MPI_MSG msg;

	msg.index = 0;
	msg.tag = tag;
	msg.src = g_rank;
	msg.dest = dest;
	msg.flag = RCL_MPI_SEND;
	msg.type = RCL_MPI_DOUBLE;
	msg.count = count;
	msg.bytes = count * sizeof(double);

	memcpy(msg.data, buffer, msg.bytes);
	rcl_mpi_msg_add(&g_msg_list, &msg);
	if( MPI_Isend(msg.data, msg.count, MPI_DOUBLE, msg.dest, msg.tag, g_mpicomm, &msg.request) == MPI_SUCCESS ) {
		msg.flag = RCL_MPI_SENT;
		return msg.count;
	} else return -1;
}
int rcl_mpi_msg_recv_double(RCL_MPI_MSG_LIST *msg_list, int tag)
{
	int msglen = 1;
	RCL_MPI_MSG msg;
	if(tag == 0) tag = MPI_ANY_TAG;
	MPI_Status status;

	MPI_Probe(MPI_ANY_SOURCE, tag, g_mpicomm, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &msglen);
	if(msglen <= 0) return 0;
	else msg.data = (double*) calloc(msglen, sizeof(double));
	if( MPI_Recv(msg.data, msglen, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, g_mpicomm, &status) == MPI_SUCCESS )
	{
		msg.flag = RCL_MPI_RECV;
		msg.type = RCL_MPI_DOUBLE;
		msg.src = status.MPI_SOURCE;
		msg.tag = status.MPI_TAG;
		msg.dest = g_rank;
		msg.count = msglen;
		msg.bytes = msglen * sizeof(double);
		rcl_mpi_msg_add(msg_list, &msg);
	}
	return msglen;
}
/// Check for finished messages and mark the corresponding messages as completed or send them if they havent been sent yet
int rcl_mpi_msg_update(RCL_MPI_MSG **msg)
{
	int i_msg = 0, flag = 0;
	MPI_Status status;
	RCL_MPI_MSG *start_msg = (*msg), *prev_msg = NULL;

	// Test all messages in list for completion:
	while((*msg) != NULL)
	{
		i_msg++;
		if( (*msg)->flag == RCL_MPI_SEND )
		{
			// Send messages which arent sent yet:
			if( MPI_Isend((*msg)->data, (*msg)->bytes, MPI_BYTE, (*msg)->dest, (*msg)->tag, g_mpicomm, &((*msg)->request)) == MPI_SUCCESS )
				(*msg)->flag = RCL_MPI_SENT;
		} else if( (*msg)->flag == RCL_MPI_SENT )
		{
			MPI_Test(&(*msg)->request, &flag, &status);
			if(flag > 0) {
				// Message finished being sent and can now be deleted:
				(*msg)->flag = RCL_MPI_COMP;
				rcl_mpi_msg_free(*msg);
			}
		}
		prev_msg = (*msg);
		(*msg) = (*msg)->next;
	}
	(*msg) = start_msg;
	return i_msg;
}
/// Send <count> doubles in <buffer> to all other ranks and use the message tag <tag> (non-blocking)
int rcl_mpi_ibcast_double(void *buffer, int count, int tag) {
	int N = 0, i_rank = 0;
	for(i_rank = 0; i_rank < g_numtasks; i_rank++)  {
		if(i_rank != g_rank) { N += rcl_mpi_isend_double(buffer, count, i_rank, tag); }
	}
	return N;
}



/// To Provide "low" level mpi routines everything will be called ...

int rcl_mpi_org_recv_int (int *buffer, int count, int source, int tag) {
  return MPI_Recv(buffer, count, MPI_INT, source, tag, g_mpicomm, &g_status);
}

int rcl_mpi_org_recv_double (double *buffer, int count, int source, int tag) {
  return MPI_Recv(buffer, count, MPI_DOUBLE, source, tag, g_mpicomm, &g_status);
}

int rcl_mpi_org_send_int (int *buffer, int count, int dest, int tag) {
  return MPI_Send(buffer, count, MPI_INT, dest, tag, g_mpicomm);
}

int rcl_mpi_org_send_double (double *buffer, int count, int dest, int tag) {
  return MPI_Send(buffer, count, MPI_DOUBLE, dest, tag, g_mpicomm);
}

MPI_Request *g_request;
MPI_Status *cg_g_status;

void rcl_init_mpi_request(int n) {
  
  g_request = malloc( n * sizeof(MPI_Request));
  cg_g_status = malloc(n * sizeof(MPI_Status));
  
}

int rcl_mpi_org_isend_double (double *buffer , int count, int dest, int tag) {
  int index_g_request = dest;
  
  if(g_rank == 0) {
    index_g_request = dest - 1;
  }
  
  if(MPI_Isend(buffer, count, MPI_DOUBLE, dest, tag, g_mpicomm, &g_request[index_g_request]) != MPI_SUCCESS) {
    rcl_mpi_printerror(g_status);
  }
  return 0;
}

int rcl_mpi_org_isend_int (int *buffer , int count, int dest, int tag) {
  int index_g_request = dest;

  if(g_rank == 0) {
    index_g_request = dest - 1;
  }

  if(MPI_Isend(buffer, count, MPI_INT, dest, tag, g_mpicomm, &g_request[index_g_request]) != MPI_SUCCESS) {
    rcl_mpi_printerror(g_status);
  }
  return 0;
}

int rcl_mpi_org_irecv_double (double *buffer, int count, int source, int tag) {
  int index_g_request = source;
  
  if(g_rank == 0) {
    index_g_request = source - 1;
  }
  if(MPI_Irecv(buffer, count, MPI_DOUBLE, source, tag, g_mpicomm, &g_request[index_g_request]) != MPI_SUCCESS) {
    rcl_mpi_printerror(g_status);
  }

  return 0;
}

int rcl_mpi_org_irecv_int (int *buffer, int count, int source, int tag) {
  int index_g_request = source;

  if(g_rank == 0) {
    index_g_request = source - 1;
  }
  if(MPI_Irecv(buffer, count, MPI_INT, source, tag, g_mpicomm, &g_request[index_g_request]) != MPI_SUCCESS) {
    rcl_mpi_printerror(g_status);
  }

  return 0;
}

int rcl_mpi_org_wait(int rank) {

  return MPI_Wait(&g_request[rank], &g_status);
}

int rcl_mpi_org_waitany(int n, int *index) {
  if(MPI_Waitany(n, g_request, index, cg_g_status) != MPI_SUCCESS) {
    printf("PROBLEMS");
  }
//  printf("WAITINDEX %d\n", *index);
  if(g_rank == 0) {
  *index += 1;
  }
  return *index;
}

int rcl_mpi_org_waitsome(int n, int *outcount, int *nodes) {
  
  if(MPI_Waitsome(n, g_request, outcount, nodes, cg_g_status) != MPI_SUCCESS) {
     printf("PROBLEMS");
  }

  return *outcount;
}
