#ifndef _RCL_MPI_H
#define _RCL_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

// This model includes mpi-Routines, which are used by the 
// function mpi_fit_pars(), which is available in the Remeis
// "isisscripts" at 
// http://www.sternwarte.uni-erlangen.de/git.public/?p=isisscripts

// Written by Thomas Dauser & Fritz Schwarm, 03/06/2016
// Few minors changed by Ashkbiz Danehkar, 08/03/2017

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


enum RCL_MPI_DATATYPE { RCL_MPI_DOUBLE, RCL_MPI_INT };
enum RCL_MPI_FLAG { RCL_MPI_NAN, RCL_MPI_SEND, RCL_MPI_SENT, RCL_MPI_COMP, RCL_MPI_RECV };
enum RCL_MPI_TAG { RCL_MPI_BESTFIT };


typedef struct RCL_MPI_MSG
{
	int			index;
	int			tag;
	int			src;
	int			dest;
	int			flag;
	int			type;
	int			count;
	int			bytes;

	MPI_Request		request;

	void			*data;
	struct	RCL_MPI_MSG	*next;
} RCL_MPI_MSG;

typedef struct RCL_MPI_MSG_LIST
{
	int			N;
	struct RCL_MPI_MSG	*head;
} RCL_MPI_MSG_LIST;


/// GENERAL
int rcl_mpi_init(void);
void rcl_mpi_finalize(void);
void rcl_mpi_printerror(MPI_Status status);
int rcl_mpi_master(void);
int rcl_mpi_client(void);
int rcl_mpi_rank_master(void);
int rcl_mpi_rank(void);
int rcl_mpi_numtasks(void);
int rcl_mpi_status_tag();
int rcl_mpi_status_src();

/// MPI FIT PARS:
int rcl_mpi_recv_state(int *buffer, int count);
int rcl_mpi_return_state(int *buffer, int count, int tag);
int rcl_mpi_dist_state(int *buffer, int count, int tag);
int rcl_mpi_recv_param(double *buffer, int count);
int rcl_mpi_return_param(double *param, int count, int tag);
int rcl_mpi_dist_param(double *buffer, int count, int tag);
int rcl_mpi_bcast_param(double *param, int count);
int rcl_mpi_bcast_param_current(double *param, int count); // Add by A. Danehkar, 08/03/2017
//int rcl_mpi_gather_conf(double *conf_min, double *conf_max);
int rcl_mpi_gather_conf(double *conf_min, double *conf_max, double *buf_min, double *buf_max, int N_par);
int rcl_mpi_iprobe_tag(int tag);

/// BLOCKING
int rcl_mpi_barrier();
int rcl_mpi_probe();
/// DOUBLES
int rcl_mpi_send_double(double *buffer, int count, int dest, int tag);
int rcl_mpi_bcast_double(double *buffer, int count);
int rcl_mpi_bcast_char(char *buffer, int count); // Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_double_current(double *buffer, int count); // Add by A. Danehkar, 08/03/2017
int rcl_mpi_bcast_char_current(char *buffer, int count); // Add by A. Danehkar, 08/03/2017
int rcl_mpi_scatter_double(double *input, int count, double *output, int root);
int rcl_mpi_recv_double(double *buffer, int count);
int rcl_mpi_gather_double(double *input, int in_count, double *output, int out_count, int root);
/// INTEGERS
int rcl_mpi_send_int(int *buffer, int count, int dest, int tag);
int rcl_mpi_bcast_int(int *buffer, int count);
int rcl_mpi_bcast_int_current(int *buffer, int count); // Add by A. Danehkar, 08/03/2017
int rcl_mpi_scatter_int(int *input, int count, int *output, int root);
int rcl_mpi_recv_int(int *buffer, int countg);
int rcl_mpi_gather_int(int *input, int count, int *output, int root);

/// NON-BLOCKING
int rcl_mpi_iprobe();
int rcl_mpi_isend_double(double *buffer, int count, int dest, int tag);
int rcl_mpi_msg_recv_double(RCL_MPI_MSG_LIST *msg_list, int tag);
int rcl_mpi_msg_update(RCL_MPI_MSG **msg);
int rcl_mpi_ibcast_double(void *buffer, int count, int tag);
int rcl_mpi_irecv_double(double *buffer, int count);

/// MESSAGE LIST STUFF
int rcl_mpi_update(void);
void rcl_mpi_msg_cp(RCL_MPI_MSG *msg_out, RCL_MPI_MSG *msg_in);
void rcl_mpi_msg_free(RCL_MPI_MSG *msg);
void rcl_mpi_msg_add(RCL_MPI_MSG_LIST *msg_list, RCL_MPI_MSG *msg);
void rcl_mpi_msg_del(RCL_MPI_MSG_LIST *msg_list);

/// CG functions
int rcl_mpi_org_recv_int(int *buffer, int count, int source, int tag);
int rcl_mpi_org_recv_double(double *buffer, int count, int source, int tag);
int rcl_mpi_org_send_int (int *buffer, int count, int dest, int tag);
int rcl_mpi_org_send_double (double *buffer, int count, int dest, int tag);
void rcl_init_mpi_request(int n);
int rcl_mpi_org_waitsome(int n, int *index, int *nodes);
int rcl_mpi_org_irecv_double (double *buffer, int count, int source, int tag);
int rcl_mpi_org_isend_double (double *buffer , int count, int dest, int tag);
int rcl_mpi_org_irecv_int (int *buffer , int count, int source, int tag);
int rcl_mpi_org_isend_int (int *buffer , int count, int source, int tag);
int rcl_mpi_org_wait(int rank);
int rcl_mpi_org_waitany(int n, int *index) ;

#ifdef __cplusplus
}
#endif

#endif
