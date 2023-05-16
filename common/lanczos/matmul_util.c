/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#include "lanczos.h"

void accum_xor(uint64 *dest, uint64 *src, uint32 n) {

	uint32 i;

	for (i = 0; i < (n & ~7); i += 8) {
		dest[i + 0] ^= src[i + 0];
		dest[i + 1] ^= src[i + 1];
		dest[i + 2] ^= src[i + 2];
		dest[i + 3] ^= src[i + 3];
		dest[i + 4] ^= src[i + 4];
		dest[i + 5] ^= src[i + 5];
		dest[i + 6] ^= src[i + 6];
		dest[i + 7] ^= src[i + 7];
	}
	for (; i < n; i++)
		dest[i] ^= src[i];
}


#ifdef HAVE_MPI

	/* The following is a functional replacement for
	   MPI_Allreduce(), but for large problems can be
	   configured to switch to a bucket accumulation
	   method that is asymptotically faster when the
	   vector is large and needs to be globally accumulated
	   and redistributed across a large number of nodes.

	   The algorithm uses the bucket strategy from the
	   paper "Global Combine on Mesh Architectures with 
	   Wormhole Routing". The implementation below is
	   based on code kindly contributed by Ilya Popovyan */

#if 1
#define GLOBAL_BREAKOVER 5000
#else
#define GLOBAL_BREAKOVER (uint32)(-1) /* turn off the fancy method */
#endif


/*------------------------------------------------------------------*/
static void global_xor_async(uint64 *send_buf, uint64 *recv_buf, 
			uint32 total_size, uint32 num_nodes, 
			uint32 my_id, MPI_Comm comm) {
	
	uint32 i;
	uint32 m, size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
	uint64 *curr_buf;
		
	/* split data */

	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	/* we expect a circular topology here */

	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
			
	/* stage 1
	   P_m sends P_{m+1} the m-th chunk of data while receiving 
	   another chunk from P_{m-1}, and does the summation op 
	   on the received chunk and (m-1)-th own chunk */

	m = my_id;
	size = chunk;
	if (my_id == num_nodes - 1)
		size += remainder;

	curr_buf = send_buf;

	for (i = 0; i < num_nodes - 1; i++) {
				
		/* asynchronously send the current chunk */

		MPI_TRY(MPI_Isend(curr_buf + m * chunk, size, 
				MPI_LONG_LONG, next_id, 97, 
				comm, &mpi_req))

		/* switch to the recvbuf after the first send */

		curr_buf = recv_buf;
				
		size = chunk;
		if ((int32)(--m) < 0) {
			m += num_nodes;
			size += remainder;
		}

		/* don't wait for send to finish, start the recv 
		   from the previous node */

		MPI_TRY(MPI_Recv(curr_buf + m * chunk, size,
				MPI_LONG_LONG, prev_id, 97, 
				comm, &mpi_status))

		/* combine the new chunk with our own */

		accum_xor(curr_buf + m * chunk,
			  send_buf + m * chunk, size);
		
		/* now wait for the send to end */

		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
	}	
		
	/* stage 2
	   P_m sends P_{m+1} m-th chunk of data, now containing 
	   a full summation of all m-th chunks in the comm,
	   while recieving another chunk from P_{m-1} and 
	   puts it to (m-1)-th own chunk */

	curr_buf = recv_buf + m * chunk;
	for (i = 0; i < num_nodes - 1; i++){
		
		/* async send to chunk the next proc in circle */

		MPI_TRY(MPI_Isend(curr_buf, size, MPI_LONG_LONG, 
				next_id, 98, comm, &mpi_req))
		
		size = chunk;
		curr_buf -= chunk;
		if (curr_buf < recv_buf) {
			curr_buf += chunk * num_nodes;			
			size += remainder;
		}		
		
		/* don't wait for send to finish, start the recv 
		   from the previous proc in circle, put the new 
		   data just where it should be in recv_buf */

		MPI_TRY(MPI_Recv(curr_buf, size, MPI_LONG_LONG,
				prev_id, 98, comm, &mpi_status))
				
		/* now wait for the send to end */

		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
	}
}

/*------------------------------------------------------------------*/
void global_xor(uint64 *send_buf, uint64 *recv_buf, 
		uint32 total_size, uint32 num_nodes, 
		uint32 my_id, MPI_Comm comm) {
	
	/* only get fancy for large buffers; even the
	   fancy method is only faster when many nodes 
	   are involved */

	if (total_size < GLOBAL_BREAKOVER || num_nodes < 2) {
		MPI_TRY(MPI_Allreduce(send_buf, recv_buf, total_size,
				MPI_LONG_LONG, MPI_BXOR, comm))
		return;
	}

	global_xor_async(send_buf, recv_buf, 
		total_size, num_nodes, my_id, comm);
}

/*------------------------------------------------------------------*/
void global_chunk_info(uint32 total_size, uint32 num_nodes, 
			uint32 my_id, uint32 *chunk_size, 
			uint32 *chunk_start) 
{
	uint32 chunk, remainder;
	
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	if (chunk_start)
		*chunk_start = my_id*chunk;
	
	if (chunk_size) {
		if (my_id == num_nodes - 1)
			*chunk_size = chunk + remainder;
		else
			*chunk_size = chunk;
	}
}

/*------------------------------------------------------------------*/
void global_xor_scatter(uint64 *send_buf, uint64 *recv_buf, 
			uint64 *scratch, uint32 total_size, 
			uint32 num_nodes, uint32 my_id, 
			MPI_Comm comm) {
	
	uint32 i;
	uint32 m, size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
    
	if (num_nodes == 1) {
		memcpy(recv_buf, send_buf, total_size * sizeof(uint64));
		return;
	}
    
	/* split data */
    
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	/* we expect a circular topology here */
    
	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
    
	/* P_m sends P_{m+1} the m-th chunk of data while receiving 
	   another chunk from P_{m-1}, and does the summation op 
	   on the received chunk and (m-1)-th own chunk */
    
	m = prev_id;
	size = chunk;
	if (m == num_nodes - 1)
		size += remainder;
        
	for (i = 0; i < num_nodes - 2; i++) {
        
		/* asynchroniously send the current chunk */
        
		MPI_TRY(MPI_Isend(send_buf + m * chunk, size, 
                          MPI_LONG_LONG, next_id, 95, 
                          comm, &mpi_req))
        
		/* switch to the recvbuf after the first send */
                
		size = chunk;
		if ((int32)(--m) < 0) {
			m += num_nodes;
			size += remainder;
		}
        
		/* don't wait for send to finish, start the recv 
		   from the previous node */
        
		MPI_TRY(MPI_Recv(scratch, size,
                         MPI_LONG_LONG, prev_id, 95, 
                         comm, &mpi_status))
        
		/* combine the new chunk with our own */
        
		accum_xor(send_buf + m * chunk, scratch, size);
		
		/* now wait for the send to end */
        
		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
	}	
    
	/* asynchronously send the current chunk */
    
	MPI_TRY(MPI_Isend(send_buf + m * chunk, size, 
			MPI_LONG_LONG, next_id, 95, 
			comm, &mpi_req))
    
	/* switch to the recvbuf after the first send */
    
	size = chunk;
	if ((int32)(--m) < 0) {
		m += num_nodes;
	        size += remainder;
	}
    
	/* don't wait for send to finish, start the recv 
	   from the previous node */
    
	MPI_TRY(MPI_Recv(recv_buf, size,
                     MPI_LONG_LONG, prev_id, 95, 
                     comm, &mpi_status))
    
	/* combine the new chunk with our own */
    
   	accum_xor(recv_buf, send_buf + m * chunk, size);
    
	/* now wait for the send to end */
    
	MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
}

/*------------------------------------------------------------------*/
void global_allgather(uint64 *send_buf, uint64 *recv_buf, 
                        uint32 total_size, uint32 num_nodes, 
                        uint32 my_id, MPI_Comm comm) {
	
	uint32 i;
	uint32 size, chunk, remainder;
	uint32 next_id, prev_id;
	MPI_Status mpi_status;
	MPI_Request mpi_req;
	uint64 *curr_buf;
    
	/* split data */
    
	chunk = total_size / num_nodes;
	remainder = total_size % num_nodes;	
	
	/* we expect a circular topology here */
    
	next_id = mp_modadd_1(my_id, 1, num_nodes);
	prev_id = mp_modsub_1(my_id, 1, num_nodes);
    
	/* P_m sends P_{m+1} m-th chunk of data, now containing 
	   a full summation of all m-th chunks in the comm,
	   while recieving another chunk from P_{m-1} and 
	   puts it to (m-1)-th own chunk */
    
	size = chunk;
	if (my_id == num_nodes - 1)
		size += remainder;

	curr_buf = recv_buf + my_id * chunk;
    
	/* put own part in place first */
	memcpy(curr_buf, send_buf, size * sizeof(uint64));
    
	for (i = 0; i < num_nodes - 1; i++){
		
		/* async send to chunk the next proc in circle */
        
		MPI_TRY(MPI_Isend(curr_buf, size, MPI_LONG_LONG, 
                          next_id, 96, comm, &mpi_req))
		
		size = chunk;
		curr_buf -= chunk;
		if (curr_buf < recv_buf) {
			curr_buf += chunk * num_nodes;			
			size += remainder;
		}		
		
		/* don't wait for send to finish, start the recv 
		   from the previous proc in circle, put the new 
		   data just where it should be in recv_buf */
        
		MPI_TRY(MPI_Recv(curr_buf, size, MPI_LONG_LONG,
                         prev_id, 96, comm, &mpi_status))
        
		/* now wait for the send to end */
        
		MPI_TRY(MPI_Wait(&mpi_req, &mpi_status))
	}
}

#endif /* HAVE_MPI */
