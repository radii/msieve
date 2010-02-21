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

#include <msieve.h>
#include <signal.h>

msieve_obj *g_curr_factorization = NULL;

/*--------------------------------------------------------------------*/
void handle_signal(int sig) {

	msieve_obj *obj = g_curr_factorization;

	printf("\nreceived signal %d; shutting down\n", sig);
	
	if (obj && (obj->flags & MSIEVE_FLAG_SIEVING_IN_PROGRESS))
		obj->flags |= MSIEVE_FLAG_STOP_SIEVING;
	else
		_exit(0);
}

/*--------------------------------------------------------------------*/
void get_random_seeds(uint32 *seed1, uint32 *seed2) {

	uint32 tmp_seed1, tmp_seed2;

	/* In a multithreaded program, every msieve object
	   should have two unique, non-correlated seeds
	   chosen for it */

#if !defined(WIN32) && !defined(_WIN64)

	FILE *rand_device = fopen("/dev/urandom", "r");

	if (rand_device != NULL) {

		/* Yay! Cryptographic-quality nondeterministic randomness! */

		fread(&tmp_seed1, sizeof(uint32), (size_t)1, rand_device);
		fread(&tmp_seed2, sizeof(uint32), (size_t)1, rand_device);
		fclose(rand_device);
	}
	else

#endif
	{
		/* <Shrug> For everyone else, sample the current time,
		   the high-res timer (hopefully not correlated to the
		   current time), and the process ID. Multithreaded
		   applications should fold in the thread ID too */

		uint64 high_res_time = read_clock();
		tmp_seed1 = ((uint32)(high_res_time >> 32) ^
			     (uint32)time(NULL)) * 
			    (uint32)getpid();
		tmp_seed2 = (uint32)high_res_time;
	}

	/* The final seeds are the result of a multiplicative
	   hash of the initial seeds */

	(*seed1) = tmp_seed1 * ((uint32)40499 * 65543);
	(*seed2) = tmp_seed2 * ((uint32)40499 * 65543);
}

/*--------------------------------------------------------------------*/
void print_usage(char *progname) {

	printf("\nMsieve v. %d.%02d\n", MSIEVE_MAJOR_VERSION, 
					MSIEVE_MINOR_VERSION);

	printf("\nusage: %s [options] [one_number]\n", progname);
	printf("\nnumbers starting with '0' are treated as octal,\n"
		"numbers starting with '0x' are treated as hexadecimal\n");
	printf("\noptions:\n"
	         "   -s <name> save intermediate results to <name>\n"
		 "             instead of the default %s\n"
	         "   -l <name> append log information to <name>\n"
		 "             instead of the default %s\n"
	         "   -i <name> read one or more integers to factor from\n"
		 "             <name> (default worktodo.ini) instead of\n"
		 "             from the command line\n"
		 "   -m        manual mode: enter numbers via standard input\n"
		 "   -mb <num> # of megabytes of memory for postprocessing \n"
		 "             (set automatically if unspecified or zero)\n"
	         "   -q        quiet: do not generate any log information,\n"
		 "             only print any factors found\n"
	         "   -d <min>  deadline: if still sieving after <min>\n"
		 "             minutes, shut down gracefully (default off)\n"
		 "   -r <num>  stop after finding <num> relations\n"
		 "   -p        run at idle priority\n"
	         "   -v        verbose: write log information to screen\n"
		 "             as well as to logfile\n"
#ifdef HAVE_CUDA
		 "   -g <num>  use GPU <num>, 0 <= num < (CUDA GPUs)>\n"
#endif
	         "   -t <num>  use at most <num> threads\n\n"
		 " elliptic curve options:\n"
		 "   -e        perform 'deep' ECM, seek factors > 15 digits\n\n"
		 " quadratic sieve options:\n"
		 "   -c        client: only perform sieving\n\n"
		 " number field sieve options:\n"
		 "   -n        use the number field sieve (80+ digits only;\n"
		 "             performs all NFS tasks in order)\n"
	         "   -nf <name> read from / write to NFS factor base file\n"
		 "             <name> instead of the default %s\n"
		 "   -np [X,Y] perform only NFS polynomial selection; if\n"
		 "             specified, only cover leading coefficients\n"
		 "             in the range from X to Y inclusive\n"
		 "   -ns [X,Y] perform only NFS sieving; if specified,\n"
		 "             handle sieve lines X to Y inclusive\n"
		 "   -nc       perform only NFS combining (all phases)\n"
		 "   -nc1 [X,Y] perform only NFS filtering. Filtering will\n"
		 "             track ideals >= X (determined automatically\n"
		 "             if 0 or unspecified) and will only use the \n"
		 "             first Y relations (or all relations, if 0 \n"
		 "             or unspecified)\n"
		 "   -nc2      perform only NFS linear algebra\n"
		 "   -ncr      perform only NFS linear algebra, restarting\n"
		 "             from a previous checkpoint\n"
		 "   -nc3 [X,Y] perform only NFS square root (compute \n"
		 "             dependency numbers X through Y, 1<=X,Y<=64)\n",
		 MSIEVE_DEFAULT_SAVEFILE, MSIEVE_DEFAULT_LOGFILE,
		 MSIEVE_DEFAULT_NFS_FBFILE);
}

/*--------------------------------------------------------------------*/
void factor_integer(char *buf, uint32 flags,
		    char *savefile_name,
		    char *logfile_name,
		    char *nfs_fbfile_name,
		    uint32 *seed1, uint32 *seed2,
		    uint32 max_relations,
		    uint64 nfs_lower,
		    uint64 nfs_upper,
		    enum cpu_type cpu,
		    uint32 cache_size1,
		    uint32 cache_size2,
		    uint32 num_threads,
		    uint32 mem_mb,
		    uint32 which_gpu) {
	
	char *int_start, *last;
	msieve_obj *obj;
	msieve_factor *factor;

	/* point to the start of the integer or expression;
	   if the start point indicates no integer is present,
	   don't try to factor it :) */

	last = strchr(buf, '\n');
	if (last)
		*last = 0;
	int_start = buf;
	while (*int_start && !isdigit(*int_start) &&
			*int_start != '(' ) {
		int_start++;
	}
	if (*int_start == 0)
		return;

	g_curr_factorization = msieve_obj_new(int_start, flags,
					savefile_name, logfile_name,
					nfs_fbfile_name,
					*seed1, *seed2, max_relations,
					nfs_lower, nfs_upper, cpu,
					cache_size1, cache_size2,
					num_threads, mem_mb, which_gpu);
	if (g_curr_factorization == NULL) {
		printf("factoring initialization failed\n");
		return;
	}

	msieve_run(g_curr_factorization);

	if (!(g_curr_factorization->flags & MSIEVE_FLAG_FACTORIZATION_DONE)) {
		printf("\ncurrent factorization was interrupted\n");
		exit(0);
	}

	/* If no logging is specified, at least print out the
	   factors that were found */

	if (!(g_curr_factorization->flags & (MSIEVE_FLAG_USE_LOGFILE |
					MSIEVE_FLAG_LOG_TO_STDOUT))) {
		factor = g_curr_factorization->factors;

		printf("\n");
		printf("%s\n", buf);
		while (factor != NULL) {
			char *factor_type;

			if (factor->factor_type == MSIEVE_PRIME)
				factor_type = "p";
			else if (factor->factor_type == MSIEVE_COMPOSITE)
				factor_type = "c";
			else
				factor_type = "prp";

			printf("%s%d: %s\n", factor_type, 
					(int32)strlen(factor->number), 
					factor->number);
			factor = factor->next;
		}
		printf("\n");
	}

	/* save the current value of the random seeds, so that
	   the next factorization will pick up the pseudorandom
	   sequence where this factorization left off */

	*seed1 = g_curr_factorization->seed1;
	*seed2 = g_curr_factorization->seed2;

	/* free the current factorization struct. The following
	   avoids a race condition in the signal handler */

	obj = g_curr_factorization;
	g_curr_factorization = NULL;
	if (obj)
		msieve_obj_free(obj);
}

#ifdef WIN32
DWORD WINAPI countdown_thread(LPVOID pminutes) {
	DWORD minutes = *(DWORD *)pminutes;

	if (minutes > 0x7fffffff / 60000)
		minutes = 0;            /* infinite */

	Sleep(minutes * 60000);
	raise(SIGINT);
	return 0;
}

#else
void *countdown_thread(void *pminutes) {
	uint32 minutes = *(uint32 *)pminutes;

	if (minutes > 0xffffffff / 60)
		minutes = 0xffffffff / 60;   /* infinite */

	sleep(minutes * 60);
	raise(SIGINT);
	return NULL;
}
#endif

/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {

	char buf[400];
	uint32 seed1, seed2;
	char *savefile_name = NULL;
	char *logfile_name = NULL;
	char *infile_name = "worktodo.ini";
	char *nfs_fbfile_name = NULL;
	uint32 flags;
	char manual_mode = 0;
	int i;
	int32 deadline = 0;
	uint32 max_relations = 0;
	uint64 nfs_lower = 0;
	uint64 nfs_upper = 0;
	enum cpu_type cpu;
	uint32 cache_size1; 
	uint32 cache_size2; 
	uint32 num_threads = 0;
	uint32 mem_mb = 0;
	uint32 which_gpu = 0;
		
	get_cache_sizes(&cache_size1, &cache_size2);
	cpu = get_cpu_type();

	if (signal(SIGINT, handle_signal) == SIG_ERR) {
	        printf("could not install handler on SIGINT\n");
	        return -1;
	}
	if (signal(SIGTERM, handle_signal) == SIG_ERR) {
	        printf("could not install handler on SIGTERM\n");
	        return -1;
	}     

	flags = MSIEVE_FLAG_USE_LOGFILE;

	i = 1;
	buf[0] = 0;
	while (i < argc) {
		if (argv[i][0] == (char)('-')) {
			switch(tolower(argv[i][1])) {
			case 'h':
			case '?':
				print_usage(argv[0]);
				return 0;

			case 'i':
			case 's':
			case 'l':
				if (i + 1 < argc && 
				    argv[i+1][0] != '-') {
					if (tolower(argv[i][1]) == 'i')
						infile_name = argv[i+1];
					else if (tolower(argv[i][1]) == 's')
						savefile_name = argv[i+1];
					else
						logfile_name = argv[i+1];
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 'm':
				if (argv[i][2] == 'b' &&
				    i + 1 < argc && isdigit(argv[i+1][0])) {
					mem_mb = atol(argv[i+1]);
					i += 2;
				}
				else {
					manual_mode = 1;
					i++;
				}
				break;

			case 'e':
				flags |= MSIEVE_FLAG_DEEP_ECM;
				i++;
				break;

			case 'n':
				if (argv[i][2] == 'p') {
					flags |= MSIEVE_FLAG_NFS_POLY;
				}
				else if (argv[i][2] == 's') {
					flags |= MSIEVE_FLAG_NFS_SIEVE;
				}
				else if (argv[i][2] == 'c') {
					if (argv[i][3] == '1')
						flags |= MSIEVE_FLAG_NFS_FILTER;
					else if (argv[i][3] == '2')
						flags |= MSIEVE_FLAG_NFS_LA;
					else if (argv[i][3] == 'r')
						flags |= 
						     MSIEVE_FLAG_NFS_LA |
						     MSIEVE_FLAG_NFS_LA_RESTART;
					else if (argv[i][3] == '3')
						flags |= MSIEVE_FLAG_NFS_SQRT;
					else if (argv[i][3] == 0)
						flags |= 
						     MSIEVE_FLAG_NFS_FILTER |
						     MSIEVE_FLAG_NFS_LA |
						     MSIEVE_FLAG_NFS_SQRT;
				}
				else if (argv[i][2] == 0) {
					flags |= MSIEVE_FLAG_NFS_POLY |
						 MSIEVE_FLAG_NFS_SIEVE |
						 MSIEVE_FLAG_NFS_FILTER |
						 MSIEVE_FLAG_NFS_LA |
						 MSIEVE_FLAG_NFS_SQRT;
				}

				if (i + 1 < argc && argv[i+1][0] != '-') {
					if (argv[i][2] == 'f') {
						nfs_fbfile_name = argv[i+1];
						i++;
					}
					else if ((argv[i][2] == 's' ||
						  argv[i][2] == 'c' ||
						  argv[i][2] == 'p') &&
						  strchr(argv[i+1], ',') != 
						  		NULL ) {
						char *tmp;
						nfs_lower = (uint64)strtod(
							argv[i+1], &tmp);
						tmp++;
						nfs_upper = (uint64)strtod(
							tmp, NULL);
						i++;
					}
				}
				i++;
				break;
					
			case 'q':
				flags &= ~(MSIEVE_FLAG_USE_LOGFILE |
					   MSIEVE_FLAG_LOG_TO_STDOUT);
				i++;
				break;
					
			case 'd':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					deadline = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 'r':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					max_relations = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
					
			case 't':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					num_threads = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
#ifdef HAVE_CUDA					
			case 'g':
				if (i + 1 < argc && isdigit(argv[i+1][0])) {
					which_gpu = atol(argv[i+1]);
					i += 2;
				}
				else {
					print_usage(argv[0]);
					return -1;
				}
				break;
#endif					
			case 'c':
				flags |= MSIEVE_FLAG_SKIP_QS_CYCLES;
				i++;
				break;

			case 'v':
				flags |= MSIEVE_FLAG_LOG_TO_STDOUT;
				i++;
				break;

			case 'p':
				set_idle_priority();
				i++;
				break;

			default:
				print_usage(argv[0]);
				return -1;
			}
		}
		else {
			if (isdigit(argv[i][0]) || argv[i][0] == '(' )
				strncpy(buf, argv[i], sizeof(buf));
			i++;
		}
	}

	get_random_seeds(&seed1, &seed2);

	if (deadline) {
#if defined(WIN32) || defined(_WIN64)
		DWORD thread_id;
		CreateThread(NULL, 0, countdown_thread, 
				&deadline, 0, &thread_id);
#else
		pthread_t thread_id;
		pthread_create(&thread_id, NULL, 
				countdown_thread, &deadline);
#endif
	}

	if (isdigit(buf[0]) || buf[0] == '(' ) {
		factor_integer(buf, flags, savefile_name, 
				logfile_name, nfs_fbfile_name,
				&seed1, &seed2,
				max_relations, 
				nfs_lower, nfs_upper, cpu,
				cache_size1, cache_size2,
				num_threads, mem_mb, which_gpu);
	}
	else if (manual_mode) {
		while (1) {
			printf("\n\nnext number: ");
			fflush(stdout);
			buf[0] = 0;
			fgets(buf, (int)sizeof(buf), stdin);
			factor_integer(buf, flags, savefile_name, 
					logfile_name, nfs_fbfile_name,
					&seed1, &seed2,
					max_relations, 
					nfs_lower, nfs_upper, cpu,
					cache_size1, cache_size2,
					num_threads, mem_mb, which_gpu);
			if (feof(stdin))
				break;
		}
	}
	else {
		FILE *infile = fopen(infile_name, "r");
		if (infile == NULL) {
			printf("cannot open input file '%s'\n", infile_name);
			return 0;
		}

		while (1) {
			buf[0] = 0;
			fgets(buf, (int)sizeof(buf), infile);
			factor_integer(buf, flags, savefile_name, 
					logfile_name, nfs_fbfile_name,
					&seed1, &seed2,
					max_relations, 
					nfs_lower, nfs_upper, cpu,
					cache_size1, cache_size2,
					num_threads, mem_mb, which_gpu);
			if (feof(infile))
				break;
		}
		fclose(infile);
	}

	return 0;
}
