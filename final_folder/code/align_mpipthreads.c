/*
 * Exact genetic sequence alignment
 * (Using brute force)
 *
 * MPI version
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2023/2024
 *
 * v1.3
 *
 * (c) 2024, Arturo Gonzalez-Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<sys/time.h>
#include<mpi.h>
#include<pthread.h>


/* Arbitrary value to indicate that no matches are found */
#define	NOT_FOUND	-1

/* Arbitrary value to restrict the checksums period */
#define CHECKSUM_MAX	65535


/* 
 * Utils: Function to get wall time
 */
double cp_Wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}

/*
 * Utils: Random generator
 */
#include "rng.c"


/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */


#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif
/*
 * Function: Increment the number of pattern matches on the sequence positions
 * 	This function can be changed and/or optimized by the students
 */
void increment_matches( int pat, unsigned long *pat_found, unsigned long *pat_length, int *seq_matches ) {
	unsigned long ind;	
	for( ind=0; ind<pat_length[pat]; ind++) {
		seq_matches[ pat_found[pat] + ind ]++;
	}
}

/*
 * Function: Fill random sequence or pattern
 */
void generate_rng_sequence( rng_t *random, float prob_G, float prob_C, float prob_A, char *seq, unsigned long length) {
	unsigned long ind; 
	for( ind=0; ind<length; ind++ ) {
		double prob = rng_next( random );
		if( prob < prob_G ) seq[ind] = 'G';
		else if( prob < prob_C ) seq[ind] = 'C';
		else if( prob < prob_A ) seq[ind] = 'A';
		else seq[ind] = 'T';
	}
}

/*
 * Function: Copy a sample of the sequence
 */
void copy_sample_sequence( rng_t *random, char *sequence, unsigned long seq_length, unsigned long pat_samp_loc_mean, unsigned long pat_samp_loc_dev, char *pattern, unsigned long length) {
	/* Choose location */
	unsigned long  location = (unsigned long)rng_next_normal( random, (double)pat_samp_loc_mean, (double)pat_samp_loc_dev );
	if ( location > seq_length - length ) location = seq_length - length;
	if ( location <= 0 ) location = 0;

	/* Copy sample */
	unsigned long ind; 
	for( ind=0; ind<length; ind++ )
		pattern[ind] = sequence[ind+location];
}

/**
 * Defining global variables, structures and functions for thread operations
 */

typedef struct pattern_recognition_args{
	unsigned long seq_start;		// IN: thread - specific, starting point of the search
	unsigned long seq_chunk_size;	// IN: thread - specific, dimension of the sequence chunk to be analyzed
	int pat_th_start;				// IN: thread - specific, starting pattern		
	int *pat_lock;					// IN: global
	int pat_n_start;				// IN: global
	int pat_n_end;					// IN: global
	unsigned long seq_length;		// IN: global
	unsigned long *pat_length;		// IN: global
	char *sequence;					// IN: global
	char **pattern;					// IN: global
	unsigned long *pat_found; 		// OUT: thread - specific
	
}pattern_recognition_args;

void *pattern_recognition(void *args){
	// Taking the function arguments
	pattern_recognition_args *params = (pattern_recognition_args*) args;
	unsigned long seq_start = params->seq_start;
	unsigned long seq_chunk_size = params->seq_chunk_size;
	int pat_th_start = params->pat_th_start;
	int pat_n_start = params->pat_n_start;				
	int pat_n_end = params->pat_n_end;
	unsigned long seq_length = params->seq_length;
	unsigned long *pat_length = params->pat_length;
	char *sequence = params->sequence;
	char **pattern = params->pattern;
	unsigned long *pat_found = params->pat_found;
	int *pat_lock = params->pat_lock;

	// Pattern search
	unsigned long lind;
	int pat;
	for(int idx=0; idx < pat_n_end - pat_n_start; idx++) {
		// Each thread starts by looking up a different thread,
		// so they don't interleave in looking for the same pattern all together
		pat = ((idx + pat_th_start) % (pat_n_end - pat_n_start)) + pat_n_start;
		
		// If the counter associated with the pattern is bigger than 0, it means
		// the pattern has already been found and it is useless to seek it
		if (pat_lock[pat] > 0){
			continue;
		}

		// For each posible starting position
		for(unsigned long start=seq_start; start < seq_start + seq_chunk_size; start++) {
			// Ensuring that the search can actually be brought out inside the sequence
			if (start > seq_length - pat_length[pat]){
				break;
			}

			// For each pattern element
			for(lind=0; lind<pat_length[pat]; lind++) {
				// Stop this test when different nucleotids are found
				if ( sequence[start + lind] != pattern[pat][lind] ) break;
			}
			// Check if the loop ended with a match
			if (lind == pat_length[pat]) {
				pat_lock[pat]++;
				pat_found[pat] = start;
				break;
			}
		}
	}
	
	return NULL;
}

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

/*
 * Function: Allocate new patttern
 */
char *pattern_allocate( rng_t *random, unsigned long pat_rng_length_mean, unsigned long pat_rng_length_dev, unsigned long seq_length, unsigned long *new_length ) {

	/* Random length */
	unsigned long length = (unsigned long)rng_next_normal( random, (double)pat_rng_length_mean, (double)pat_rng_length_dev );
	if ( length > seq_length ) length = seq_length;
	if ( length <= 0 ) length = 1;

	/* Allocate pattern */
	char *pattern = (char *)malloc( sizeof(char) * length );
	if ( pattern == NULL ) {
		fprintf(stderr,"\n-- Error allocating a pattern of size: %lu\n", length );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* Return results */
	*new_length = length;
	return pattern;
}


/*
 * Function: Regenerate a sample of the sequence
 */
void generate_sample_sequence( rng_t *random, rng_t random_seq, float prob_G, float prob_C, float prob_A, unsigned long seq_length, unsigned long pat_samp_loc_mean, unsigned long pat_samp_loc_dev, char *pattern, unsigned long length ) {
	/* Choose location */
	unsigned long  location = (unsigned long)rng_next_normal( random, (double)pat_samp_loc_mean, (double)pat_samp_loc_dev );
	if ( location > seq_length - length ) location = seq_length - length;
	if ( location <= 0 ) location = 0;

	/* Regenerate sample */
	rng_t local_random = random_seq;
	rng_skip( &local_random, location );
	generate_rng_sequence( &local_random, prob_G, prob_C, prob_A, pattern, length);
}


/*
 * Function: printf usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<seq_length> <prob_G> <prob_C> <prob_A> <pat_rng_num> <pat_rng_length_mean> <pat_rng_length_dev> <pat_samples_num> <pat_samp_length_mean> <pat_samp_length_dev> <pat_samp_loc_mean> <pat_samp_loc_dev> <pat_samp_mix:B[efore]|A[fter]|M[ixed]> <long_seed>\n");
	fprintf(stderr,"\n");
}



/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	/* 0. Default output and error without buffering, forces to write immediately */
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	/* 1. Read scenary arguments */
	/* 1.0. Init MPI before processing arguments */
	MPI_Init( &argc, &argv );
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* 1.1. Check minimum number of arguments */
	if (argc < 15) {
		fprintf(stderr, "\n-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 1.2. Read argument values */
	unsigned long seq_length = atol( argv[1] );
	float prob_G = atof( argv[2] );
	float prob_C = atof( argv[3] );
	float prob_A = atof( argv[4] );
	if ( prob_G + prob_C + prob_A > 1 ) {
		fprintf(stderr, "\n-- Error: The sum of G,C,A,T nucleotid probabilities cannot be higher than 1\n\n");
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	prob_C += prob_G;
	prob_A += prob_C;

	int pat_rng_num = atoi( argv[5] );
	unsigned long pat_rng_length_mean = atol( argv[6] );
	unsigned long pat_rng_length_dev = atol( argv[7] );
	
	int pat_samp_num = atoi( argv[8] );
	unsigned long pat_samp_length_mean = atol( argv[9] );
	unsigned long pat_samp_length_dev = atol( argv[10] );
	unsigned long pat_samp_loc_mean = atol( argv[11] );
	unsigned long pat_samp_loc_dev = atol( argv[12] );

	char pat_samp_mix = argv[13][0];
	if ( pat_samp_mix != 'B' && pat_samp_mix != 'A' && pat_samp_mix != 'M' ) {
		fprintf(stderr, "\n-- Error: Incorrect first character of pat_samp_mix: %c\n\n", pat_samp_mix);
		show_usage( argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	unsigned long seed = atol( argv[14] );

#ifdef DEBUG
	/* DEBUG: printf arguments */
	if ( rank == 0 ) {
		printf("\nArguments: seq_length=%lu\n", seq_length );
		printf("Arguments: Accumulated probabilitiy G=%f, C=%f, A=%f, T=1\n", prob_G, prob_C, prob_A );
		printf("Arguments: Random patterns number=%d, length_mean=%lu, length_dev=%lu\n", pat_rng_num, pat_rng_length_mean, pat_rng_length_dev );
		printf("Arguments: Sample patterns number=%d, length_mean=%lu, length_dev=%lu, loc_mean=%lu, loc_dev=%lu\n", pat_samp_num, pat_samp_length_mean, pat_samp_length_dev, pat_samp_loc_mean, pat_samp_loc_dev );
		printf("Arguments: Type of mix: %c, Random seed: %lu\n", pat_samp_mix, seed );
		printf("\n");
	}
#endif // DEBUG

	/* 2. Initialize data structures */
	/* 2.1. Skip allocate and fill sequence */
	rng_t random = rng_new( seed );
	rng_skip( &random, seq_length );

	/* 2.2. Allocate and fill patterns */
	/* 2.2.1 Allocate main structures */
	int pat_number = pat_rng_num + pat_samp_num;
	unsigned long *pat_length = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	char **pattern = (char **)malloc( sizeof(char*) * pat_number );
	if ( pattern == NULL || pat_length == NULL ) {
		fprintf(stderr,"\n-- Error allocating the basic patterns structures for size: %d\n", pat_number );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 2.2.2 Allocate and initialize ancillary structure for pattern types */
	int ind;
	unsigned long lind;
	#define PAT_TYPE_NONE	0
	#define PAT_TYPE_RNG	1
	#define PAT_TYPE_SAMP	2
	char *pat_type = (char *)malloc( sizeof(char) * pat_number );
	if ( pat_type == NULL ) {
		fprintf(stderr,"\n-- Error allocating ancillary structure for pattern of size: %d\n", pat_number );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_NONE;

	/* 2.2.3 Fill up pattern types using the chosen mode */
	switch( pat_samp_mix ) {
	case 'A':
		for( ind=0; ind<pat_rng_num; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		break;
	case 'B':
		for( ind=0; ind<pat_samp_num; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		break;
	default:
		if ( pat_rng_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		}
		else if ( pat_samp_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		}
		else if ( pat_rng_num < pat_samp_num ) {
			int interval = pat_number / pat_rng_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_RNG;
				else pat_type[ind] = PAT_TYPE_SAMP;
		}
		else {
			int interval = pat_number / pat_samp_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_SAMP;
				else pat_type[ind] = PAT_TYPE_RNG;
		}
	}

	/* 2.2.4 Generate the patterns */
	for( ind=0; ind<pat_number; ind++ ) {
		if ( pat_type[ind] == PAT_TYPE_RNG ) {
			pattern[ind] = pattern_allocate( &random, pat_rng_length_mean, pat_rng_length_dev, seq_length, &pat_length[ind] );
			generate_rng_sequence( &random, prob_G, prob_C, prob_A, pattern[ind], pat_length[ind] );
		}
		else if ( pat_type[ind] == PAT_TYPE_SAMP ) {
			pattern[ind] = pattern_allocate( &random, pat_samp_length_mean, pat_samp_length_dev, seq_length, &pat_length[ind] );
#define REGENERATE_SAMPLE_PATTERNS
#ifdef REGENERATE_SAMPLE_PATTERNS
			rng_t random_seq_orig = rng_new( seed );
			generate_sample_sequence( &random, random_seq_orig, prob_G, prob_C, prob_A, seq_length, pat_samp_loc_mean, pat_samp_loc_dev, pattern[ind], pat_length[ind] );
#else
			copy_sample_sequence( &random, sequence, seq_length, pat_samp_loc_mean, pat_samp_loc_dev, pattern[ind], pat_length[ind] );
#endif
		}
		else {
			fprintf(stderr,"\n-- Error internal: Paranoic check! A pattern without type at position %d\n", ind );
			MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
		}
	}
	free( pat_type );

	/* Avoid the usage of arguments to take strategic decisions
	 * In a real case the user only has the patterns and sequence data to analize
	 */
	argc = 0;
	argv = NULL;
	pat_rng_num = 0;
	pat_rng_length_mean = 0;
	pat_rng_length_dev = 0;
	pat_samp_num = 0;
	pat_samp_length_mean = 0;
	pat_samp_length_dev = 0;
	pat_samp_loc_mean = 0;
	pat_samp_loc_dev = 0;
	pat_samp_mix = '0';

	/* 2.3. Other result data and structures */
	int pat_matches = 0;

	/* 2.3.1. Other results related to patterns */
	unsigned long *pat_found;
	pat_found = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	if ( pat_found == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux pattern structure for size: %d\n", pat_number );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	
	/* 3. Start global timer */
	MPI_Barrier( MPI_COMM_WORLD );
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	/* 2.1. Allocate and fill sequence */
	char *sequence = (char *)malloc( sizeof(char) * seq_length );
	if ( sequence == NULL ) {
		fprintf(stderr,"\n-- Error allocating the sequence for size: %lu\n", seq_length );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}
	random = rng_new( seed );
	generate_rng_sequence( &random, prob_G, prob_C, prob_A, sequence, seq_length);

	/* 3. Other results related to the main sequence */
	int *seq_matches;
	// The "+1" will be used to to communicate pat_matches
	seq_matches = (int *)malloc(sizeof(int) * (seq_length + 1));
	if ( seq_matches == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux sequence structures for size: %lu\n", (seq_length + 1));
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	/* 4. Initialize ancillary structures */
	for( ind=0; ind<pat_number; ind++) {
		pat_found[ind] = (unsigned long)NOT_FOUND;			// pat_found[i] means that pattern i has been found starting at position pat_found[i]
	}
	for( lind=0; lind<seq_length; lind++) {
		seq_matches[lind] = 0;						// seq_matches[i] means that in position i seq_matches[i] patterns have been found
	}

	/* 5. Search for each pattern */
	/* 5.1. Variable initialization */
	/* 5.1.1. Rank variables */
	int pat;
	int start_idx = (pat_number / size) * rank;
	int end_idx = (pat_number / size) * (rank + 1);
	if (rank == size - 1){
		end_idx = pat_number;
	}
#ifdef DEBUG
	/* DEBUG: printf sequence and patterns */
	int flag = 0;
	if (rank != 0){
		MPI_Recv(&flag, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	printf("-----------------\n");
	printf("IN RANK %d\n", rank);
	printf("Sequence: ");
	for( lind=0; lind<seq_length; lind++ ) 
		printf( "%c", sequence[lind] );
	printf("\n-----------------\n");
	printf("Patterns: %d\n", end_idx - start_idx);
	int debug_pat;
	for( debug_pat=start_idx; debug_pat<end_idx; debug_pat++ ) {
		printf( "Pat[%d]: ", debug_pat );
		for( lind=0; lind<pat_length[debug_pat]; lind++ ) 
			printf( "%c", pattern[debug_pat][lind] );
		printf("\n");
	}
	printf("-----------------\n\n");
	if (rank != size - 1){
		MPI_Send(&flag, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
	}
#endif

	/* 5.1.2. Thread management */
	int num_threads = NUM_THREADS;
	pthread_t** thread_handles = malloc(sizeof(pthread_t*) * num_threads);
	pattern_recognition_args* args = malloc(sizeof(pattern_recognition_args) * num_threads);

	/* 5.1.3. Input parameters */
	unsigned long* seq_start = (unsigned long*) malloc(sizeof(unsigned long) * num_threads);
	unsigned long* seq_chunk_size = (unsigned long*) malloc(sizeof(unsigned long) * num_threads);
	int* pat_th_start = (int*) malloc(sizeof(int) * num_threads);

	/* 5.1.4. Control constructs */
	// pat_lock[i] gets incremented when a thread discovers pattern i.
	// Whenever a thread is beginning to check for a thread, it checks that it has not
	// already been found through this list.
	// Note that this list does not guarantee through sync constructs that no two threads
	// will look for the same pattern, but it does not matter as long as the pattern is found.
	// This is more time-efficient than using synchronization constructs, since they are not 
	// stricly necessary.
	int *pat_lock = (int*) calloc(sizeof(int), pat_number);

	/* 5.1.5. Return parameters */
	unsigned long **thread_pat_found = (unsigned long**) malloc(sizeof(unsigned long*)*num_threads);
	
	/* 5.2. Spawning threads */
	for (int i = 0; i < num_threads; i++){
		/* 5.2.1 Initializing thread handles */
		thread_handles[i] = malloc(sizeof(pthread_t));

		/* 5.2.2 Initializing function parameters */
		seq_start[i] = ((unsigned long)i) * seq_length/((unsigned long)num_threads);
		pat_th_start[i] = i * (end_idx - start_idx)/num_threads;
		if (i == num_threads - 1){
			seq_chunk_size[i] = seq_length - seq_start[i];
		} else {
			seq_chunk_size[i] = seq_length/((unsigned long)num_threads);
		}

		thread_pat_found[i] = (unsigned long*) malloc(sizeof(unsigned long)*pat_number);
		for(ind = 0; ind < pat_number; ind++) {
			thread_pat_found[i][ind] = (unsigned long)NOT_FOUND;
		}

		args[i] = (pattern_recognition_args) {
			seq_start[i],
			seq_chunk_size[i],
			pat_th_start[i],
			pat_lock,
			start_idx,
			end_idx,
			seq_length,
			pat_length,
			sequence,
			pattern,
			thread_pat_found[i]
			};

		/* 5.2.3 Creating threads */
		pthread_create(thread_handles[i], NULL, pattern_recognition, (void*) &args[i]);
	}

	/* 5.3. Joining threads */
	// Threads are joined in reverse order since the first threads
	// could take more, given they might have to match longer sequences
	for (int i = num_threads - 1; i >= 0; i--){
		pthread_join(*thread_handles[i], NULL);

		/* 5.4 Gathering partial results and assemblating them */
		/* 5.4.1. Pattern matches */
		for (pat = start_idx; pat < end_idx; pat++){
			if (pat_found[pat] == (unsigned long)NOT_FOUND && thread_pat_found[i][pat] != (unsigned long)NOT_FOUND){
				pat_found[pat] = thread_pat_found[i][pat];
				pat_matches++;

				/* 5.4.2. Sequence matches */
				increment_matches(pat, thread_pat_found[i], pat_length, seq_matches);
			}
		}
	}

	/* 6. Exchanging information */
	/* 6.0.1. This will be needed later for exchanging seq_matches: splitting the sequence in chunks of size INT_MAX */
	unsigned long index;
	int quantity;
	// Calculating the number of necessary MPI_Sends (with the assumptions seq_length is not
	// bigger than INT_MAX^2)
	int rounds = (seq_length + 1)/INT_MAX;
	if ((seq_length + 1)%INT_MAX != 0){
		rounds++;
	}
	
	/* 6.0.2. Creating requests: one for pat_found and one for each chunk of seq_matches */
	MPI_Request comm_req;

	if (rank == 0){
		/* 6.1. Receiving partial information at rank 0 */
		/* 6.1.1. Pattern position receving (pat_found)*/
		// Creating the information about each process:
		// Elements per thread and displacement in the final array
		int *recvcounts = (int*) malloc(sizeof(int) * size);
		int *displs = (int*) malloc(sizeof(int) * size);
		for (int i = 0; i < size; i++){
			start_idx = (pat_number / size) * i;
			end_idx = (pat_number / size) * (i + 1);
			if (i == size - 1){
				end_idx = pat_number;
			}
			recvcounts[i] = end_idx - start_idx;
			displs[i] = start_idx;
		}

		MPI_Igatherv(MPI_IN_PLACE, recvcounts[0], MPI_UNSIGNED_LONG, pat_found, recvcounts, displs, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD, &comm_req);

		/* 6.1.2. Matches in the sequence (seq_matches) and pat_matches */
		seq_matches[seq_length] = pat_matches;
		for (int round = 0; round < rounds; round++){
				index = round*INT_MAX;
				if (round == rounds - 1){
					quantity = seq_length + 1 - index;
				} else {
					quantity = INT_MAX;
				}
				MPI_Reduce(MPI_IN_PLACE, &seq_matches[index], quantity, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			}

		/* 6.1.3. Waiting for results */
		MPI_Wait(&comm_req, MPI_STATUS_IGNORE);

		/* 6.1.4. Reporting pat_matches on the original variable */
		pat_matches = seq_matches[seq_length];

		/* 6.1.5. Freeing up temporary structures */
		free(recvcounts);
		free(displs);
	
	} else {
		/* 6.2. Sending data from non-master processes */
		/* 6.2.1. Gathering pat_found */
		MPI_Igatherv(&pat_found[start_idx], end_idx-start_idx, MPI_UNSIGNED_LONG, NULL, NULL, NULL, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD, &comm_req);
		
		/* 6.2.2. Sending pat_matches through seq_matches */
		seq_matches[seq_length] = pat_matches;
		
		/* 6.2.3.  MPI_Ireduce supports sending INT_MAX elements per call. Splitting up calls */
		for (int round = 0; round < rounds; round++){
			index = round*INT_MAX;
			if (round == rounds - 1){
				quantity = seq_length + 1 - index;
			} else {
				quantity = INT_MAX;
			}
			MPI_Reduce(&seq_matches[index], NULL, quantity, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		/* 6.2.4. Waiting MPI calls */
		MPI_Wait(&comm_req, MPI_STATUS_IGNORE);
		}
	}

	/* 7. Check sums */
	unsigned long checksum_matches = 0;
	unsigned long checksum_found = 0;
	for( ind=0; ind < pat_number; ind++) {
		if ( pat_found[ind] != (unsigned long)NOT_FOUND )
			checksum_found = ( checksum_found + pat_found[ind] ) % CHECKSUM_MAX;
	}
	for( lind=0; lind < seq_length; lind++) {
		checksum_matches = ( checksum_matches + seq_matches[lind] ) % CHECKSUM_MAX;
	}

#ifdef DEBUG
	/* DEBUG: Write results */
	if (rank == 0){
		printf("-----------------\n");
		printf("Found start:");
		for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
			printf( " %lu", pat_found[debug_pat] );
		}
		printf("\n");
		printf("-----------------\n");
		printf("Matches:");
		for( lind=0; lind<seq_length; lind++ ) 
			printf( " %d", seq_matches[lind] );
		printf("\n");
		printf("-----------------\n");
	}
#endif // DEBUG

	/* Free local resources */	
	free( sequence );
	free( seq_matches );
	for (int j = 0; j < num_threads; j++){
		free(thread_handles[j]);
		free(thread_pat_found[j]);
	}	
	free(thread_pat_found);
	free(thread_handles);
	free(args);
	free(seq_start);
	free(seq_chunk_size);
	free(pat_th_start);
	free(pat_lock);

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

	/* 8. Stop global time */
	MPI_Barrier( MPI_COMM_WORLD );
	ttotal = cp_Wtime() - ttotal;

	/* 9. Output for leaderboard */
	if ( rank == 0 ) {
		printf("\n");
		/* 9.1. Total computation time */
		printf("Time: %lf\n", ttotal );

		/* 9.2. Results: Statistics */
		printf("Result: %d, %lu, %lu\n\n", 
				pat_matches,
				checksum_found,
				checksum_matches );
	}
				
	/* 10. Free resources */	
	int i;
	for( i=0; i<pat_number; i++ ) free( pattern[i] );
	free( pattern );
	free( pat_length );
	free( pat_found );

	/* 11. End */
	MPI_Finalize();
	return 0;
}
