/*
 *  Needleman-Wunsch alignment for short reads
 *
 *  Qinhu Wang
 *  Northwest A&F University
 *  wangqinhu@nwafu.edu.cn
 *
 *  Jul 5, 2014
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#define MAXCOL     100                 // max column each line
#define MAXSEQ     100                 // max sequene length
#define MAXFNL     256                 // max filename length
#define MATCH      4                   // match score
#define MISMATCH  -3                   // mismatch score
#define GAP       -4                   // gap score

void usage() {
	
	printf("\nNeedleman-Wunsch Alignment\n\n");
	printf("usage:\n\n");
	printf("nwalign <options>\n\n");
	printf("  -h  help\n");
	printf("  -i  input sequence file 1\n");
	printf("  -j  input sequence file 2\n");
	printf("  -a  input sequence a, directly in command line\n");
	printf("  -b  input sequence b, directly in command line\n");
	printf("  -s  output the alignment score only\n");
    printf("  -r  output the raw alignment\n");
    printf("  -f  output the full alignment\n");
	exit(1);

}

/* define nucleotude */
int isnuc (char nuc) {

	switch ( nuc ) {

		case 'A': return 1;
		case 'C': return 1;
		case 'G': return 1;
		case 'T': return 1;
		case 'U': return 1;
		case 'a': return 1;
		case 'c': return 1;
		case 'g': return 1;
		case 't': return 1;
		case 'u': return 1;
		default:  return 0;

	}

}

/* convert seq to UPPER case */
void upper (char seq[]) {

	int		i;

	i = 0;

	while ( seq[i] != '\0' ) {
		seq[i] = toupper( seq[i]);
		i++;
	}

}

/* read a fasta file and pass the sequence to an array */
void readseq (char filename[], char seq[]) {

	FILE	*file;                           // hold fasta file
	int		i;                               // index for seq[]
	char	line[ MAXCOL ];                  // store each line
	int		j;                               // index for line[]

	file = fopen( filename, "r" );

	if ( file == NULL ) {

		printf("Error: file %s does not exist!\n", filename);
		exit(1);

	}

	i = 0;
	while ( fgets( line, sizeof line, file ) ) {

		/* if a fasta header line was found */
		if ( line[0] == '>' ) continue;
		
		/* now try to parse the fasta body */
		for ( j = 0; j < strlen(line); j++ ) {
			if ( isnuc( line[j] ) ) seq[i++] = line[j];
		}

	} // end of while
	seq[i] = '\0';

	fclose( file );

} // end readseq

/* scoring scheme */
int score (char a, char b) {
	
	int		score;                              // final score

	if ( a == '-' || b == '-' )
		score = GAP;
	else
		score = ( a == b ) ? MATCH : MISMATCH;
	
	return score;

} // end score

/* max3: return the max score */
int max3 (int d, int l, int u) {
	
	int		max;                             // max score

	max = ( d > l ) ? d : l;
	max = ( u > max ) ? u : max;

	return max;
	
} // end max3

/* mss: return the max score source */
char mss (int d, int l, int u) {

	int		max;                             // max score
	char	mss;                             // max score source
	
	max = max3( d, l, u );

	if (max == d) mss = '\\';
	if (max == l) mss = '|';
	if (max == u) mss = '_';

	return mss;

} // end of mss

/* revers a string: used for the last step of back tracing */
void revseq (char seq[]) {

	int		i;                               // index for seq and revseq
	int		l;                               // seq length

	l = strlen( seq );
	
	char	hold[ l ];                       // hold seq
	
	strcpy(hold, seq);

	for ( i = 0; i < l; i++ ) {
		seq[i] = hold[l-1-i]; 
	}

} // end of revseq

/* get the alignment string */
void align_str (char aln1[], char aln2[], char aln[]) {

	int		i, l;                            // index and length for alignment

	l = strlen( aln1 );

	if ( l != strlen( aln2 ) ) {
	
		printf("Error, seq1 and seq2 are not equal in length!\n");
		exit(1);
	
	}

	// get the alignment string
	for ( i = 0; i < l; i++ ) {

		if ( aln1[i] == aln2[i] ) {
			aln[i] = ':';
		} else {
			aln[i] = ' ';
		}
	
	}
	
	aln[i] = '\0';

} // end of align_str

/* needleman-wunsch alignment function */
int nwalign (char seq1[], char seq2[], char bts1[], char bts2[], char aln[]) {

	int		mat[ MAXSEQ + 1 ][ MAXSEQ + 1 ]; // score matrix
	int		i, j;                            // matrix indices
	int		m, n;                            // sequences length
	char	btm[ MAXSEQ + 1 ][ MAXSEQ + 1 ]; // back trace mark
	int		d, u, l;                         // source score
	int 	bi;                              // back trace index, forward

	/* initialization the matrix */
	m = strlen( seq1 );
	n = strlen( seq2 );

	for ( i = 0; i <= m; i++ ) {
		mat[i][0] = GAP * i;
		btm[i][0] = '|';
	}
	for ( j = 0; j <= n; j++ ) {
		mat[0][j] = GAP * j;
		btm[0][j] = '_';
	}

	/* filling the matrix */
	for ( i = 1; i <= m; i ++ )
		for ( j = 1; j <= n; j++ ) {

			// diagonal
			d = mat[i-1][j-1] + score( seq1[i-1], seq2[j-1] ); 
			// up
			u = mat[i-1][j] + score( seq1[i-1], '-' ); 
			// left
			l = mat[i][j-1] + score( '-', seq2[j-1] );
			// find the max score
			mat[i][j] = max3( d, u, l );
			// mark the trace path
			btm[i][j] = mss( d, u, l );

		}

	/* back trace  */

	/*
	// print the back trace matrix
	for ( i = 0; i <= m; i++ ) {
		for ( j = 0; j <= n; j++ ) printf("%2d %c\t", mat[i][j], btm[i][j]);
		printf("\n");
	}
	*/
	
	// back trace
	i = m;
	j = n;
	bi = 0;
	while ( ( i != 0 ) || ( j != 0) ) {
		if ( btm[i][j] == '\\' ) {
			bts1[bi] = seq1[i - 1];
			bts2[bi] = seq2[j - 1];
			i--;
			j--;
		} else if ( btm[i][j] == '_')  {
			bts1[bi] = '-';
			bts2[bi] = seq2[j - 1];
			j--;
		} else if ( btm[i][j] == '|' ) {
			bts1[bi] = seq1[i - 1];
			bts2[bi] = '-';
			i--;
		} else {
			break;
		}
		bi++;
	}
	bts1[bi] = '\0';
	bts2[bi] = '\0';
	revseq(bts1);
	revseq(bts2);

	// get aligment string
	align_str( bts1, bts2, aln );
	
	return mat[m][n];

} // end of nwalign

/* format alignment output */
void print_align (char aln1[], char aln2[], char aln[]) {

	int		i,l;                             // index and length for alignment
	int		j, k, m, n;                      // vars for format control

	l = strlen( aln1 );

	// output alignment header
	printf("1        10        20        30        40        50\n");
	printf("....:....|....:....|....:....|....:....|....:....|\n\n");

	// format output
	m = l;
	n = 0;
	for ( i = 0; i < l; i = i + 50 ) {
		// decide the length of the char to print
		if ( m > 50 ) {
			k = 50;
		} else {
			k = m % 50;
		}
		// print seqs/alignment
		// seq1
		for ( j = 0; j < k; j++ ) {
			printf("%c", aln1[i+j]);
		}
		printf("\n");
		// alingment
		for ( j = 0; j < k; j++ ) {
			printf("%c", aln[i+j]);
		}
		printf("\n");
		// seq2
		for ( j = 0; j < k; j++ ) {
			printf("%c", aln2[i+j]);
		}
		printf("\n\n");
		// recal controls
		m = m - 50;
		n++;
	}

} // end of print_align

int main (int argc, char *argv[]) {

	char	seq1[ MAXSEQ ];                  // input seq1
	char	seq2[ MAXSEQ ];                  // input seq2
	char	bts1[ MAXSEQ ];                  // back trace seq1
	char	bts2[ MAXSEQ ];                  // back trace seq2
	char	aln[ MAXSEQ ];                   // the one line alignment string
	int		nwas;                            // needleman-wunsch alginment score
	int		c;                               // options
	char	infile1[ MAXFNL ];               // input file1
	char	infile2[ MAXFNL ];               // input file2
	int	    format;                          // output controls

	// options
	while ( ( c = getopt( argc, argv, "ha:b:i:j:srf" )  ) != EOF )
		
		switch (c) {
			case 'h':
				usage();
				break;
			case 'a':
				strcpy(seq1, optarg);
				break;
			case 'b':
				strcpy(seq2, optarg);
				break;
			case 'i':
				strcpy(infile1, optarg);
				break;
			case 'j':
				strcpy(infile2, optarg);
				break;
			case 's':
				format = 1;
				break;
			case 'r':
				format = 2;
				break;
			case 'f':
				format = 3;
				break;

		}

	if ( argc == 1 ) {
	
		usage();

	}

	// check if sequences were passed by file
	if ( ( infile1[0] != '\0' ) && ( infile2[0] != '\0' ) ) {

		readseq( infile1, seq1 );
		readseq( infile2, seq2 );

	}

	// check if sequences are ready
	if ( ( seq1[0] == '\0' ) || ( seq2[0] == '\0' ) ) {

		printf("Error: insufficent input sequences found!\n");
		usage();

	}

	// covert sequences to UPPER case
	upper( seq1 );
	upper( seq2 );

	// do needleman-wunsch alignment
	nwas = nwalign( seq1, seq2, bts1, bts2, aln );

	// output alignment
	if ( format == 1 ) {

		printf("%d", nwas);

	} else if ( format == 2 ) {

		printf("%d\n", nwas);
		puts(bts1);
		puts(aln);
		puts(bts2);

	} else if ( format == 3 ) { 

		printf("Needleman-Wunsch Alignment:\n\n");
		printf("Score:\t%d\n\n", nwas);
		print_align( bts1, bts2, aln );
	
	} else {
	
		printf("Unknown output format!\n\n");
		usage();

	}

}
