/**
    MAW: Minimal Absent Words
    Copyright (C) 2014 Alice Heliou and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "mawdefs.h"

int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq    = NULL;         // the sequence in memory
        unsigned char * seq_id = NULL;         // the sequence id in memory
        unsigned char ** seqs    = NULL;         // the sequences in memory
        unsigned char ** seqs_id = NULL;         // the sequences id in memory
        unsigned int    num_seqs = 0;           // the total number of sequences considered
	char *          alphabet;               // the alphabet
	unsigned int    i, j;

	double ** D;
	double * buf;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 9 )
        {
                usage ();
                return ( 1 );
        }
        else
        {
                if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet = ( char * ) DNA;
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet = ( char * ) PROT;
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
                        return ( 1 );
                }

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;
        }


	double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	c = fgetc( in_fd );
	do
	{
		unsigned int seq_id_len = 0;
		unsigned int max_alloc_seq_id = 0;
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				if( c == ' ' ) continue;
				else	seq_id[ seq_id_len++ ] = c;
			}
			seq_id[ seq_id_len ] = '\0';
			
		}

		INT max_alloc_seq = 0;
		INT seq_len = 0;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ seq_len++ ] = c;
			}
			else
			{
				if ( strchr ( IUPAC, c ) )
				{
					seq[ seq_len++ ] = 'N';
				}
				else if ( c == ' ' )
				{
				}
				else
				{
					fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
					return ( 1 );
				}
			}

		}
		if( seq_len != 0 )
		{
			num_seqs++;
			if ( seq_len >= max_alloc_seq )
			{
				seq = ( unsigned char * ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq += ALLOC_SIZE;
			}
			seq[ seq_len ] = '\0';

			if ( sw . c )
			{
				if ( sw . K > seq_len ) sw . K = seq_len;
 
				seq = ( unsigned char * ) realloc ( seq,   ( seq_len + seq_len + 1 ) * sizeof ( unsigned char ) );
				unsigned char * seq2 = ( unsigned char * ) strdup ( ( const char * ) seq );
				strcat ( ( char * ) seq, ( char * ) seq2 );
				seq_len += seq_len;
				free ( seq2 );

				#ifdef _USE_32
				fprintf( stderr, "Processing circular sequence %s of length %d.\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing circular sequence %s of length %ld.\n", ( char * ) seq_id, seq_len );
				#endif
			}
			else
			{
				#ifdef _USE_32
				fprintf( stderr, "Processing linear sequence %s of length %d.\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing linear sequence %s of length %ld.\n", ( char * ) seq_id, seq_len );
				#endif
			}
			seq[ seq_len ] = '\0';

			seqs = ( unsigned char ** ) realloc ( seqs,   ( num_seqs ) * sizeof ( unsigned char * ) );
			seqs[num_seqs - 1] = ( unsigned char * ) calloc ( seq_len + 1,   sizeof ( unsigned char ) );
			strcat ( ( char * ) seqs[ num_seqs - 1 ], ( char * ) seq );
			seqs[num_seqs - 1][ seq_len ] = '\0';
			seqs_id = ( unsigned char ** ) realloc ( seqs_id,   ( num_seqs ) * sizeof ( unsigned char * ) );
			seqs_id[num_seqs - 1] = ( unsigned char * ) calloc ( seq_id_len + 1,   sizeof ( unsigned char ) );
			strcat ( ( char * ) seqs_id[ num_seqs - 1 ], ( char * ) seq_id );
			seqs_id[num_seqs - 1][ seq_id_len ] = '\0';
 
		}
		free ( seq );
		seq = NULL;
		free ( seq_id );
		seq_id = NULL;
		
	} while( c != EOF );
	
#if 0
	for ( i = 0; i < NmawX; i ++ )
	  fprintf( stderr, "<%c, %d, %d>\n", ( char ) mawX[i] . letter, mawX[i] . pos, mawX[i] . size );
	
	fprintf( stderr, "-----------\n" );
	
	for ( i = 0; i < NmawY; i ++ )
	  fprintf( stderr, "<%c, %d, %d>\n", ( char ) mawY[i] . letter, mawY[i] . pos, mawY[i] . size );
#endif
	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( num_seqs == 1 )
	{
		fprintf ( stderr, " Error: the number of input sequences must be greater than 1!\n" );
		return ( 1 );
	}

	/* 2d dynamic memory allocation of the Distance Matrix */
	D = ( double ** ) malloc ( ( num_seqs ) * sizeof ( double * ) );
   	if ( D == NULL )
    	{
      		fprintf ( stderr, " ERROR: DM matrix could not be allocated!!!\n" );
      		return ( 1 );
    	} 
	buf =  ( double * ) calloc ( ( ( size_t ) ( num_seqs ) ) * ( ( size_t ) ( num_seqs ) ), sizeof ( double ) ); 
   	if ( buf == NULL )
    	{
      		fprintf ( stderr, " ERROR: DM matrix could not be allocated!!!\n" );
      		return ( 1 );
    	} 
        for ( i = 0; i < num_seqs; ++ i ) 	D[i] = &buf[( size_t ) i * ( size_t ) ( num_seqs ) ];
	
	fprintf( stderr, "Computing minimal absent words and making the comparison.\n" );
	#pragma omp parallel for
	for ( int i = 0; i < num_seqs; i++ )	
	{
		TMaw * mawX = NULL; 
		unsigned int NmawX = 0;	
		compute_maw ( seqs[i], seqs_id[i], sw, &mawX, &NmawX );

		for ( int j = i; j < num_seqs; j++ )	
		{
			TMaw * mawY = NULL;
			unsigned int NmawY = 0;

			compute_maw ( seqs[j], seqs_id[j], sw, &mawY, &NmawY );

			double XY_distance = 0;
			if ( i == j )	D[i][j] = 0.0;
			else    
			{        
				//edit_distance ( seqs[i], seqs[j], &XY_distance );
				maw_seq_comp ( seqs[i], mawX, &NmawX, seqs[j], mawY, &NmawY, &XY_distance );
				D[i][j] = D[j][i] = XY_distance;
			}
			free ( mawY );
		}
		free ( mawX );
	}

	double end = gettime();

	if ( ! ( out_fd = fopen ( sw . output_filename, "w" ) ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	fprintf( out_fd, "%d\n", num_seqs );

	for ( i = 0; i < num_seqs; i++ )	
	{
		fprintf( out_fd, "%s\t", seqs_id[i] );
		for ( j = 0; j < num_seqs; j++ )	
		{
			fprintf( out_fd, "%lf\t", D[i][j]);
		}
		fprintf( out_fd, "\n");
	}

	fclose ( out_fd );

        fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs.\n", num_seqs, ( end - start ) );
	
	for ( i = 0; i < num_seqs; i ++ )
		free ( seqs[i] );
	free ( seqs);
	for ( i = 0; i < num_seqs; i ++ )
		free ( seqs_id[i] );
	free ( seqs_id );

        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );
	free ( buf );
        free ( D );

	return ( 0 );
}
