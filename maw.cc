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
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq    = NULL;         	// the sequence in memory
        unsigned char * seq_id = NULL;         	// the sequence id in memory
        unsigned int    num_seqs = 0;           // the total number of sequences considered
	char *          alphabet;               // the alphabet
	unsigned int    i, j;
	unsigned int    k, K, r;	

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

		k       = sw . k;
		K       = sw . K;
		r       = sw . r;

                input_filename          = sw . input_filename;
                output_filename         = sw . output_filename;
        }


	double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	if ( strncmp (input_filename, "-", 1) == 0) {
		in_fd = stdin;
	}
	else if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
	c = fgetc( in_fd );
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			unsigned int max_alloc_seq_id = 0;
			unsigned int seq_id_len = 0;
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				seq_id[ seq_id_len++ ] = c;
			}
			seq_id[ seq_id_len ] = '\0';
			
		}

		unsigned int max_alloc_seq = 0;
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

			if ( r )
			{
				unsigned char * rc_seq = ( unsigned char * ) calloc ( ( seq_len + 1 ) , sizeof( unsigned char ) );
				RevComStr ( seq, rc_seq , seq_len );
				rc_seq[seq_len] = '\0';
				seq = ( unsigned char * ) realloc ( seq,   ( 2 * seq_len + 2 ) * sizeof ( unsigned char ) );
				strcat ( ( char * ) seq, DEL_STR );
				strcat ( ( char * ) seq, ( char * ) rc_seq );
				free ( rc_seq );

				#ifdef _USE_32
				fprintf( stderr, "Processing both strands of sequence %s of length %d\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing both strands of sequence %s of length %ld\n", ( char * ) seq_id, seq_len );
				#endif
			}
			else
			{
				#ifdef _USE_32
				fprintf( stderr, "Processing 5'->3' strand of sequence %s of length %d\n", ( char * ) seq_id, seq_len );
				#endif

				#ifdef _USE_64
				fprintf( stderr, "Processing 5'->3' strand of sequence %s of length %ld\n", ( char * ) seq_id, seq_len );
				#endif
			}
			compute_maw ( seq, seq_id, sw );
		}
		free ( seq );
		seq = NULL;
		free ( seq_id );
		seq_id = NULL;
		
	} while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	double end = gettime();

        fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs\n", num_seqs, ( end - start ) );
	
        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );

	return ( 0 );
}
