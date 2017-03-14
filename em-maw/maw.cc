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
#include "stream.h"
#include "utils.h"
#include <omp.h>

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
#endif



int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name
        unsigned char * seq_id = NULL;         	// the sequence id in memory
	unsigned char * seq_idbis = NULL;       // the sequence id in memory with no specials character
        unsigned int    num_seqs = 0;           // the total number of sequences considered
	char *          alphabet;               // the alphabet
	unsigned int    i, j;
	unsigned int    k, K, r;
	INT ram_use;

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
        	ram_use			= sw . ram_use;
        }
	if (ram_use > (1<<20))
		cout << "Memory allowed "<< (ram_use>>20) <<"MiB"<<endl;
	else if (ram_use > (1<<10) )
		cout << "Memory allowed "<< (ram_use>>10) << "KB" <<endl;
	else    cout << "Memory allowed "<< ram_use <<"Bytes" << endl;
	double start = gettime();
    	INT n;
	/* Read the (Multi)FASTA file in memory */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

    unsigned char seq1strand_fname[strlen( output_filename)+20] ;
    unsigned char seqfinal_fname[strlen(output_filename)+20] ;
    sprintf((char*)seq1strand_fname, "%s_r%d_1strand.txt", output_filename, sw.r);
    sprintf((char*)seqfinal_fname, "%s_r%d_final.txt", output_filename, sw.r);
    
	char c;
	unsigned char last_c;
	c = fgetc( in_fd );
	last_c=c;
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
			unsigned int seq_idbis_len=0;
			bool space=false;
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id )
				{
					seq_id = ( unsigned char * ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					if (!space)
						seq_idbis = ( unsigned char * ) realloc ( seq_idbis,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id += ALLOC_SIZE;
				}
				if (c==' ') space=true;
				seq_id[ seq_id_len++ ] = c;
				if (strchr ( "|&#!;<>(){}$~[]?", c )) c='_'; 
				if (c=='\'' ||c=='\"'||c=='\\') c='_';
				if (space) continue;
				seq_idbis[ seq_idbis_len ++ ] = c;
			}
			last_c=c;
			seq_id[ seq_id_len ] = '\0';
			seq_idbis[ seq_idbis_len ] = '\0';
		}	

	    
	    INT seq_len = 0;
            stream_writer<unsigned char> * fseq_inw=new stream_writer<unsigned char>((const char*)seq1strand_fname, ram_use);
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
                
                if( strchr ( alphabet, c ) )
                {
                    seq_len++;
		    last_c=c;
                    fseq_inw->write(c);
                }
                else
                {
                    if ( strchr ( IUPAC, c ) )
                    {
                        seq_len++;
			last_c='N';
                        fseq_inw->write('N');
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
                
		delete (fseq_inw);

		stream_reader<unsigned char> * fseq_inr = new stream_reader<unsigned char> ((const char*)seq1strand_fname,ram_use/2);
                stream_writer<unsigned char> *fseq_final = new stream_writer<unsigned char>((const char*)seqfinal_fname, ram_use/2);
                unsigned char cb;
		for (INT i=0; i<seq_len; i++)
                {
                    cb=fseq_inr->read();
                    fseq_final->write(cb);
                }
		delete(fseq_inr);
                
		if ( r )
                {
		   stream_reader<unsigned char> *fseq_bis= new stream_reader<unsigned char>((const char*) seq1strand_fname, ram_use/2);
			
		   fseq_bis->goto_end(seq_len);
		   fseq_final->write('$');                    
                    for (INT i=0; i<seq_len; i++)
                    {
                        cb=fseq_bis->read_reverse();
                        fseq_final->write(RevComChar(cb));
                    }
		    delete(fseq_bis);
                }
		delete(fseq_final);
                remove((const char*) seq1strand_fname);
                if (r)
                {
                    n=seq_len*2+1;
#ifdef _USE_32
                    fprintf( stderr, "Processing both strands of sequence %s of length %d\n", ( char * ) seq_id, seq_len );
#endif
                    
#ifdef _USE_64
                    fprintf( stderr, "Processing both strands of sequence %s of length %ld\n", ( char * ) seq_id, seq_len );
#endif
                }
                else
                {
                    n=seq_len;
#ifdef _USE_32
                    fprintf( stderr, "Processing 5'->3' strand of sequence %s of length %d\n", ( char * ) seq_id, seq_len );
#endif
                    
#ifdef _USE_64
                    fprintf( stderr, "Processing 5'->3' strand of sequence %s of length %ld\n", ( char * ) seq_id, seq_len );
#endif
                }
                
                if (sw.c==1)
		{
                    /* Run pSAscan */
                    char sa_fname[strlen(input_filename)+strlen((const char *) seq_idbis)+20] ;
		    sprintf(sa_fname, "%s_%s_r%d_SA.sa5", input_filename, seq_idbis, sw.r);
		    char commandesa[ strlen(sa_fname) + 1000 ];
		    sprintf(commandesa, "./pSAscan-0.1.0/src/psascan %s -m %ldL -o %s", seqfinal_fname, ram_use>>20, sa_fname);
		    int outsa=system(commandesa);

		    std::cout<<"BWT construction"<<std::endl; 
		    /*Construct BWT*/
                    char bwt_fname[strlen(input_filename)+strlen((const char*)seq_idbis)+20] ;
		    sprintf(bwt_fname, "%s_%s_r%d_BWT.bwt5", input_filename, seq_idbis, sw.r);
                    compute_bwt((char*)seqfinal_fname,sa_fname,bwt_fname,ram_use,n );

                    /* Run LCPscan */
                    char lcp_fname[strlen(input_filename)+strlen((const char*)seq_idbis)+20] ;
		    sprintf(lcp_fname, "%s_%s_r%d_LCP.lcp5", input_filename, seq_idbis, sw.r);
                    char commande[strlen(sa_fname) + strlen(lcp_fname) + 1000];
		    sprintf(commande,"./EM-SparsePhi-0.1.0/src/construct_lcp_parallel -m %ldL -o %s -s %s %s", ram_use>>20, lcp_fname, sa_fname,seqfinal_fname);
 		    int out=system(commande);
		}
		compute_maw (n,last_c,seq_idbis,seq_id, sw,seqfinal_fname );
        	remove((const char*) seqfinal_fname);
            }
	    else
	    {
		delete(fseq_inw);remove((const char*) seq1strand_fname);
	    }
        
	   free ( seq_id );
	   free ( seq_idbis );
	   seq_id = NULL;
           seq_idbis = NULL;
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

