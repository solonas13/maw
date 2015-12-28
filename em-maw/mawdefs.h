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

#include <sdsl/bit_vectors.hpp>
#include "stream.h"
#include "utils.h"
#include "uint40.h"
#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet

using namespace sdsl;
using namespace std;

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

struct TSwitch
{
    char               * alphabet;
    char               * input_filename;
    char               * output_filename;
    unsigned int         k;
    unsigned int         K;
    unsigned int         r;
    unsigned int         c;
    long		 ram_use;
    unsigned int         total_length;
};

struct SApairs
{
    INT  pos;
    INT  SApos;
    inline int compare(SApairs * s2)
    {
        if(SApos > s2->SApos ) return 1;
        else
        {
            if (SApos==s2->SApos) return 0;
            else return -1;
        }
    }
};

struct BWTpairs
{
    INT  pos;
    unsigned char  bwt;
    inline int compare(BWTpairs * s2)
    {
        if(pos > s2->pos ) return 1;
        else
        {
            if (pos==s2->pos) return 0;
            else return -1;
        }
    }
};

unsigned char RevComChar ( unsigned char c );
unsigned int compute_bwt( char* seq_fname, char* sa_fname, char* bwt_fname, long ram_use, INT n );
static void swap (SApairs * a, INT i, INT j);
static void swap (BWTpairs * a, INT i, INT j);
static INT partition( SApairs * a, INT l, INT r );
static INT partition( BWTpairs * a, INT l, INT r );
static void quickrec(SApairs * a,INT l, INT r);
static void quickrec(BWTpairs * a, INT l, INT r);
double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, INT iLen );
unsigned int compute_maw ( INT n, unsigned char c, unsigned char * file_id, unsigned char * seq_id, struct TSwitch sw );
unsigned char Mapping( int a );
int RevMapping ( unsigned char b );
unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );


unsigned int GetBefore (
                        stream_reader<unsigned char> * fBWT,
                        unsigned char c,
                        INT n,
                        int sigma,
                        stream_reader<uint40> * fSA,
                        stream_reader<uint40> * fLCP,
                        char * Before_fname,
			char * Before_fname2,
                        char * Beforelcp_fname,
			char * Beforelcp_fname2,
                        long ram_use);
unsigned int GetMaws(
                                unsigned char * seq_id,
                                stream_reader<uint40> * SA,
                                INT n,
                                int sigma,
                                stream_reader<uint40> * LCP,
                                char * Before_fname,
                                char * Beforelcp_fname,
                                unsigned int k,
                                unsigned int K,
                                char * out_file,
				char * fout,
				int r,
				long ram_use );
