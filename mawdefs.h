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
#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

using namespace sdsl;
using namespace std;

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   unsigned int         k;
   unsigned int         K;
   unsigned int         r;
   unsigned int         total_length;
 };

double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, unsigned int iLen );
unsigned int compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw );
unsigned char Mapping( int a );
int RevMapping ( unsigned char b );

unsigned int GetBefore (
				unsigned char * seq,
                                int n,
                                int sigma,
				int * SA,
                                int * LCP,
                                bit_vector * Before,
                                bit_vector * Beforelcp );
unsigned int GetMaws(
				unsigned char * seq,
				unsigned char * seq_id,
				int * SA,
				int n,
				int sigma,
				int * LCP,
				bit_vector * Before,
				bit_vector * Beforelcp,
				unsigned int k,
				unsigned int K,
				char * out_file );
