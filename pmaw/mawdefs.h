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
#include "stack.h"
#include <list>
#define ALLOC_SIZE              1048576
#define DEL                     '$'
#define DEL_STR                 "$"

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define IUPAC                   "ACGTUWSMKRYBDHVN"          	//IUPAC alphabet
#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace sdsl;
using namespace std;

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_maw_filename;
   INT                  k;
   INT                  K;
   INT	                r;
   INT                  total_length;
   INT                  threads;
 };

class Maw{
	public :
 	char letter;
	int start;
	int end;

	Maw(char letter, int start, int end);
};

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
double gettime( void );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );
INT RevComStr ( unsigned char * str, unsigned char * str2, INT iLen );
unsigned char Mapping( INT a );
INT RevMapping ( unsigned char b );
INT compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw );
INT Empty_stack(
                        TStack * lifo_att,
                        TStack * lifo_pos,
                        TStack * lifo_set,
                        INT ** TSetletter,
                        INT sigma,
                        INT n,
                        INT k,
                        INT K,
                        unsigned char* seq,
			std::stringstream & buff_maw,
                        INT * Tx,
                        INT * Tleft,
                        INT * Tright,
                        INT * LCP,
                        INT * SA);


INT Get_min_abs_w(
        unsigned char * seq,
        unsigned char * seq_id,
        INT * SA,
        INT n,
        INT sigma,
        INT * LCP,        
	rrr_vector<>::select_1_type max_sel,
        rrr_vector<>::rank_1_type max_rank,
        INT nb_max,
        INT k,
        INT K,
        char * fd_maw,
        INT r, 
        INT threads );

INT GetMax(INT n, INT sigma, INT* LCP, bit_vector * max_loc);
bool Next_pos(TStack * lifo_att, INT sigma, INT * Tx, INT * Tleft, INT * Tright, INT ** TSetletter, INT n, INT k, INT * LCP, TStack * lifo_pos, TStack * lifo_set, INT * SA, unsigned char * seq);
