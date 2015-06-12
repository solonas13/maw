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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "mawdefs.h"

static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "min-length",              required_argument, NULL, 'k' },
   { "max-length",              required_argument, NULL, 'K' },
   { "circular",                required_argument, NULL, 'c' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL, 0   }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> k                              = 2;
   sw -> K				= 10;
   sw -> c                              = 0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:k:K:c:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           args ++;
           break;

         case 'K':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> K = val;
           args ++;
           break;

         case 'c':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> c = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 5 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: sc-maw <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'\n"
                     "                                      for protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename with the distance matrix.\n" );
   fprintf ( stdout, "  -k, --min-length          <int>     The min length for maws.\n");
   fprintf ( stdout, "  -K, --max-length          <int>     The max length for maws.\n");
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -c, --circular            <int>     `1' for circular sequence comparison or\n"
                     "                                      `0' otherwise (default: 0).\n" );
 }

