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
#include <sys/time.h>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include "mawdefs.h"
#include <bitset>

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                           // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors
#include "stackinfile.h"
#include "triplet.h"
using namespace sdsl;
using namespace std;

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};


bool BWTcompare(BWTpairs l, BWTpairs r){return l.pos<r.pos;}
bool SAcompare(SApairs l, SApairs r ){return l.SApos<r.SApos;}

unsigned int compute_bwt( char* seq_fname, char* sa_fname, char* bwt_fname, long ram_use, INT n )
{
    if (ram_use > 4*sizeof(uint40)+n*sizeof(unsigned char))
    {
	FILE *fseq = fopen(seq_fname, "r");
	unsigned char * Seq =new unsigned char [n];
	INT n_read=fread(Seq, sizeof(unsigned char), n, fseq);
	if (n_read!=n)
	{
		std::cout<<"Problem while reading the sequence during BWT computation"<<std::endl;
    }
	fclose(fseq);
	long ram_left= ram_use - n*sizeof(unsigned char);
	stream_reader<uint40>* fSA =  new stream_reader <uint40> (sa_fname, ram_left/2);
	stream_writer<unsigned char> * fBWT = new stream_writer <unsigned char>(bwt_fname, ram_left/2);
	INT pos;
	for (INT i=0; i<n; i++)
	{
		pos=fSA->read();
		if (pos>0)
			fBWT->write(Seq[pos-1]);
		else
			fBWT->write('$');
	}
	delete(fSA);
	delete(fBWT);
	delete [] Seq;
    } 
    else 
    {
	FILE * fSA = fopen(sa_fname, "r");
	FILE * fBWT = fopen(bwt_fname, "w");
	stream_reader<unsigned char> * fseq = new stream_reader<unsigned char>(seq_fname, ram_use/2);
    	INT elems = ram_use / (8*sizeof(uint40));
    	uint40 * SAvalue = new uint40 [elems];
    	SApairs * sapairs = new SApairs [elems];
    	unsigned char * BWTvalue = new unsigned char [elems];
    	BWTpairs * bwtpairs= new BWTpairs [elems];
    	INT pos=0;
        unsigned char c;
    	INT nbloop=n/elems +1;
    	INT sum=0;
	bool b=false;
    	for (INT i=0; i<nbloop; i++)
    	{
        	if (i==nbloop-1)
       		{
            		elems=n-i*elems;
        	}
        	INT n_read=fread(SAvalue, sizeof(uint40), elems, fSA);
                if (n_read!=elems)
                {
                    std::cout<<"Problem while reading the sequence during BWT computation"<<std::endl;
                }
        	for (INT j=0; j<elems; j++)
        	{
            		sapairs[j].pos=j;
            		sapairs[j].SApos=SAvalue[j];
        	}
       		std::sort(sapairs,sapairs+elems,SAcompare);
        	pos=0;
                fseq->goto_set();
        	for (INT j=0; j<elems; j++)
        	{
	    		if (sapairs[j].SApos==0)
	    		{
				bwtpairs[j].pos=sapairs[j].pos;
				bwtpairs[j].bwt='$';
				sum++;
	    		}
			else
			{
	    	  	    while (pos<sapairs[j].SApos)
	    		    {
            			c=fseq->read();
				pos++;
	    		    }
           		    if (pos==sapairs[j].SApos)
            		    {
	        		bwtpairs[j].pos=sapairs[j].pos;
                		bwtpairs[j].bwt=c;
				sum++;
            		    }
			}
        	}
        	std::sort(bwtpairs, bwtpairs+elems, BWTcompare);
       		for (INT j=0; j<elems; j++)
        	{
	    		BWTvalue[j]=bwtpairs[j].bwt;
        	}
       	 	fwrite(BWTvalue, sizeof(unsigned char), elems, fBWT);
	}
	delete[] SAvalue;
        delete[] sapairs;
        delete[]  BWTvalue;
        delete[] bwtpairs;

     	delete(fseq);
     	fclose(fSA);
     	fclose(fBWT);
    }
}

unsigned int compute_maw (INT n,unsigned char c,unsigned char * file_id,  unsigned char * seq_id, struct TSwitch sw )
{
   int sigma;

	if      ( ! strcmp ( "DNA", sw . alphabet ) )   sigma = strlen ( ( char * ) DNA );
    else if ( ! strcmp ( "PROT", sw . alphabet ) )  sigma = strlen ( ( char * ) PROT );

    long ram_use =sw.ram_use;    
    char str[100] ;
    sprintf(str, "%s_%s_r%d_BWT.bwt5",sw.input_filename, file_id, sw.r);

    stream_reader<unsigned char> *fBWT = new stream_reader<unsigned char>(str,ram_use/6);
    sprintf(str, "%s_%s_r%d_SA.sa5",sw.input_filename, file_id,sw.r);
    stream_reader<uint40> *fSA = new stream_reader<uint40>(str,ram_use/6);
    sprintf(str, "%s_%s_r%d_LCP.lcp5",sw.input_filename, file_id,sw.r);
    stream_reader<uint40> *fLCP = new stream_reader<uint40>(str,ram_use/6);
    
    char Before_fname[100];
    sprintf(Before_fname, "%s_Before_step1.txt",sw.output_filename);
    char Beforelcp_fname[100];
    sprintf(Beforelcp_fname, "%s_Beforelcp_step1.txt", sw.output_filename);
    char Before_fname2[100];
    sprintf(Before_fname2, "%s_Before_step2.txt",sw.output_filename);
    char Beforelcp_fname2[100];
    sprintf(Beforelcp_fname2, "%s_Beforelcp_step2.txt", sw.output_filename);
    GetBefore ( fBWT,c, n , sigma, fSA, fLCP, Before_fname, Before_fname2, Beforelcp_fname,Beforelcp_fname2, ram_use/2);
    remove(Before_fname);
    remove(Beforelcp_fname);
    char fout[100];
    sprintf(fout, "%s_compressed.txt", sw.output_filename);
    GetMaws( seq_id, fSA, n, sigma, fLCP, Before_fname2, Beforelcp_fname2, sw . k, sw . K, sw . output_filename,fout, sw.r,sw.f, ram_use/2 );
    remove(Before_fname2);
    remove(Beforelcp_fname2);
    delete ( fSA );
    delete ( fBWT );
    delete( fLCP );

	return ( 1 );
}


unsigned char Mapping( int a )
{
	char c = DEL;
        switch ( a )
	{
            case 0:
                c = 'A';
                break;
            case 1:
                c = 'C';
                break;
            case 2:
                c = 'G';
                break;
            case 3:
                c = 'T';
                break;
            case 4:
                c = 'N';
                break;
            case 5:
                c = 'R';
                break;
            case 6:
                c = 'D';
                break;
            case 7:
                c = 'Q';
                break;
            case 8:
                c = 'E';
                break;
            case 9:
                c = 'H';
                break;
            case 10:
                c = 'I';
                break;
            case 11:
                c = 'L';
                break;
            case 12:
                c = 'K';
                break;
            case 13:
                c = 'M';
                break;
            case 14:
                c = 'F';
                break;
            case 15:
                c = 'P';
                break;
            case 16:
                c = 'S';
                break;
            case 17:
                c = 'W';
                break;
            case 18:
                c = 'Y';
                break;
            case 19:
                c = 'V';
                break;
        }
	return ( c );
}

int RevMapping ( unsigned char b )
{
	int a = -1;
        switch ( b )
	{
            case 'A':
                a = 0;
                break;
            case 'C':
                a = 1;
                break;
            case 'G':
                a = 2;
                break;
            case 'T':
                a = 3;
                break;
            case 'N':
                a = 4;
                break;
            case 'R':
                a = 5;
                break;
            case 'D':
                a = 6;
                break;
            case 'Q':
                a = 7;
                break;
            case 'E':
                a = 8;
                break;
            case 'H':
                a = 9;
                break;
            case 'I':
                a = 10;
                break;
            case 'L':
                a = 11;
                break;
            case 'K':
                a = 12;
                break;
            case 'M':
                a = 13;
                break;
            case 'F':
                a = 14;
                break;
            case 'P':
                a = 15;
                break;
            case 'S':
                a = 16;
                break;
            case 'W':
                a = 17;
                break;
            case 'Y':
                a = 18;
                break;
            case 'V':
                a = 19;
                break;
        }
	return ( a );
}

/* computes the reverse complement of str */
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, INT iLen )
{
   INT i = 0;
   while ( iLen -- )
    {
      switch ( str[iLen] )
       {
         case 'A':
           str2[i++] = 'T';
           break;
         case 'C':
           str2[i++] = 'G';
           break;
         case 'G':
           str2[i++] = 'C';
           break;
         case 'T':
           str2[i++] = 'A';
           break;
         case 'N':
           str2[i++] = 'N';
           break;
         default:
           return ( 0 );
       }
    }
   return ( 1 );
}

/* computes the reverse complement of char */
unsigned char RevComChar ( unsigned char c )
{
    switch ( c )
    {
        case 'A':
        return 'T';
        break;
        case 'C':
        return 'G';
        break;
        case 'G':
        return 'C';
        break;
        case 'T':
        return 'A';
        break;
        case 'N':
        return 'N';
        break;
    }
    
    return '\0' ;
}



unsigned int GetBefore (
                        stream_reader<unsigned char>  * fBWT,
                        unsigned char c,
                        INT n,
                        int sigma,
                        stream_reader<uint40>  * fSA,
                        stream_reader<uint40>  * fLCP,
                        char * Before_fname,
			char * Before_fname2,
                        char * Beforelcp_fname,
			char * Beforelcp_fname2,
			long ram_use)
{
    INT hm = 0;
    INT k = 0;
    INT lcp;
    INT mem;
    INT proxa;
    INT proxb;
    INT width_to_print=floor(log10(n))+2;
    // ram_use/4
    TStackinfile lifo_lcp;
    StackNew ( &lifo_lcp, sizeof( INT ),"streambis_stacklcp.txt", ram_use/8 );
    TStackinfile lifo_mem;
    StackNew ( &lifo_mem, sizeof( INT ),"streambis_stackmem.txt", ram_use/16 );
    TStackinfile lifo_rem;
    StackNew ( &lifo_rem, sizeof( INT ),"streambis_stackrem.txt", ram_use/16 );
    // ram_use/4
    TStackinfile lifo_interval;
    StackNew ( &lifo_interval, sizeof( bitset<8> ),"streambis_intervallcp.txt", ram_use/8 );
    TStackinfile lifo_intervalmem;
    StackNew ( &lifo_intervalmem, sizeof( bitset<8> ),"streambis_intervalmem.txt", ram_use/16 );
    TStackinfile lifo_intervalrem;
    StackNew ( &lifo_intervalrem, sizeof( bitset<8> ),"streambis_intervalrem.txt", ram_use/16 );
    // ram_use/2
    stream_writer<char> * fBefore1 = new stream_writer<char> (Before_fname, ram_use/4);
    stream_writer<char> * fBeforelcp1 = new stream_writer<char> (Beforelcp_fname, ram_use/4);	


    lcp = 0;
    StackPush(&lifo_lcp, &lcp);
    
    fLCP->goto_set();
    INT * LCPmem= new INT [2];
    fLCP->goto_set();
    fSA->goto_set();
    fBWT->goto_set();
    
    bitset<8> interval;
    bitset<8> interval_mem;
    bitset<8> interval_lcp;
    

    interval_lcp[sigma-1-RevMapping(c)]=1;
    StackPush(&lifo_interval, &interval_lcp);
    bitset<8> * Beforemem = new bitset<8> [3];
    if( ( Beforemem == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for Before.\n" );
        return ( 0 );
    }
    bitset<8>* Beforelcpmem = new bitset<8> [3];
    if( ( Beforelcpmem == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for BeforeLCP.\n" );
        return ( 0 );
    }
    INT * SAmem= new INT [2];
    SAmem[0]=0;
    SAmem[1]=fSA->read();
    LCPmem[0]=0;
    LCPmem[1]=fLCP->read();
    unsigned char * BWTmem = new unsigned char [2];
    BWTmem[0]=' ';
    BWTmem[1]=fBWT->read();
    // First pass : top-down
    for ( INT i = 0; i < n; i++ )
    {
        if (i>0)
        {
            SAmem[0]=SAmem[1];
            SAmem[1]=fSA->read();
            LCPmem[0]=LCPmem[1];
            LCPmem[1]=fLCP->read();
            BWTmem[0]=BWTmem[1];
            BWTmem[1]=fBWT->read();
        }
        // first we update the interval table
        // we empty the interval that corresponds to a higher lcp value
	if ( i > 0 && LCPmem[1] < LCPmem[0])
        {
            StackPop(&lifo_lcp,&lcp);
            StackPop(&lifo_interval, &interval_lcp);
            
            while(!StackEmpty(&lifo_lcp)&&lcp>LCPmem[1])
            {
                StackPop(&lifo_lcp,&mem);
                StackPop(&lifo_interval, &interval_mem);
                if (mem <=LCPmem[1])
                {
                        if (mem!=LCPmem[1]){interval=interval_lcp;}     //initialisation of the next intervals if it hasn't been open
                        Beforemem[0]=interval_lcp;
                        Beforelcpmem[0]=interval_lcp;
			
                        if (mem==LCPmem[1]){ Beforelcpmem[0]=interval_mem; interval=interval_mem;}
                }
                lcp=mem;
                interval_lcp=interval_mem;
            }
            StackPush(&lifo_lcp,&lcp);
            StackPush(&lifo_interval, &interval_lcp);
            
        }
	StackPop(&lifo_interval, &interval_lcp);
	StackPush(&lifo_interval, &interval_lcp);
	if (LCPmem[0]==LCPmem[1])
		interval=interval_lcp;
        // we update those having a lower lcp
        if ( SAmem[1] - 1 >= 0 )
        {
            k = RevMapping(BWTmem[1]);
        }
        else
            k = - 1;
        if ( k != -1 )
        {
            while(!StackEmpty(&lifo_lcp))
            {
                StackPop(&lifo_lcp,&lcp);
                StackPop(&lifo_interval, &interval_lcp);
                if (interval_lcp[sigma-k-1]==1)
		{
		    StackPush(&lifo_mem, &lcp);
                    StackPush(&lifo_intervalmem, &interval_lcp);
		    break;
		}
                interval_lcp[sigma-k-1]=1;
                
                StackPush(&lifo_mem, &lcp);
                StackPush(&lifo_intervalmem, &interval_lcp);
            }

            while (!StackEmpty(&lifo_mem))
            {        
                StackPop(&lifo_mem,&lcp);
                StackPop(&lifo_intervalmem, &interval_lcp);
                
                StackPush(&lifo_lcp, &lcp);
                StackPush(&lifo_interval, &interval_lcp);
            }
            interval[sigma-1-k]=1;
        }
       
        if ( ( i - 1 ) >= 0 && SAmem[0] - 1 >= 0 )
            if (i>0 && LCPmem[1]>0 && RevMapping(BWTmem[0])!=-1) // in this case we also add the letter preceding the last suffix
            {
                interval[sigma-1-RevMapping(BWTmem[0])]=1;
            }
        Beforelcpmem[1]=interval;
        if (k!=-1)
        {
            Beforemem[1][sigma-k-1]=1;
            Beforemem[2][sigma-k-1]=1;
            Beforelcpmem[1][sigma-k-1]=1;
            Beforelcpmem[2][sigma-k-1]=1;
        }
        StackPop(&lifo_lcp,&lcp);
	StackPop(&lifo_interval, &interval_lcp);
        if(lcp!=LCPmem[1])  // no duplicates
        {
            StackPush(&lifo_lcp,&lcp);
            StackPush(&lifo_interval, &interval_lcp);
            lcp=LCPmem[1];
        }
        StackPush(&lifo_lcp,&lcp);
        StackPush(&lifo_interval, &interval);
        fBefore1->write(char(Beforemem[0].to_ulong()));
        fBeforelcp1->write(char(Beforelcpmem[0].to_ulong()));
        Beforemem[0]=Beforemem[2];
        Beforelcpmem[0]=Beforelcpmem[2];
        fBefore1->write(char(Beforemem[1].to_ulong()));
        fBeforelcp1->write(char(Beforelcpmem[1].to_ulong()));
        Beforemem[1]=bitset<8>( 0<<7);
        Beforelcpmem[1]=bitset<8>(0<<7);
        Beforemem[2]=bitset<8>(0<<7);
        Beforelcpmem[2]=bitset<8>(0<<7);
	interval.reset();
	interval_mem.reset();
	interval_lcp.reset();
    }
    
    delete [] SAmem;
    delete [] BWTmem;
    fBefore1->write(char(Beforemem[0].to_ulong()));
    fBeforelcp1->write(char(Beforelcpmem[0].to_ulong()));
    INT result;
    delete(fBefore1);
    delete(fBeforelcp1);

    stream_reader<char>* fBefore2 = new stream_reader<char>(Before_fname, ram_use/8);
    stream_reader<char>* fBeforelcp2 = new stream_reader<char>(Beforelcp_fname, ram_use/8);
    stream_writer<char>* fBefore3 = new stream_writer<char>(Before_fname2, ram_use/8);
    stream_writer<char>* fBeforelcp3 = new stream_writer<char>(Beforelcp_fname2, ram_use/8);
    fBefore2->goto_end(2*n+1);
    fBeforelcp2->goto_end(2*n+1);
    //we empty the interval table
    while(!StackEmpty(&lifo_lcp))
    {
        StackPop(&lifo_lcp,&lcp);
        StackPop(&lifo_interval,&interval_lcp);
    }
    lcp=0;
    StackPush(&lifo_interval,&interval);
    StackPush(&lifo_lcp,&lcp);
    
    long  file_index;
    for (INT i=n-1; i>-1; i--)
    {
	interval.reset();
    	interval_lcp.reset();
    	interval_mem.reset();
        LCPmem[1]=LCPmem[0];
        if (i==n-1)
           LCPmem[0]=fLCP->read_reverse(true);
        else
           LCPmem[0]=fLCP->read_reverse();
        StackPop(&lifo_lcp,&lcp);
        StackPop(&lifo_interval,&interval_lcp);
        proxa=LCPmem[0]+1;   //proxb is the lcp-value that is just higher than LCP[i]
 	
	while(!StackEmpty(&lifo_lcp) && lcp>LCPmem[0])
        {
	    StackPush(&lifo_rem,&lcp);
            StackPush(&lifo_intervalrem,&interval_lcp);
            
	    StackPop(&lifo_lcp,&mem);
            StackPop(&lifo_interval,&interval_mem);
            
            if (mem <LCPmem[0])            //initialisation of the interval if it hasn't been open
            {
                interval=interval_lcp;
                proxa=lcp;
            }
            if (mem ==LCPmem[0]){proxa=lcp;interval=interval_mem;}
            lcp=mem;
            interval_lcp=interval_mem;
        }
        StackPush(&lifo_lcp,&lcp);
        StackPush(&lifo_interval, &interval_lcp);
        if (LCPmem[1]==LCPmem[0])
        	interval=interval_lcp; 
    //read Before and Beforelcp
        Beforemem[1]=bitset<8>(fBefore2->read_reverse());
        Beforelcpmem[1]=bitset<8>(fBeforelcp2->read_reverse());      
        Beforemem[0]=bitset<8>(fBefore2->read_reverse());
        Beforelcpmem[0]=bitset<8>(fBeforelcp2->read_reverse());
        // we update the lower intervals
        for (int k=0; k<sigma;k++)
        {
            if (Beforemem[0][sigma-1-k]==1)
            {
                while(!StackEmpty(&lifo_lcp))
                {
                    StackPop(&lifo_lcp,&lcp);
                    StackPop(&lifo_interval, &interval_lcp);
                    if (interval_lcp[sigma-k-1]==1)
		    {
			StackPush(&lifo_lcp, &lcp);
                    	StackPush(&lifo_interval, &interval_lcp);
		    	break;
		    }
                    interval_lcp[sigma-k-1]=1;
                    
                    StackPush(&lifo_mem, &lcp);
                    StackPush(&lifo_intervalmem, &interval_lcp);
                    
                }
                while (!StackEmpty(&lifo_mem))
                {
                    StackPop(&lifo_mem,&lcp);
                    StackPop(&lifo_intervalmem, &interval_lcp);
                    
                    StackPush(&lifo_lcp, &lcp);
                    StackPush(&lifo_interval, &interval_lcp);
                }
                interval[sigma-1-k]=1;
            }
        }
        
         Beforelcpmem[0]|=interval;
        
	if (i<n-1)
        {
            if (LCPmem[1]>LCPmem[0])
            {
                StackPop(&lifo_rem,&lcp);  // this lcp is the one that is just higher than LCP[i]
                StackPop(&lifo_intervalrem, &interval_mem);
                Beforemem[0]|=interval_mem;
                if (lcp==proxb) Beforemem[1]|=interval_mem;
                if (lcp==LCPmem[1]) Beforelcpmem[1]|=interval_mem;
		while(!StackEmpty(&lifo_rem))
                {
                    StackPop(&lifo_rem,&lcp);
                    StackPop(&lifo_intervalrem, &interval_mem);
                    if (lcp==proxb) Beforemem[1]|=interval_mem;
                    if (lcp==LCPmem[1]) Beforelcpmem[1]|=interval_mem;
                }
            }
            else
	    {
		if (LCPmem[1]==LCPmem[0])
            	{ // then nothing to do with proxb
                    Beforelcpmem[1]|=interval;
                }   
                else
                { // then proxb=LCPmem[0]
                    StackPop(&lifo_interval, &interval_lcp); // this interval is of LCP-depth LCPmem[1]
                    StackPush(&lifo_interval, &interval_lcp);
                    Beforelcpmem[1]|=interval_lcp;
                    Beforemem[1]|=interval;
            	} 
            }
        }
        proxb=proxa;
        StackPop(&lifo_lcp,&lcp);
	StackPop(&lifo_interval, &interval_lcp);
        if(lcp!=LCPmem[0])  // no duplicates
        {
            StackPush(&lifo_lcp,&lcp);
	    StackPush(&lifo_interval,&interval_lcp);
            lcp=LCPmem[0];
        }
        StackPush(&lifo_lcp,&lcp);
	StackPush(&lifo_interval, &interval);

	fBefore3->write(char(Beforemem[1].to_ulong()));
	fBefore3->write(char(Beforemem[0].to_ulong()));
	fBeforelcp3->write(char(Beforelcpmem[1].to_ulong()));
	fBeforelcp3->write(char(Beforelcpmem[0].to_ulong()));
    }
    fBefore3->write(fBefore2->read_reverse());
    fBeforelcp3->write(fBeforelcp2->read_reverse());
    delete(fBefore3);
    delete(fBefore2);
    delete(fBeforelcp3);
    delete(fBeforelcp2);
    StackDispose(&lifo_lcp);
    StackDispose(&lifo_mem);
    StackDispose(&lifo_rem);
    StackDispose(&lifo_interval);
    StackDispose(&lifo_intervalmem);
    StackDispose(&lifo_intervalrem);
    remove ("streambis_stacklcp.txt");
    remove ("streambis_stackmem.txt");
    remove ("streambis_stackrem.txt");
    remove ("streambis_intervallcp.txt");
    remove ("streambis_intervalmem.txt");
    remove ("streambis_intervalrem.txt");
    delete [] Beforemem;
    delete [] Beforelcpmem;
    delete [] LCPmem;
    
    return ( 1 );
}


unsigned int GetMaws(  unsigned char * seq_id, stream_reader<uint40>  * fSA, INT n, int sigma, stream_reader<uint40>  * fLCP, char * Before_fname, char * Beforelcp_fname, unsigned int k, unsigned int K, char * out_file,char * out_file_compressed, int r,int f, long ram_use )
{
    FILE * out_fd;
    INT width_to_print=floor(log10(n))+2;
    INT lcp = 0;
    INT mem;
    
    fLCP->goto_set();
    INT * LCPmem= new INT [2];
    
    bool bmem=false;
    stream_writer<int> * fmem= new stream_writer<int> ("streambis_mem.txt", ram_use/2);
    TStackinfile lifo_lcp;
    StackNew(&lifo_lcp, sizeof(INT), "streambis_stack.txt", ram_use/2);
    StackPush(&lifo_lcp,&lcp);

 //we compute a bitvector that contains a `1', if an identical row has already been seen => to avoid duplicates.
    long ftell_mem;
    for ( INT i = 0; i < n; i++ )
    {
        LCPmem[0]=fLCP->read();
	StackPop(&lifo_lcp,&lcp);
        while(!StackEmpty(&lifo_lcp)&&lcp>LCPmem[0])
        {
	    StackPop(&lifo_lcp,&mem);
            if ( mem == LCPmem[0] )
            {
                bmem=true;
                fmem->write(1);
            }
            lcp = mem;
        }
	StackPush(&lifo_lcp,&lcp);
	lcp=LCPmem[0];
	StackPush(&lifo_lcp,&lcp);
        if (!bmem)
	    fmem->write(0);
        else
            bmem=false;
    }
    
    StackDispose(&lifo_lcp);
   delete(fmem);
   remove ("streambis_stack.txt");
   stream_reader<int> * fmemr = new stream_reader<int>("streambis_mem.txt", ram_use/4);
   stream_reader<char> * fBefore = new stream_reader<char> (Before_fname, ram_use/4);
   stream_reader<char> * fBeforelcp = new stream_reader <char> (Beforelcp_fname, ram_use/4);
   stream_writer<triplet> * fout= new stream_writer<triplet> (out_file_compressed, ram_use/4,1);
   if (fout==NULL)
	return 1;
   fBefore->goto_end(2*n+1);
   fBeforelcp->goto_end(2*n+1); 
	/* Print a separator */
    triplet T;
    changetriplet(&T,'>',0,0);
    fout->write(T);
    INT offset_out=fout->getpos();
    bitset<8> * Beforemem = new bitset<8> [2];
    if( ( Beforemem == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for Before.\n" );
        return ( 0 );
    }
    bitset<8>* Beforelcpmem = new bitset<8>[2];
    if( ( Beforelcpmem == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory for BeforeLCP.\n" );
        return ( 0 );
    }
    fSA->goto_set();
    fLCP->goto_set();
    fmemr->goto_set();
    INT SAmem;
    
    // la ligne -1
    Beforemem[0]=bitset<8>(fBefore->read_reverse());
    Beforelcpmem[0]=bitset<8>(fBeforelcp->read_reverse());
    LCPmem[1]=fLCP->read();
    mem=fmemr->read();
    INT nbmaw=0;
    for ( INT i = 0; i < n; i++ )
    {
   	    Beforemem[0]=bitset<8>(fBefore->read_reverse());
            Beforelcpmem[0]=bitset<8>(fBeforelcp->read_reverse());
            Beforemem[1]=bitset<8>(fBefore->read_reverse());
            Beforelcpmem[1]=bitset<8>(fBeforelcp->read_reverse());
	
            SAmem=fSA->read();
            mem=fmemr->read();
            LCPmem[0]=LCPmem[1];
            if (i<n-1)
            {
               LCPmem[1]=fLCP->read();
            }
            
            
        	for( int l = 0; l < sigma; l++ )
        	{
		bool B1 = (
                           Beforemem[0][sigma-1-l]==0 &&
                           Beforelcpmem[0][sigma-1-l]==1 &&
                           SAmem + LCPmem[0] < n &&
                           (r==0 || SAmem + LCPmem[0] < (n-1)/2 || SAmem >(n-1)/2 ) &&
                           LCPmem[0] + 2 <= K	&&
                           LCPmem[0] + 2 >= k	);
                
                bool B2 = (
                           i < n - 1 &&
                           Beforemem[1][sigma-1-l]==0 &&
                           Beforelcpmem[1][sigma-1-l]==1  &&
                           SAmem + LCPmem[1] < n &&
                           (r==0 || SAmem + LCPmem[1] < (n-1)/2 || SAmem >(n-1)/2  ) &&
                           LCPmem[1] + 2 <= K &&
                           LCPmem[1] + 2 >= k &&
                           mem == 0 );


		    	if ( B1 )
		    	{
                    		uint40 start = SAmem;
                    		uint40 size = SAmem+ LCPmem[0] + 1 - start;
				changetriplet( &T, Mapping(l),(uint40) start,size);
				fout->write(T);	
				nbmaw++;			

		    	}
		    	else if ( B2 )
		    	{
                   		uint40 start = SAmem;
                    		uint40 size = SAmem + LCPmem[1] + 1 - start;
				changetriplet (&T, Mapping(l),(uint40) start,size);
				fout->write(T);	
				nbmaw++;			
		    	}
        	}
    	}
    

   	delete(fout);
	delete(fmemr);
   	delete(fBefore);
   	delete(fBeforelcp);
    	delete [] Beforemem;
    	delete [] Beforelcpmem;
   	delete [] LCPmem;
 	remove("streambis_mem.txt");

	cout<< nbmaw << " minimal absent words have been found."<<endl;

	if (f==1)
	{
	   if ( ! ( out_fd = fopen ( out_file, "a") ) )
	   {
		fprintf ( stderr, " Error: Cannot open file %s!\n", out_file );
		return ( 1 );
	   }
	    fprintf(out_fd,">%s\n",(char *) seq_id); 
	    stream_reader<triplet> * fout_r= new stream_reader<triplet> (out_file_compressed, ram_use/2);
	    fout_r->goto_pos(offset_out);
	    stream_reader<unsigned char> * fseq= new stream_reader<unsigned char> ("seq.txt",ram_use/2);
            unsigned char c=' ';
	    char * maw= new char[K+1];
	    INT j=-1;
	    while (! fout_r->empty())
	    {	
		j++;
	        T=fout_r->read();
		maw[0]=T.c;
	        for (INT i=0; i<T.size; i++)
	        {
		    c=fseq->getValue(T.start+i);
	 	    maw[1+i]=c;
	        }
		maw[T.size+1]='\0';
	        fprintf(out_fd,"%s\n",maw);
	    }
	    cout <<"maw printed"<<endl;
	    delete [] maw;
	    delete(fout_r);
	    remove(out_file_compressed);
	    delete(fseq);
	    fprintf( out_fd, "\n" );
	    if ( fclose ( out_fd ) )
	    {
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	    }
        }
	return ( 1 );
}
