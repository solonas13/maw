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
#include <iostream>
#include <sstream>  
#include <tgmath.h> 
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

#include <list>

#include <omp.h>

#include <sys/resource.h>
#include <unistd.h>

#include "mawdefs.h"
#include "lcp.h"

#include <divsufsort.h>                                         // include header for suffix sort
#include <divsufsort64.h>                                       // include header for suffix sort
#include <sdsl/bit_vectors.hpp>					// include header for bit vectors

#undef _LOG_WRITE
//#undef _LOG_OUT

using namespace sdsl;
using namespace std;


Maw::Maw(char letter, int start, int length)
{ this->letter=letter; this->start=start; this->length=length;}

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

INT compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw )
{
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT n = strlen ( ( char * ) seq );
	INT sigma=0;


    if      ( ! strcmp ( "DNA", sw . alphabet ) )   sigma = strlen ( ( char * ) DNA );
    else if ( ! strcmp ( "PROT", sw . alphabet ) )  sigma = strlen ( ( char * ) PROT );
#ifdef _LOG_WRITE
    std::cout<< "sigma "<< sigma<<std::endl;
    std::cout<< sw.alphabet<<std::endl;
#endif
    /* Compute the suffix array */
#ifdef _LOG_WRITE
    cout << "Compute SA"<<endl;
#endif

    SA = ( INT * ) malloc( (size_t ) ( n ) * sizeof(INT));
    if( ( SA == NULL) )
    {
            fprintf(stderr, " Error: Cannot allocate memory.\n" );
            return ( 0 );
    }
#ifdef _USE_64
    if( divsufsort64( seq, SA, n  ) != 0 )
    {
            fprintf(stderr, " Error: Cannot allocate memory.\n" );
            exit( EXIT_FAILURE );
    }
#endif
#ifdef _USE_32
    if( divsufsort( seq, SA, n  ) != 0 )
    {
            fprintf(stderr, " Error: Cannot allocate memory.\n" );
            exit( EXIT_FAILURE );
    }
#endif

#ifdef _LOG_WRITE
cout << "Compute invSA"<<endl;
#endif
    /*Compute the inverse SA array */
    invSA = (INT *) calloc( ( size_t )  ( n ) , sizeof( INT ) );
    if( ( invSA == NULL) )
    {
            fprintf(stderr, " Error: Cannot allocate memory.\n" );
            return ( 0 );
    }
    #pragma omp parallel for
    for ( INT i = 0; i < n; i ++ )
    {
            invSA [SA[i]] = i;
    }
#ifdef _LOG_WRITE
    cout << "Compute LCP"<<endl;
#endif
    /* Compute the LCP array */
    LCP = LCParray( seq, n, SA, invSA );
#ifdef _LOG_WRITE
    cout << "Free invSA"<<endl;
#endif
    free ( invSA );

    bit_vector * max_loc= new bit_vector[1];
    max_loc[0]=bit_vector(n,0);
#ifdef _LOG_WRITE
    cout << "Compute max"<< endl;
#endif
    GetMax(n,  sigma, LCP, max_loc);

    rrr_vector<> rrrmax(max_loc[0]);
    rrr_vector<>::select_1_type max_sel(&rrrmax);
    rrr_vector<>::rank_1_type max_rank(&rrrmax);

    INT nb_max=rrrmax.size();
    delete[] max_loc;

#ifdef _LOG_WRITE
    cout << "Compute Maws"<<endl;    
#endif

    Get_min_abs_w( seq, seq_id,SA, n , sigma, LCP,max_sel, max_rank, nb_max,sw.k, sw.K, sw.output_maw_filename, sw . r, sw . threads );

    free ( SA );
    free ( LCP );
    sdsl::util::clear(rrrmax);
    sdsl::util::clear(max_sel);
    sdsl::util::clear(max_rank);

    return ( 1 );
}


unsigned char Mapping( INT a )
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
                c = 'N';
                break;
            case 4:
                c = 'T';
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

INT RevMapping ( unsigned char b )
{
	INT a = -1;
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
            case 'N':
                a = 3;
                break;
            case 'T':
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
INT RevComStr ( unsigned char * str, unsigned char * str2, INT iLen )
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


// store in a bit_vector the position of each local maximum
INT GetMax(INT n, INT sigma, INT * LCP, bit_vector * max_loc){
    TStack lifo_lcp;
    StackNew (&lifo_lcp, sizeof(INT));
    INT lcp;
    bool up;
    up=false;
    lcp = 0;
    StackPush(&lifo_lcp, &lcp);

    for (INT i =0; i<n+1; i++)
    {
        StackTop(&lifo_lcp,&lcp);
        if (lcp>LCP[i] && up)
        {
            max_loc[0][i-1]=1;
            up=false;
        }
        if (lcp<LCP[i]&& !up){up=true;}

        while(!StackEmpty(&lifo_lcp)&&lcp>LCP[i])
        {
            StackPop(&lifo_lcp,&lcp);
        }
        lcp=LCP[i];
        StackPush(&lifo_lcp,&lcp);
	}
    StackDispose(&lifo_lcp);
    return (0);
}

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
        INT threads )
{

    int nb_threads= threads;

    omp_set_num_threads( nb_threads ); 
#ifdef _LOG_WRITE
    printf("There are %d threads\n", nb_threads);
#endif    
    int l=(int) (log(n)/log(sigma-1));
    if (l>32){l=32;}
    INT lcpmax=0;

    INT c=r;
    int T=l;
    INT * Count= ( INT * ) calloc( (size_t) l, sizeof( INT ) );
    if( ( Count == NULL) )
    {
      fprintf(stderr, " Error: Cannot allocate memory.\n" );
      return ( 0 );
    }
    
    for (INT i=0; i<n; i++)
    {
        if (LCP[i]<l && seq[SA[i]+LCP[i]]!='$'){Count[LCP[i]]+=1;}
	    if (LCP[i]>lcpmax){lcpmax=LCP[i];}
    }

#ifdef _LOG_WRITE
    cout  << "max of LCP value :"<< lcpmax <<endl;
#endif    
    INT V2=Count[0];

#ifdef _LOG_WRITE
    cout << "Count[0] " << Count[0]<< endl;
#endif
    INT V=V2-1;
    for (INT i=0; i<l; i++)
    {
        if (Count[i]<V+1+r){T=i;break;}
        V=V*V2;
    }
    free ( Count);
#ifdef _LOG_WRITE
    cout << "lower bound of MAW length :"<< T+1 << endl;
#endif    
    if (T==0 || T==1)
    {
        if (k<2)   
            k=2;
        T=0; V=V2-1;
    }
    
    else{T=T-1; V/V2;}
    if (T>0 && k<T+2){k=T+2;}  
    if (K>lcpmax+2){K=lcpmax+2;}
    
#ifdef _LOG_WRITE
    cout << "k="<<k<<" K="<<K<<endl;
#endif    
    while (V>=V2*nb_threads && V>=sqrt(n))
    {
        T=T-1; V=V/V2;
    }
#ifdef _LOG_WRITE
    cout << "T "<<T<<"   V "<<V<<endl;
#endif
    INT * Cut=( INT * ) calloc( (size_t)(2*V+3+r), sizeof( INT ) );
    INT y=0;
    Cut[y++]=0;
    if (T>0)
    {
        for (INT i=0; i<n; i++)
        {
            if (LCP[i]==T){Cut[y++]=i;}
        }
    }
    Cut[y]=n+1;

#ifdef _LOG_WRITE
    cout << "number of independent parts "<< y<<endl;
#endif
    char * tmp;
    tmp = ( char * ) calloc( (size_t) 100, sizeof( char ) );

    sprintf( tmp, "%s", fd_maw);
    std::ofstream file_maw;
    //file_maw.open ( tmp );
    file_maw . open ( tmp, std::fstream::out | std::fstream::app );		
    free(tmp);
    file_maw << ">"<< seq_id << "\n";

    list<Maw> *  Table= new list<Maw> [nb_threads];
    double starttime=omp_get_wtime();
    #pragma omp parallel for
    for (INT z=0; z<y; z++)
    {
        TStack lifo_pos;
        TStack lifo_set;
        TStack lifo_att;
        TStack lifo_mem_pos;
        char * maw;
        INT * Setletter;
        INT * Setletter1;
        int id=omp_get_thread_num();

        StackNew (&lifo_pos, sizeof(INT));
        StackNew (&lifo_set, sizeof(INT)*(sigma+2));
        StackNew (&lifo_att, sizeof(INT));
        StackNew (&lifo_mem_pos, sizeof(INT));
        maw = ( char * ) calloc( ( (size_t)(K+1) ) , sizeof( char ) );
        INT x, left;
        INT right;
        INT pos=-1;
        INT pos2=-1;
        INT lx; INT lleft;INT lright;
        INT lcp; INT lcpleft; INT p; INT pos3=-1;

        StackPush(&lifo_pos, &pos);
        pos=Cut[z]; StackPush(&lifo_att, &pos);

        Setletter = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
        if( ( Setletter == NULL) )
        {
          fprintf(stderr, " Error: Cannot allocate memory.\n" );
        }
        
        StackPush(&lifo_set, &*Setletter);
        INT lpos;
        bool Min_right;
        INT leftint;
        INT rightint;
        Setletter1 = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    	if( ( Setletter1 == NULL) )
    	{
      	  fprintf(stderr, " Error: Cannot allocate memory.\n" );
    	}
        INT mini; INT maxi;	
        if (z==0){mini=0;}
        else{mini=max_rank(Cut[z]);}
        maxi=max_rank(Cut[z+1]);

        for (INT t=mini; t<maxi; t++)
        {
            x=max_sel(t+1);
                left=x-1;
                right=x+1;
                for(INT i=0; i<sigma+2 ;i++)
                {
                    Setletter[i]=0;
                }
                Min_right=false;

                while (!Min_right&& x+1<n+1 && x+1>0)
                {
                    if (!StackEmpty(&lifo_att))
                    {
                        StackTop(&lifo_att,&pos2);
                    }

                    if (left==-1){lcpleft=-1;}
                    else{lcpleft=LCP[left];}
                    if (right==n+1){lright=-1;}
                    else{lright=LCP[right];}
                    if (SA[x]>0)
                    {
                        lx=RevMapping(seq[SA[x]-1]);
                    }
                    else{lx=-1;}
                    lcp=LCP[x];
                    Setletter[sigma+1]=lcp;

                    StackTop(&lifo_pos, &pos);

                    if (pos==-1){lpos=-1;}
                    else{lpos=LCP[pos];}

        //1rst step
                    while(lcp<lpos && pos+1>0)
                    {

                        Empty_stack(&lifo_att, &lifo_pos, &lifo_set, &Setletter, sigma,n,k,K,seq, &x, &left, &right, LCP,SA,Table);

                        Setletter[sigma+1]=LCP[x];
                        StackTop(&lifo_pos, &pos);
                        if (pos==-1){lpos=-1;}
                        else{lpos=LCP[pos];}

                        if (!StackEmpty(&lifo_att))
                        {
                            StackPop(&lifo_att,&pos2);
                            if(pos2>Cut[z]&&(lpos<LCP[pos2] || (lpos==LCP[pos2] && pos<pos2)))
                            {
                                p=1;
                                while (pos2-p>=0 && LCP[pos2-p]<=LCP[pos2] && (pos2-p==0 || LCP[pos2-p-1]<=LCP[pos2-p])&&pos!=pos2-p)
                                {
                                    p++;
                                }
                                for (INT i=0; i<sigma+1 ; i++){Setletter1[i]=0;}
                                p=p-1;
                                while (p>0)
                                {
                                    pos3=pos2-p;
                                    p=p-1;
                                    if(LCP[pos3]==LCP[pos3+1])
                                    {
                                        Setletter1[0]=1;
                                    }
                                    else{Setletter1[0]=0;}
                                    Setletter1[sigma+1]=LCP[pos3];
                                    StackPush(&lifo_set, &*Setletter1);
                                    StackPush(&lifo_pos,&pos3);
                                }
                                Setletter1[0]=0;
                                Setletter1[sigma+1]=LCP[pos2];
                                StackPush(&lifo_set, &*Setletter1);
                                StackPush(&lifo_pos,&pos2);
                            }
                           else{if (pos2+1>Cut[z]&&lpos!=LCP[pos2]) {StackPush(&lifo_att,&pos2);}}
                        }

                        StackTop(&lifo_pos, &pos);
                        if (pos==-1){lpos=-1;}
                        else{lpos=LCP[pos];}
                   }

        //2nd step
                if (pos==-1){lpos=-1;}
                else{lpos=LCP[pos];}
                if (left<0){lcpleft=-1;}
                else{lcpleft=LCP[left];}
                if(left >=0 && SA[left]>0){lleft=RevMapping(seq[SA[left]-1]);}
                else{lleft=-1;}
                leftint=left;
                rightint=right;

                if (lcp>lcpleft && lcp>lright && lcp>lpos&& x+1<n+1 )
                {

                    if (x+1==right && lcp+2 <=K && lcp+2>=k)
                    {
                    // longest proper suffix of MAWS : seq[SA[x] ... SA[x]+lcp], MAWS of length lpos+2
                        if (SA[x]+lcp<n)
                        {
                            for (INT i=1; i<sigma+1; i++)
                            {
                                if( (Setletter[i]==1 || lleft==i-1) && i-1!=lx)
                                {
                                    maw[0] = Mapping(i-1);
                                    memcpy( &maw[1], &seq[SA[x]], lcp+1 );
                                    maw[lcp + 2] = '\0';
                                    if ( ! ( strstr ( maw, "$" ) ) )
                                    {
#ifdef _LOG_OUT
                                        Table[id].push_back(Maw(maw[0],SA[x],lcp+1));
#endif
                                    }
                                }
                            }
                        }

                    // longest proper suffix of MAWS : seq[SA[pos-1] ... SA[pos-1]+lpos], MAWS of length lpos+2
                        if (x>0 && SA[x-1]+lcp<n && lx!=lleft && lx>=0 && Setletter[lx+1]==0 )
                        {
                            maw[0] = Mapping(lx);
                            memcpy( &maw[1], &seq[SA[x-1]], lcp+1 );
                            maw[lcp + 2] = '\0';
                            if ( ! ( strstr ( maw, "$" ) ))
                            {
#ifdef _LOG_OUT
                                Table[id].push_back(Maw(maw[0],SA[x-1],lcp+1));
#endif
                            }
                        }
                    }
                    else   // x=left+1
                    {
                    // longest proper suffix of MAWS : seq[SA[x] ... SA[x]+lcp], MAWS of length lpos+2
                        if (SA[x]+lcp<n && lleft !=lx && lleft>=0 && Setletter[lleft+1]==0 && lcp+2 <=K && lcp+2>=k)
                        {
                            maw[0] = Mapping(lleft);
                            memcpy( &maw[1], &seq[SA[x]], lcp+1 );
                            maw[lcp + 2] = '\0';
                            if ( ! ( strstr ( maw, "$" ) ))
                            {
#ifdef _LOG_OUT
                                Table[id].push_back(Maw(maw[0],SA[x],lcp+1));				    
#endif
                            }
                        }
                    // longest proper suffix of MAWS : seq[SA[pos-1] ... SA[pos-1]+lpos], MAWS of length lpos+2
                        if (x>0 && SA[x-1]+lcp<n && lcp+2 <=K && lcp+2>=k)
                        {
                            for (INT i=1; i<sigma+1; i++)
                            {
                                if( (Setletter[i]==1 || lx==i-1) && i-1!=lleft)
                                {
                                    maw[0] = Mapping(i-1);
                                    memcpy( &maw[1], &seq[SA[x-1]], lcp+1 );
                                    maw[lcp + 2] = '\0';
                                    if ( ! ( strstr ( maw, "$" ) ))
                                    {
#ifdef _LOG_OUT
                                        Table[id].push_back(Maw(maw[0],SA[x-1],lcp+1));
#endif
                                    }
                                }
                            }
                        }
                    }

                    if (x>=0 && SA[x]>0 && RevMapping(seq[SA[x]-1])>=0)
                    {
                        Setletter[RevMapping(seq[SA[x]-1])+1]=1;
                    }
                }

    // 3rd step
                if (!StackEmpty(&lifo_pos))
                {
                    StackPop(&lifo_pos,&pos);
                    if (!StackEmpty(&lifo_pos))
                    {
                        StackTop(&lifo_pos,&pos3);
                    }
                    StackPush(&lifo_pos,&pos);
                }
                if (x>=0 && x+1<n+1&& left+1>0  && LCP[x]+2<k && right+1<n+2)
                {
                    if (x>=0 && SA[x]>0 && RevMapping(seq[SA[x]-1])>=0)
                    {
                        Setletter[RevMapping(seq[SA[x]-1])+1]=1;
                    }
                    if (LCP[left]<=LCP[x] && LCP[left]>=LCP[right] && LCP[right]<=LCP[x])
                    {
                        x=left; left=left-1;
                    }
                    else
                    {
                        if (LCP[right]>LCP[x])
                        {
                            Min_right=true;
                        }
                        if (LCP[right]+2>=k)
                        {
                            if(pos+1>0 && LCP[pos]==LCP[x])
                            {
                                StackTop(&lifo_set, &*Setletter1);
                                Setletter[0]=1;
                                if (pos!=x)
                                {
                                    StackPush(&lifo_pos, &x);
                                }
                                if (((pos3+1>0 && LCP[pos3]==LCP[pos])) && Setletter1[0]==0)
                                {
                                    StackPop(&lifo_set, &*Setletter1);
                                }
                                StackPush(&lifo_set, &*Setletter);
                            }
                            else
                            {
                                Setletter[0]=0;
                                StackPush(&lifo_pos, &x);
                                StackPush(&lifo_set, &*Setletter);
                            }
                        }

                        x=right; right=right+1;
                    }
                }
                else
                {
                    Min_right=Next_pos(&lifo_att, sigma,&x, &left, &right,&Setletter,n,k, LCP, &lifo_pos, &lifo_set,SA, seq);
                }
            }
	}
        lcp=-1;
        StackTop(&lifo_pos, &pos);
        if (pos>-1){lpos=LCP[pos];}
        else {lpos=-1;}
        if (mini<maxi && x>=0 && x<n && SA[x]>0 && RevMapping(seq[SA[x]-1])>=0){Setletter[RevMapping(seq[SA[x]-1])+1]=1;}


        while(mini < maxi && lcp < lpos && pos+1 > 0)
        {
            while (left+1>0 && LCP[left]>=lpos)
            {
                StackPop(&lifo_pos, &left);
                StackPush(&lifo_mem_pos,&left);
            }

            if (!StackEmpty(&lifo_att))
            {
                StackTop(&lifo_att,&pos2);
                if (LCP[pos2]<=LCP[pos]&&pos2>left){left=pos2;}
            }

            while (!StackEmpty(&lifo_mem_pos))
            {
                StackPop(&lifo_mem_pos, &pos);
                StackPush(&lifo_pos,&pos);
            }
            Empty_stack(&lifo_att, &lifo_pos, &lifo_set, &Setletter, sigma,n,k,K,seq, &x, &left, &right, LCP,SA,Table);
            StackTop(&lifo_pos, &pos);
            if (pos==-1){lpos=-1;}
            else{lpos=LCP[pos];}

            if (!StackEmpty(&lifo_att))
            {
                StackPop(&lifo_att,&pos2);
                if (pos2+1>Cut[z]&&(lpos<LCP[pos2] || (lpos==LCP[pos2] && pos<pos2)))
                {
                    p=1;
                    while (pos2-p>=0 && LCP[pos2-p]<=LCP[pos2] && (pos2-p==0 || LCP[pos2-p-1]<=LCP[pos2-p])&&pos!=pos2-p)
                    {
                        pos3=pos2-p;
                        p++;
                    }
                    for (INT i=0; i<sigma+1 ; i++){Setletter1[i]=0;}
                    Setletter1[sigma+1]=LCP[pos2];
                    p=p-1;
                    while (p>0)
                    {
                        pos3=pos2-p;
                        p=p-1;
                        if(LCP[pos3]==LCP[pos3+1])
                        {
                            Setletter1[0]=1;
                        }
                        else{Setletter1[0]=0;}
                        Setletter1[sigma+1]=LCP[pos3];
                        StackPush(&lifo_set, &*Setletter1);
                        StackPush(&lifo_pos,&pos3);
                    }
                    Setletter1[0]=0;
                    Setletter1[sigma+1]=LCP[pos2];
                    StackPush(&lifo_set, &*Setletter1);
                    StackPush(&lifo_pos,&pos2);
                }
                else{if (pos2+1>0 && lpos!=LCP[pos2]) {StackPush(&lifo_att,&pos2);}}
            }
        }

        free(Setletter);
        free(Setletter1);
        StackDispose(&lifo_pos);
        StackDispose(&lifo_set);
        StackDispose(&lifo_att);
        StackDispose(&lifo_mem_pos);
        free(maw);

	Setletter=NULL;
	Setletter1=NULL;
	maw=NULL;

    }

    double endtime=omp_get_wtime();
    free (Cut);
    char* mawb = ( char * ) calloc( ( (size_t)(K) ) , sizeof( char ) );
    int nb_maw=0;
    for (int i=0; i< nb_threads; i++)
    {
	while (!Table[i].empty())
	{
	   file_maw<<Table[i].front().letter;
	   memcpy( &mawb[0], &seq[Table[i].front().start], Table[i].front().length );
           mawb[Table[i].front().length ] = '\0';
	   file_maw << mawb<< "\n";
           Table[i].pop_front();
	   nb_maw++;
	}
    }
    delete[] Table;
    file_maw << "\n";
    file_maw.close();
    free(mawb);
    #ifdef _LOG_WRITE
    std::cout << nb_maw << std::endl;
    printf("Computation time %.6g \n",endtime-starttime);
    #endif 
    return ( 0 );
}


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
                        INT * Tx,
                        INT * Tleft,
                        INT * Tright,
                        INT * LCP,
                        INT *SA,
			list<Maw> * Table)
{
    INT * Setletterbis;
    INT * Setletter1;
    INT * Setletter;
    Setletter = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    Setletterbis = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    Setletter1 = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    TStack lifo_mem_set;
    TStack lifo_mem_pos;
    StackNew (&lifo_mem_set, sizeof(INT)*(sigma+2));
    StackNew (&lifo_mem_pos, sizeof(INT));
    INT pos, pos2, lpos,pos3;

    int id=omp_get_thread_num();


    for (INT i=0; i< sigma+2; i++)
    {
        Setletter[i]=TSetletter[0][i];
    }
    INT x, left, right;
    x=Tx[0];
    left=Tleft[0];
    right=Tright[0];
    INT t=0;
    INT leftint;
    INT rightint;
    leftint=left;
    rightint=right;

    char * maw;
    maw = ( char * ) calloc( ((size_t)( K + 1 )) , sizeof( char ) );
    if( ( maw == NULL) )
    {
        fprintf(stderr, " Error: Cannot allocate memory.\n" );
        return ( 0 );
    }

    StackPop(&lifo_pos[0], &pos);
    if (!StackEmpty(&lifo_pos[0]))
    {
        StackTop(&lifo_pos[0],&pos3);
    }
    else{pos3=-1;}
    StackPush(&lifo_pos[0], &pos);
    if (pos+1>0){lpos=LCP[pos];}
    else{lpos=-1;}


    Setletter[0]=0;
    Setletter[sigma+1]=lpos;

    for (INT i=1; i<sigma+1; i++)
    {
        Setletterbis[i]=Setletter[i];
    }

    if (lpos+2>=k){
        StackPop(&lifo_set[0], &*Setletter1);
        if (pos<x){rightint=x;}
        if (pos<right-2)
        { // in this case Setletter is the right extrem set of the considered interval
            if (Setletter1[0]==1 || (pos3<0 ||(pos3>=0 && LCP[pos3]!=LCP[pos]))) // in this case we keep Setletter1 and add Setletter
            {
                StackPush(&lifo_set[0], &*Setletter1);
                StackPush(&lifo_mem_set, &*Setletter);
            }
            else
            {
                t=0;
                for (INT i=0; i< sigma+1; i++)
                {
                    t=t+Setletter1[i];
                }
                if (t>0) // the set has already been added
                {
                    Setletter1[0]=1;
                    StackPush(&lifo_mem_set,&*Setletter1);
                    for (INT i=1; i<sigma+1; i++)
                    {
                        Setletterbis[i]=Setletterbis[i]+Setletter1[i]-Setletterbis[i]*Setletter1[i];
                    }
                }
                else
                {
                    StackPush(&lifo_mem_set,&*Setletter);
                }
            }
        }
        else
        {
            StackPush(&lifo_set[0], &*Setletter1);
            if (Setletter1[0]==1)
            {
                for(INT i=0; i<sigma+1 ;i++)
                {
                    Setletter1[i]=0;
                }
                StackPush(&lifo_mem_set, &*Setletter1);
            }
        }
    }

    StackPop(&lifo_set[0], &*Setletter1);
    StackPop(&lifo_pos[0],&pos);
    while (!StackEmpty(&lifo_pos[0]) && pos+1>0 && LCP[pos]==lpos)
    {
        StackPush(&lifo_mem_pos,&pos);
        if (pos>=0 && pos<n && SA[pos]>0 && RevMapping(seq[SA[pos]-1])>=0)
        {
            Setletterbis[RevMapping(seq[SA[pos]-1])+1]=1;
        }
        StackPop(&lifo_pos[0],&pos);
    }
    StackPush(&lifo_pos[0],&pos);

    while (!StackEmpty(&lifo_set[0]) && Setletter1[sigma+1]==lpos)
    {
        StackPush(&lifo_mem_set,&*Setletter1);
        for (INT i=1; i<sigma+1; i++)
        {
            Setletterbis[i]=Setletterbis[i]+Setletter1[i]-Setletterbis[i]*Setletter1[i];
        }
        StackPop(&lifo_set[0],&*Setletter1);
    }
    StackPush(&lifo_set[0],&*Setletter1);

    if (lpos+2>=k)
    {
        StackTop(&lifo_mem_pos, &pos);
        StackTop(&lifo_mem_set,&*Setletter1);

        if (pos>left+2)
        {// in this case Setletter is the left extrem set of the considered interval
            if (Setletter1[0]==1)
            {
                StackPush(&lifo_mem_set,&*Setletter);
            }
            else
            {
		t=0;
                for (INT i=0; i< sigma+1; i++)
                {
                    t=t+Setletter1[i];
                }
                if (t==0)
                {
                    StackPop(&lifo_mem_set,&*Setletter1);
                    StackPush(&lifo_mem_set,&*Setletter);
                }
                else
                {
                    t=0;
                    if (pos3<0 || (pos3>=0 &&LCP[pos3]!=LCP[pos]))
                    {
                        StackPush(&lifo_mem_set,&*Setletter);
                    }
                }
            }
        }
        else
        {
            if(Setletter1[0]==1)
            {
                for(INT i=0; i<sigma+1 ;i++)
                {
                    Setletter1[i]=0;
                }
                StackPush(&lifo_mem_set,&*Setletter1);
            }
        }
    }
             // this way we have k+1 sets and k positions, so for each postion we have the preceding and the following sets
    StackPop(&lifo_mem_set, &*Setletter1); // this is the first set

    for (INT i=0; i<sigma+1; i++)
    {
        Setletter[i]=Setletterbis[i];
    }
    if (x<=pos-1)
    {
        if (x>=0 && SA[x]>0 && RevMapping(seq[SA[x]-1])>=0)
        {
            Setletter1[RevMapping(seq[SA[x]-1])+1]=1;
            Setletterbis[RevMapping(seq[SA[x]-1])+1]=1;
        }
        leftint=x;
    }
    else
    {
        if(pos+1>0 && left>pos-1){left=pos-1;}
        if (left+1>0 && LCP[left]>lpos)
        {
            StackTop(&lifo_pos[0],&left);
        }

        if (!StackEmpty(&lifo_att[0]))
        {
            StackTop(&lifo_att[0],&pos2);
            if (pos2+1>0 && pos+1>0 && LCP[pos2]<=LCP[pos]&&pos2>left){left=pos2;}
        }

        if (left >=0 && SA[left]>0 && RevMapping(seq[SA[left]-1])>=0)
        {
            Setletter1[RevMapping(seq[SA[left]-1])+1]=1;
            Setletterbis[RevMapping(seq[SA[left]-1])+1]=1;
        }
        leftint=left;
    }


if (lpos+2>=k&& lpos+2<=K){
    while(!StackEmpty(&lifo_mem_pos))
    {
        StackPop(&lifo_mem_pos, &pos);

                // longest proper suffix of MAWS : seq[SA[pos-1] ... SA[pos-1]+lpos], MAWS of length lpos+2
        if (pos>0 && LCP[pos-1]!=LCP[pos] && lpos>0 && SA[pos-1]+lpos<n&& Setletter1[0]==0 &&pos<n)
        {
            for (INT i=1; i<sigma+1; i++)
            {
                if(Setletterbis[i]==1 && Setletter1[i]==0)
                {
                    maw[0] = Mapping(i-1);
                    memcpy( &maw[1], &seq[SA[pos-1]], lpos+1 );
                    maw[lpos + 2] = '\0';
                    if ( ! ( strstr ( maw, "$" ) ))
                    {
#ifdef _LOG_OUT
Table[id].push_back(Maw(maw[0],SA[pos-1],lpos+1));
#endif
                    }
                }
            }
        }

        if (!StackEmpty(&lifo_mem_set)){
        StackPop(&lifo_mem_set, &*Setletter1);
        }

        // longest proper suffix of MAWS : seq[SA[pos] ... SA[pos]+lpos], MAWS of length lpos+2

        if (pos+1 <n+1 && SA[pos]+lpos<n )
        {
            for (INT i=1; i<sigma+1; i++)
            {
                if(Setletterbis[i]==1 && Setletter1[i]==0 && (SA[pos]==0 || (SA[pos]>0 && i-1!=RevMapping(seq[SA[pos]-1])) ) )
                {
                    maw[0] = Mapping(i-1);
                    memcpy( &maw[1], &seq[SA[pos]], lpos+1 );
                    maw[lpos + 2] = '\0';
                    if ( ! ( strstr ( maw, "$" ) ))
                    {
#ifdef _LOG_OUT
Table[id].push_back(Maw(maw[0],SA[pos],lpos+1));
#endif
                    }
                }
            }
        }

        if (pos>=0 &&pos<n && SA[pos]>0 && RevMapping(seq[SA[pos]-1])>=0)
        {
            Setletter1[RevMapping(seq[SA[pos]-1])+1]=1;
        }

        if (right<pos+1){right=pos+1;}

    }
}


    if (x<=left && x>0){left=x-1;}
    if (x>=right && x+1 <n+1){right=x+1;}


    for (INT i=0; i< sigma+2; i++)
    {
        TSetletter[0][i]=Setletter[i];
    }
    Tx[0]=x;
    Tleft[0]=left;
    Tright[0]=right;
    free(Setletter);
    free(Setletter1);
    free(Setletterbis);
    StackDispose(&lifo_mem_pos);
    StackDispose(&lifo_mem_set);
    free(maw);

    return 0;
	}


bool Next_pos(TStack * lifo_att, INT sigma, INT * Tx,INT * Tleft, INT * Tright,INT ** TSetletter, INT n, INT k, INT * LCP, TStack * lifo_pos, TStack * lifo_set, INT * SA, unsigned char * seq)
{

    INT * Setletterbis;
    INT * Setletter1;
    INT * Setletter;
    INT * Setletterg;
    Setletterg = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    Setletter = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    Setletterbis = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    Setletter1 = ( INT * ) calloc ( ( sigma+2 ), sizeof( INT ) );
    TStack lifo_mem_set;
    TStack lifo_mem_pos;
    StackNew (&lifo_mem_set, sizeof(INT)*(sigma+2));
    StackNew (&lifo_mem_pos, sizeof(INT));

    //bool print=0;
    bool Min_right=false;

    for (INT i=0; i< sigma+2;i++)
    {
        Setletter[i]=TSetletter[0][i];
        Setletterg[i]=Setletter[i];
    }
    INT x; INT left; INT right;

    x=Tx[0];
    left=Tleft[0];
    right=Tright[0];

    Setletterbis[sigma+1]=LCP[x];
    Setletter[sigma+1]=LCP[x];
    Setletterg[sigma+1]=LCP[x];

    INT rightc; INT leftb; INT xc; INT xb;
    INT pos; INT pos2; INT pos3;
    bool B_eq=false;

    StackPop(&lifo_pos[0],&pos);
    if (!StackEmpty(&lifo_pos[0]))
    {
        StackTop(&lifo_pos[0],&pos3);
    }
    else{pos3=-1;}
    StackPush(&lifo_pos[0],&pos);

/* if there is an equality ,LCP[x]=LCP[left] or LCP[x]=LCP[right], we need to keep in memory these positions */
    if ((left+1>0 && LCP[x]==LCP[left] && x>pos)|| (right<n+1 && LCP[x]==LCP[right]))
    {
        leftb=left; rightc=right; xc=x; xb=x;
        for (INT i=1; i<sigma+1; i++)
        {
            Setletterbis[i]=0;
        }
        Setletterbis[0]=1;

        if (left+1>0 && LCP[x]==LCP[left] && x>pos) //case 1
        {
            B_eq=true;
            if (right==x+1){Setletterg[0]=1;}
            else {Setletterg[0]=0;}
            StackPush(&lifo_mem_pos, &x);
            StackPush(&lifo_mem_set, &*Setletterg);

            if (x>=0 && SA[x]>0 && RevMapping(seq[SA[x]-1])>=0)
            {
                Setletterg[RevMapping(seq[SA[x]-1])+1]=1;
            }
            xb=left;
            leftb=left-1;

	    while(leftb+1>0 && xb>pos && LCP[xb]==LCP[leftb])
            {
                StackPush(&lifo_mem_pos, &xb);
                StackPush(&lifo_mem_set, &*Setletterbis);

                if (xb>=0 && SA[xb]>0 && RevMapping(seq[SA[xb]-1])>=0)
                {
                    Setletterg[RevMapping(seq[SA[xb]-1])+1]=1;
                }
                xb=leftb;
                leftb=leftb-1;
            }

            pos2=xb;
            if (right==x+1 && (pos==-1 || (pos+1>0 && LCP[xb]!=LCP[pos]))){Setletterbis[0]=0;}
            else {Setletterbis[0]=1;}
			if (xb>=0){
            StackPush(&lifo_mem_set, &*Setletterbis);}

            StackPop(&lifo_mem_set, &*Setletterbis);

            if (pos+1>0 && LCP[pos]==LCP[pos2])
            {
                if (pos!=pos2) // we don't want redundancy in lifo_pos
                {
                    StackPush(&lifo_pos[0], &pos2);
                }
                StackPop(&lifo_set[0],&*Setletter1);
		if (Setletter1[0]==1 ||(Setletter1[0]==0&& ((pos3+1>0 && LCP[pos]!=LCP[pos3])|| pos3<0)))
		{
			StackPush(&lifo_set[0], &*Setletter1);
		}
            }
            else
            {
                StackPush(&lifo_pos[0],&pos2);
                StackPush(&lifo_set[0], &*Setletterbis);
            }

            while (!StackEmpty(&lifo_mem_pos))
            {
                StackPop(&lifo_mem_pos,&pos2);
                StackPush(&lifo_pos[0],&pos2);
            }

            while (!StackEmpty(&lifo_mem_set))
            {
                StackPop(&lifo_mem_set, &*Setletterbis);
                StackPush(&lifo_set[0], &*Setletterbis);
            }
        }

        if (right<n+1 && LCP[x]==LCP[right]) // case 2
        {
            StackPop(&lifo_pos[0], &pos);
            if (!StackEmpty(&lifo_pos[0]))
			{
				StackTop(&lifo_pos[0],&pos3);
			}
            else{pos3=-1;}
            StackPush(&lifo_pos[0], &pos);
            if (x==right-1 && (pos==-1|| (pos+1>0&&LCP[x]!=LCP[pos])))
            {
                Setletter[0]=0;
            }
            else{Setletter[0]=1;}

            if (pos+1>0 && LCP[x]==LCP[pos])
            {
                if (x!=pos){StackPush(&lifo_pos[0],&x);}

                StackPop(&lifo_set[0], &*Setletter1);

				if (Setletter1[0]==1 ||(Setletter1[0]==0&& ((pos3+1>0 && LCP[pos]!=LCP[pos3])|| pos3<0)))
				{
					StackPush(&lifo_set[0], &*Setletter1);
				}
				if (Setletter[0]==1 && !(x==pos && left==pos3 && B_eq && Setletter1[0]==1 ))
				{
					StackPush(&lifo_set[0], &*Setletter);
				}
            }
            else
            {
                StackPush(&lifo_pos[0],&x);
                StackPush(&lifo_set[0],&*Setletter);
            }

            if (x>=0 && SA[x] >0 && RevMapping(seq[SA[x]-1])>=0)
            {
                Setletter[RevMapping(seq[SA[x]-1])+1]=1;
            }

            xc=right;
            rightc=right+1;

            for (INT i=1; i<sigma+1; i++)
            {
                Setletterbis[i]=0;
            }
            Setletterbis[0]=1;

            if (xc>pos)
            {
				while(rightc<n+1 && LCP[xc]==LCP[rightc])
				{
					StackPush(&lifo_pos[0],&xc);
					StackPush(&lifo_set[0],&*Setletterbis);

					if (xc>=0 && SA[xc]>0 && RevMapping(seq[SA[xc]-1])>=0)
					{
						Setletter[RevMapping(seq[SA[xc]-1])+1]=1;
					}
					xc=rightc;
					rightc=rightc+1;
                }

                StackPush(&lifo_pos[0],&xc);
                if (x<right-1)
                {
                    Setletterbis[0]=0;
                }
                else{Setletterbis[0]=1;}

                StackPush(&lifo_set[0],&*Setletterbis);
            }
        }

        for (INT i=0; i<sigma+1; i++)
        {
            Setletter[i]=Setletter[i]+Setletterg[i]-Setletter[i]*Setletterg[i];
        }
        left=leftb;
        if (rightc<n+1 &&  LCP[rightc]>LCP[xc])   // case 2b
        {
            Min_right=true;
            if(((left>0 && LCP[left-1]<=LCP[left]&& LCP[left]<LCP[rightc] )|| (left==0 && LCP[left]<LCP[rightc]) )&& LCP[left]<LCP[x]) 
            {
                //case 2c
                StackTop(&lifo_pos[0],&pos);
                if (!(pos+1>0 && LCP[pos]==LCP[left]&&left<=pos)){
                StackPush(&lifo_att[0],&left);
                }
            }
        }
        left=leftb;
        right=rightc;
        x=xc;
    }
    else
    {
        if (((right<n+1 && left+1>0 && LCP[left]>=LCP[right] )|| (right==n+1) )&&left+1>0 && LCP[left]<=LCP[x])
        { // case 3
            x=left;
            left=left-1;
        }
        else
        {
            if (right<n+1 && LCP[right]>LCP[x])
            { // case 5
                Min_right=true;

                if(((left>0 && LCP[left-1]<=LCP[left]&& LCP[left]<LCP[right] )|| (left==0 && LCP[left]<LCP[right]) )&& LCP[left]<LCP[x] && !B_eq)
                { // case 5b
                    if (!(pos+1>0 && LCP[pos]==LCP[left]&&left<=pos))
                    {
                        StackPush(&lifo_att[0],&left);
                    }
                }
                Setletter[sigma+1]=LCP[x];

                if(pos+1>0 && LCP[pos]==LCP[x])
                {
                    StackTop(&lifo_set[0], &*Setletter1);
                    Setletter[0]=1;
                    if (pos!=x)
                    {
                        StackPush(&lifo_pos[0], &x);
                    }
                    if (((pos3+1>0 && LCP[pos3]==LCP[pos])) && Setletter1[0]==0)
                    {
                        StackPop(&lifo_set[0], &*Setletter1);
                    }
                    StackPush(&lifo_set[0], &*Setletter);
                }
                else
                {
                    Setletter[0]=0;
                    StackPush(&lifo_pos[0], &x);
                    StackPush(&lifo_set[0], &*Setletter);
                }
            }
            if (!((pos==-1 ||(pos+1>0&&LCP[pos]==LCP[x])) && right==n+1))
            { // case 4
                x=right;
                right=right+1;
            }
        }

    }


    if ((left==-1 ||(left+1>0 && LCP[left]>LCP[x]))&& right==n+1 )
    {
        left=pos;
        Min_right=true;
    }

    for (INT i=0; i< sigma+2; i++)
    {
        TSetletter[0][i]=Setletter[i];
    }
    Tx[0]=x;
    Tleft[0]=left;
    Tright[0]=right;
    free(Setletter);
    free(Setletter1);
    free(Setletterbis);
    free(Setletterg);
    StackDispose(&lifo_mem_pos);
    StackDispose(&lifo_mem_set);
    return Min_right;
	}



