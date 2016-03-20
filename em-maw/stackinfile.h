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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include<algorithm>
#include<string>
#include<bitset>

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif
struct TStackinfile
 {
	INT elemSize;
	std::FILE * file;
	INT bufelems;
	INT filled;
	INT offset;
	void * buffer;
 };

void StackNew( TStackinfile *s,INT elemSize, std::string fname, INT bufsize=(4<<20) );
void StackDispose ( TStackinfile *s );
bool StackEmpty ( const TStackinfile *s );
void StackPush ( TStackinfile *s, const void *elemAddr );
void StackPop( TStackinfile *s,void *elemAddr );

void flush(TStackinfile *s);
void refill(TStackinfile *s);
