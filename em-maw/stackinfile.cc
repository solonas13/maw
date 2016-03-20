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

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>
#include <iostream>
#include "stackinfile.h"
#include<algorithm>
#include<string>
#include<bitset>


void StackNew( TStackinfile *s,INT elemSize, std::string fname, INT bufsize)
{
	assert(elemSize>0);
	s->elemSize = elemSize;
	s->file=std::fopen(fname.c_str(), "w+");
	s->bufelems=bufsize/elemSize;
	if (s->bufelems%2==1) s->bufelems-=1;
	s->filled=0;
	s->offset=0;
	s->buffer=malloc(s->bufelems*elemSize);
}

void StackDispose ( TStackinfile *s)
{
	free( s->buffer);
	std::fclose(s->file);
}

bool StackEmpty(const TStackinfile *s)
{
	return s->offset==0 && s->filled==0;
}


void StackPush( TStackinfile *s,const void *elemAddr)
{
	void *destAddr = (char *)s->buffer + s->filled *s->elemSize;
	memcpy(destAddr, elemAddr, s->elemSize);
	s->filled++;
	if (s->filled==s->bufelems)
		flush(s);
}

void StackPop( TStackinfile *s, void *elemAddr)
{
	const void *sourceAddr;	
	assert(!StackEmpty(s));
	s->filled--;
	sourceAddr = (const char *)s->buffer + s->filled *s->elemSize;
	memcpy(elemAddr,sourceAddr,s->elemSize);
	if (s->offset>0 && s->filled==0)
		refill(s);
}

void flush(TStackinfile *s){
	if (s->filled==s->bufelems)
	{
		INT size=s->bufelems/2;
		std::fwrite(s->buffer, s->elemSize, size, s->file);
		s->filled=size;
		s->offset+=size;
		void* destAddr;
		const void* elemAddr;
		memcpy((char *)s->buffer, (const char*)s->buffer+size*s->elemSize,size*s->elemSize);
	}
	else
	{  
		std::fwrite(s->buffer, s->elemSize, s->filled, s->file);
		s->offset+=s->filled;
		s->filled=0;
	}
}

void refill(TStackinfile *s){
	INT size=s->bufelems/2;
	if (s->offset < size)
	{
		std::fseek(s->file, 0, SEEK_SET);
		s->filled = (INT)std::fread(s->buffer, s->elemSize, s->offset, s->file);
		std::fseek(s->file, 0, SEEK_SET);
	}
	else
	{
		std::fseek(s->file, -size*s->elemSize, SEEK_CUR);
		s->filled = (INT)std::fread(s->buffer, s->elemSize, size, s->file);
		std::fseek(s->file, -size*s->elemSize, SEEK_CUR);
	}
	s->offset-=s->filled;
}
