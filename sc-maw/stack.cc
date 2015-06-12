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
#include "stack.h"

#define kInitialAllocationSize 4

void StackNew( TStack *s,int elemSize)
{
	assert(elemSize>0);
	s->elemSize = elemSize;
	s->logLength = 0;
	s->allocLength = kInitialAllocationSize;
	s->elems = malloc( kInitialAllocationSize * elemSize );
	assert( s->elems != NULL );
}

void StackDispose ( TStack *s)
{
	free(s->elems);
}

bool StackEmpty(const TStack *s)
{
	return s->logLength ==0;
}


void StackPush( TStack *s,const void *elemAddr)
{
	void *destAddr;
	if (s->logLength == s->allocLength){
		s->allocLength *= 2;
		s->elems = realloc(s->elems, s->allocLength *s->elemSize);
		assert(s->elems != NULL);
	}
	destAddr = (char *)s->elems + s->logLength *s->elemSize;
	memcpy(destAddr, elemAddr,s->elemSize);
	s->logLength++;
}

void StackPop( TStack *s, void *elemAddr)
{
	const void *sourceAddr;	
	assert(!StackEmpty(s));
	s->logLength--;
	sourceAddr = (const char *)s->elems + s->logLength *s->elemSize;
	memcpy(elemAddr,sourceAddr,s->elemSize);
}
