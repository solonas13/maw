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

struct TStack
 {
	void *elems;
	int elemSize;
	int logLength;
	int allocLength;
 };

void StackNew( TStack *s,int elemSize );
void StackDispose ( TStack *s );
bool StackEmpty ( const TStack *s );
void StackPush ( TStack *s, const void *elemAddr );
void StackPop( TStack *s,void *elemAddr );
