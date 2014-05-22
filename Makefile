MF=     Makefile
 
CC=     g++
 
CFLAGS= -g -msse3 -O3 -fomit-frame-pointer -funroll-loops 
 
LFLAGS= -I $(HOME)/sdsl/include -L $(HOME)/sdsl/lib -lsdsl -ldivsufsort -ldivsufsort64 -lm -lz
 
EXE=    maw
 
SRC=    maw.cc input.cc stack.cc functions.cc lcp.cc
 
HD=     mawdefs.h stack.h lcp.h Makefile
 
# 
# No need to edit below this line 
# 
 
.SUFFIXES: 
.SUFFIXES: .cc .o 
 
OBJ=    $(SRC:.cc=.o) 
 
.cc.o: 
	$(CC) $(CFLAGS)-c $(LFLAGS) $< 
 
all:    $(EXE) 
 
$(EXE): $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS) 
 
$(OBJ): $(MF) $(HD) 
 
clean: 
	rm -f $(OBJ) $(EXE) *~
