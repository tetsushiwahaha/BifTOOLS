#
# TOOLS installed path. Change it appropriately. 
#
TOOLSDIR = $(HOME)/Bifurcation/TOOLS
INCLUDEFLAGS = -I. -I.. -I$(HOME)/include -I$(TOOLSDIR)/include 
#CFLAGS = -g -DEXMAT -DDEBUG  -DSHOWDELTA $(INCLUDEFLAGS) 
#CFLAGS = -g -DDEBUG -DSHOWDELTA $(INCLUDEFLAGS) 
#CFLAGS = -g -DDEBUG -DSHOWDELTA -DEXMAT -DDEBUG_H -DDEBUG_SYSVAR $(INCLUDEFLAGS) 
CFLAGS = -O -DSHOWDELTA 
#
# -O2 option of clang on FreeBSD generates inappropriate binaries apparently.
#

#
# Some platform has only gcc.
# CC = gcc
#
CC = clang
LIBS  = -lm
#
# Change the binary name if any.
#
OBJ = nn8

$(OBJ): BifMain.o BifTools.o BifInit.o BifFunc.o BifNewton.o
	$(CC) $(CFLAGS) $(LIBFLAGS) $(INCLUDEFLAGS) -o $(OBJ) BifMain.o BifTools.o \
	BifInit.o BifFunc.o $ BifNewton.o $(LIBS)

.c.o:	
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<

BifMain.o: $(TOOLSDIR)/BifMain.c
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<
BifTools.o: $(TOOLSDIR)/BifTools.c
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<
BifInit.o:$(TOOLSDIR)/BifInit.c
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<
BifNewton.o:$(TOOLSDIR)/BifNewton.c
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<
BifFunc.o:./BifFunc.c
	$(CC) $(CFLAGS) -c $(INCLUDEFLAGS) $<

v: 
	vi $(OBJ).c

clean:
	\rm Tools.o *.o *.bak  $(OBJ)

