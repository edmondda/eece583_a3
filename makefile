# This makefile is written for gcc running under Solaris on a SPARCstation.
# To compile on other systems, you may have to change:
# (1) CC to the name of your C compiler.
# (2) LIB_DIR to point at your directory of X11 libraries (libX11.a, etc.)
# (3) On many systems, the -R/usr/openwin/lib option on the LIB line
#     will have to be removed.
# (4) X11_INCLUDE to point at the directory containing "x11/xlib.h" etc.
# (5) FLAGS should be changed to whatever options turn on maximum optimization
#     in your compiler.

CC = gcc

# Need X11 and math libraries.  -R option is sometimes needed under Solaris to
# tell the linker where shared object libraries are.

#LIB = -lX11 -lm -R/usr/openwin/lib
#LIB = -lX11 -lm
LIB = -lX11 -lXau -lws2_32

# Path to the X libraries.

#LIB_PATH = -LC:/cygwin64/lib

#X11_INCLUDE = -IC:/cygwin64/usr/include/

# Flags for lots of warnings, debugging, and full optimization,
# respectively are listed below.

#FLAGS = -Wall -Wpointer-arith -Wcast-qual -Wstrict-prototypes -O -D__USE_FIXED_PROTOTYPES__ -ansi -pedantic -Wmissing-prototypes -Wshadow -Wcast-align -D_POSIX_SOURCE
#FLAGS = -g
FLAGS = -O2

EXE = assignment3

OBJ = assignment3.o graphics.o

SRC = assignment3.c graphics.c

H = graphics.h


$(EXE): $(OBJ)
	$(CC) $(FLAGS) $(OBJ) -o $(EXE) $(LIB_PATH) $(LIB)

graphics.o: graphics.c $(H)
	$(CC) -c $(FLAGS) $(X11_INCLUDE) graphics.c

ass3.o: assignment3.c $(H)
	$(CC) -c $(FLAGS) assignment3.c
