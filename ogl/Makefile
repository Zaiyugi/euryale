CXX = g++

LIB = libOGL.a

OFILES = Matrix.o Vector.o Utility.o Camera.o
AFILES = $(OFILES)

INCLUDES := -I ./ \

LINKS := \
	-L/usr/local/lib \

#CFLAGS = -g
UNAME	:= $(shell uname)
CFLAGS = -g -lm -lGL -lglut -lGLU -Wall -std=c++0x
OFLAGS = -c $(INCLUDES)

ifeq ($(UNAME), Darwin)
CFLAGS = -g -lm -framework OpenGL -framework GLUT
endif

all: $(OFILES) 
	ar rv $(LIB) $?;
	cp $(LIB) ../lib

.C.o:
	${CXX} $(CFLAGS) $< $(OFLAGS)

debug:
	make 'DFLAGS = /usr/lib/debug/malloc.o'

clean:
	-rm *.o *~ ${PROJECT}

