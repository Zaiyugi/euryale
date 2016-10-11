
.PHONY: all clean archive gl math
.SUFFIXES: .cpp .C .h .o .cu .cuh

DATE := $(shell date +%F)
UNAME := $(shell uname)

ROOTDIR := .
SRCDIR = src
INCLUDE = include
OBJDIR := obj
LOCAL_LIBDIR = ./lib

HT := /usr
HDSO := /usr/lib64

VPATH = src:include

OFILES = \
	Utility.o Matrix.o LinearAlgebra.o \
	Image.o OIIOFiles.o SPHForce.o SPHEngine.o \

OBJS = $(patsubst %, $(OBJDIR)/%, $(OFILES))

GL_DIR = $(ROOTDIR)/ogl

VR_LIB = $(LOCAL_LIBDIR)/libVR.a

DPA_INCLUDE_DIR := /group/dpa/include

LIBRARY_DIRS := \
	-L$(LOCAL_LIBDIR) \
	-L/group/dpa/lib \
	-L/usr/lib \

LIBRARIES := \
	-lVR \
	-lOGL \
	-lOpenImageIO \

INCLUDES := \
	-I. \
	-I./include \
	-I$(DPA_INCLUDE_DIR) \

GL_INCLUDES := \
	-I. \
	-I$(GL_DIR) \

CXX = g++

CFLAGS = -Wall -g -O3 -fopenmp -std=c++0x
COMPILE_FLAGS = -Wall -g -O3 -fopenmp -std=c++0x
OFLAGS = -c $(INCLUDES)
GLFLAGS = -lGL -lglut -lGLU -lGLEW

all: $(OBJDIR) $(OBJS)
	ar rv $(VR_LIB) $(OBJS)
	make gl

$(OBJDIR):
	@if [ ! -d $(OBJDIR) ]; then \
		echo "-----------------------------"; \
		echo "ERROR: Object directory does not exist"; \
		echo "Creating directory: $(OBJDIR)"; \
		echo "-----------------------------"; \
		mkdir $(OBJDIR); \
	fi

$(OBJDIR)/%.o: $(SRCDIR)/%.C $(INCLUDE)/%.h
	$(CXX) $< $(CFLAGS) $(OFLAGS) -o $@ 

.cpp: $(OBJS)
	$(CXX) $(COMPILE_FLAGS) $(INCLUDES) $@.cpp -o $@ $(LIBRARY_DIRS) $(GLFLAGS) $(LIBRARIES)

gl:
	cd $(GL_DIR); make all

clean:
	-rm $(OBJDIR)/*.o $(LIB)/*.a
	cd $(GL_DIR); make clean

archive:
	-rm *.tar.gz
	tar -czvf $@.tar.gz *
