#MATLAB = /usr/local/MATLAB/R2011b
MATLAB = /Applications/MATLAB_R2015b.app
#OSTYPE = linux-gnu
OSTYPE = maci64
CC = gcc
CMEX = $(MATLAB)/bin/mex
DEBUG = -O
MEXSUFFIX = 

LDFLAGS = 
INCLUDES = -lgsl -lgslcblas -lm
CFLAGS = $(INCLUDES)
PROGRAM = make_snake
FUNCTION = $(PROGRAM)
MAKEFILE = Makefile
HDRS = 
SRCS = snake.c image.c pnpoly.c spline.c
MAINSRC = $(PROGRAM).c
HEADER = ./include.sh
#LIBS = $(LOCAL_LIBS) $(SYS_LIBS) -L/home/etudiant/lib -I/home/etudiant/include
LIBS = $(LOCAL_LIBS) $(SYS_LIBS) -L/Users/cathy/local/lib -I/Users/cathy/local/include


all: snake.h $(FUNCTION)

$(FUNCTION): $(MAKEFILE)

sources: $(SRCS) $(HDRS)

$(FUNCTION): $(MAINSRC) $(SRCS) $(HDRS) Makefile
	$(CMEX) $(DEBUG) $(CFLAGS) $(LDFLAGS) $(SRCS) $(MAINSRC) $(LIBS) -output $(FUNCTION)

delmex:
	rm -f $(FUNCTION)

clean: delmex
	rm -f *~ core .depend.$(OSTYPE)

tags: sources
	etags -t $(SRCS) $(HDRS)

.PHONY: all sources delmex full
.PHONY: clean tags
