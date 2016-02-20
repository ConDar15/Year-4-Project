#Compiler and basic flags
CC=gcc
CFLAGS=-std=c11

#Directories used
OBJDIR=obj
OUTDIR=out
LIBDIR=lib

#Library linking options
MPFRLIB=-lgmp -lmpfr
MATHLIB=$(MPFRLIB) -lm
UTILLIB=-L$(LIBDIR) -lutil

#Used to compile a given file with the main program
EXE=-D COMPILE_MAIN

#The lists of files that are used or created
INFILES =bisect_root.c exact_root.c newton_inv_sqrt.c newton_sqrt.c
OUTFILES=$(addprefix out/, $(INFILES:.c=.out))
OBJECTS=$(addprefix obj/, $(INFILES:.c=.o))
LU=$(LIBDIR)/libutil.a

#Defualt option to build all the files
all: $(OBJECTS) $(LU) $(OUTFILES)

#Cleans the workspace, best used for a fresh start
clean:
	rm $(OBJDIR)/*.o $(OUTDIR)/*.out $(LIBDIR)/*.a 

#The following two are used to built the utilities library
$(OBJDIR)/util.o: utilities.c utilities.h
	$(CC) $(CFLAGS) -c $< -o $@
	
$(LU): $(OBJDIR)/util.o
	ar -cr $@ $(OBJDIR)/util*.o	
	
#How to compile a c file to an object file	
$(OBJECTS): $(notdir $(@:.o=.c)) $(notdir $(@:.o=.h)) utilities.h
	$(CC) $(CFLAGS) -c $(notdir $(@:.o=.c)) -o $@

#The following are how to compile each of the OUTFILES
#	The accompanying entry is a shorthand for the first
$(OUTDIR)/bisect_root.out: bisect_root.c bisect_root.h utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) bisect_root.c \
		$(MATHLIB) $(UTILLIB) -o $(OUTDIR)/bisect_root.out
bisect_root: $(OUTDIR)/bisect_root.out

$(OUTDIR)/exact_root.out: exact_root.c exact_root.h utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) exact_root.c \
		$(MPFRLIB) $(UTILLIB) -o $(OUTDIR)/exact_root.out
exact_root: $(OUTDIR)/exact_root.out

$(OUTDIR)/newton_inv_sqrt.out: newton_inv_sqrt.c newton_inv_sqrt.h \
						   utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) newton_inv_sqrt.c \
		$(MATHLIB) $(UTILLIB) -o $(OUTDIR)/newton_inv_sqrt.out
newton_inv_sqrt: $(OUTDIR)/newton_inv_sqrt.out

$(OUTDIR)/newton_sqrt.out: newton_sqrt.c newton_sqrt.h exact_root.h \
					 	   utilities.h $(OBJDIR)/exact_root.o $(LU)
	$(CC) $(CFLAGS) $(EXE) newton_sqrt.c $(OBJDIR)/exact_root.o\
		$(MATHLIB) $(UTILLIB) -o $(OUTDIR)/newton_sqrt.out
newton_sqrt: $(OUTDIR)/newton_sqrt.out
