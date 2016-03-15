#Compiler and basic flags
CC=gcc
CFLAGS=-std=c11

#Directories used
OBJDIR=obj
OUTDIR=out
LIBDIR=lib
TSTDIR=test

#Library linking options
MPFRLIB=-lgmp -lmpfr
MATHLIB=$(MPFRLIB) -lm
UTILLIB=-L$(LIBDIR) -lutil

#Used to compile a given file with the main program
EXE=-D COMPILE_MAIN

#The lists of files that are used or created
INFILES  =bisect_root.c exact_root.c newton_inv_sqrt.c newton_sqrt.c \
		  geometric_trig.c geometric_inv_trig.c taylor_trig.c \
		  taylor_inv_trig.c cordic_trig.c
TESTFILES=test_newton_sqrt.c
OUTFILES=$(addprefix $(OUTDIR)/, $(INFILES:.c=.out))
OBJECTS =$(addprefix $(OBJDIR)/, $(INFILES:.c=.o))
TESTS   =$(addprefix $(TSTDIR)/, $(TESTFILES:.c=.out))
LU=$(LIBDIR)/libutil.a

#Defualt option to build all the files
all: $(OBJECTS) $(LU) $(OUTFILES) $(TESTS)

#Cleans the workspace, best used for a fresh start
clean:
	rm $(OBJDIR)/*.o $(OUTDIR)/*.out $(LIBDIR)/*.a 

#The following few are used to built the utilities library
$(OBJDIR)/util.o: utilities.c utilities.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/util_trig.o: trig_utilities.c trig_utilities.h
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/util_test.o: $(TSTDIR)/test_utilities.c $(TSTDIR)/test_utilities.h
	$(CC) $(CFLAGS) -c $< -o $@
	
$(LU): $(OBJDIR)/util.o $(OBJDIR)/util_trig.o $(OBJDIR)/util_test.o
	ar -cr $@ $(OBJDIR)/util*.o	
	
#How to compile a c file to an object file	
$(OBJECTS): $(notdir $(@:.o=.c)) $(notdir $(@:.o=.h)) utilities.h\
		    trig_utilities.h
	$(CC) $(CFLAGS) -c $(notdir $(@:.o=.c)) -o $@

#The following are how to compile each of the OUTFILES
#	The accompanying entry is a shorthand for the first

## SQUARE ROOT FILES ##
$(OUTDIR)/bisect_root.out: bisect_root.c bisect_root.h utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) bisect_root.c \
		$(MATHLIB) $(UTILLIB) -o $@
bisect_root: $(OUTDIR)/bisect_root.out

$(OUTDIR)/exact_root.out: exact_root.c exact_root.h utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) exact_root.c \
		$(MPFRLIB) $(UTILLIB) -o $@
exact_root: $(OUTDIR)/exact_root.out

$(OUTDIR)/newton_inv_sqrt.out: newton_inv_sqrt.c newton_inv_sqrt.h \
						   	   utilities.h $(LU)
	$(CC) $(CFLAGS) $(EXE) newton_inv_sqrt.c \
		$(MATHLIB) $(UTILLIB) -o $@
newton_inv_sqrt: $(OUTDIR)/newton_inv_sqrt.out

$(OUTDIR)/newton_sqrt.out: newton_sqrt.c newton_sqrt.h exact_root.h \
					 	   utilities.h $(OBJDIR)/exact_root.o $(LU)
	$(CC) $(CFLAGS) $(EXE) newton_sqrt.c $(OBJDIR)/exact_root.o\
		$(MATHLIB) $(UTILLIB) -o $@
newton_sqrt: $(OUTDIR)/newton_sqrt.out

## TRIGONOMETRIC FILES ##
$(OUTDIR)/geometric_trig.out: geometric_trig.c trig_utilities.h \
							  geometric_trig.h $(LU)
	$(CC) $(CFLAGS) $(EXE) geometric_trig.c $(UTILLIB)\
		$(MATHLIB) -o $@
geometric_trig: $(OUTDIR)/geometric_trig.out

$(OUTDIR)/geometric_inv_trig.out: geometric_inv_trig.c trig_utilities.h \
							  	  geometric_inv_trig.h $(LU)
	$(CC) $(CFLAGS) $(EXE) geometric_inv_trig.c $(UTILLIB)\
		$(MATHLIB) -o $@
geometric_inv_trig: $(OUTDIR)/geometric_inv_trig.out

$(OUTDIR)/taylor_trig.out: taylor_trig.c trig_utilities.h \
						   taylor_trig.h $(LU)
	$(CC) $(CFLAGS) $(EXE) taylor_trig.c $(UTILLIB)\
		$(MATHLIB) -o $@
geometric_trig: $(OUTDIR)/taylor_trig.out

$(OUTDIR)/taylor_inv_trig.out: taylor_inv_trig.c trig_utilities.h \
						   	   taylor_inv_trig.h $(LU)
	$(CC) $(CFLAGS) $(EXE) taylor_inv_trig.c $(UTILLIB)\
		$(MATHLIB) -o $@
geometric_trig: $(OUTDIR)/taylor_inv_trig.out

$(OUTDIR)/cordic_trig.out: cordic_trig.c trig_utilities.h trig_fixed.h \
						   cordic_trig.h $(LU)
	$(CC) $(CFLAGS) $(EXE) cordic_trig.c $(UTILLIB) $(MPFRLIB) -o $@
cordic_trig: $(OUTDIR)/cordic_trig.out

## TESTING FILES ##
$(TSTDIR)/test_newton_sqrt.out: $(TSTDIR)/test_newton_sqrt.c \
				$(TSTDIR)/test_utilities.h \
				$(OBJDIR)/newton_sqrt.o $(OBJDIR)/exact_root.o \
				$(LU)
	$(CC) $(CFLAGS) $< $(OBJDIR)/newton_sqrt.o $(OBJDIR)/exact_root.o \
		$(MPFRLIB) $(UTILLIB) -I. -o $@
test_newton_sqrt: $(TSTDIR)/test_newton_sqrt.out
