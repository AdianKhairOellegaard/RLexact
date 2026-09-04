OBJECTS = RLexact.o RLio.o RLtables.o RLhamil.o RLsymm.o Diagonal.o nr.o RLutil.o  RLcross.o RLexpect.o RLlancz.o RLmatrix.o RLsparse.o regc.o RLhamil.dipole.o
LIBS = -lm 
FLAGS =  -I. -I/usr/include/g++-3 -O0 -std=c++11 -gdwarf-2 -ansi -g -Wno-deprecated-declarations#-fsanitize=address,undefined -fno-omit-frame-pointer
OPTFLAGS = -I. -I/usr/include/g++-3 -O0 -gdwarf-2 -ansi -std=c++11 -Wno-deprecated-declarations
# FLAGS =  -I. -I/usr/include/g++-3 -O0 -gdwarf-2
# OPTFLAGS =  -I. -I/usr/include/g++-3 -O3 -gdwarf-2
#FLAGS =  -p -I. -I/usr/include/g++-3 -O1 -gdwarf-2 -ansi -pedantic
EXECUTABLE = RLexact
COMPILER = mpic++
RLexact:	$(OBJECTS)
		$(COMPILER) $(FLAGS) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

RLtables.o:	RLexact.h Functions.h RLtables.C 
		$(COMPILER) $(OPTFLAGS) -c RLtables.C

RLexact.o:	RLexact.h Functions.h RLexact.C
		$(COMPILER) $(FLAGS) -c RLexact.C

regc.o:		RLexact.h Functions.h regc.cpp 
		$(COMPILER) $(OPTFLAGS) -c regc.cpp

RLio.o:		RLexact.h Functions.h RLio.C
		$(COMPILER) $(FLAGS) -c RLio.C

RLlancz.o:	RLexact.h Functions.h RLlancz.C
		$(COMPILER) $(OPTFLAGS) -c RLlancz.C

RLmatrix.o:	RLexact.h Functions.h RLmatrix.C
		$(COMPILER) $(OPTFLAGS) -c RLmatrix.C

RLsymm.o:	RLexact.h Functions.h RLsymm.C
		$(COMPILER) $(OPTFLAGS) -c RLsymm.C

RLhamil.o:	RLexact.h Functions.h RLhamil.C
		$(COMPILER) $(OPTFLAGS) -c RLhamil.C

RLhamil_dipole.o:	RLexact.h Functions.h RLhamil.dipole.C
		$(COMPILER) $(OPTFLAGS) -c RLhamil.dipole.C

RLutil.o:	RLexact.h Functions.h RLutil.C
		$(COMPILER) $(OPTFLAGS) -c RLutil.C

Diagonal.o:	RLexact.h Functions.h Diagonal.C
		$(COMPILER) $(OPTFLAGS) -c Diagonal.C

nr.o:		RLexact.h Functions.h cnr.h nr.C
		$(COMPILER) $(OPTFLAGS) -c nr.C

RLcross.o:	RLexact.h Functions.h RLcross.C
		$(COMPILER) $(OPTFLAGS) -c RLcross.C

RLexpect.o:	RLexact.h Functions.h RLexpect.C
		$(COMPILER) $(OPTFLAGS) -c RLexpect.C

RLsparse.o:	RLexact.h Functions.h RLsparse.C
		$(COMPILER) $(OPTFLAGS) -c RLsparse.C

clean:
		rm -f *.o RLexact
