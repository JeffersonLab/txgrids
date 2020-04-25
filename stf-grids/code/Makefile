CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELINCS = $(shell apfel-config --cppflags)
APFELLIBS = $(shell apfel-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELINCS)

CLIBS += $(LHAPDFLIBS) $(APFELLIBS)

install : all
all : PDFGrid StructureFunctionGrid SumRules

PDFGrid: PDFGrid.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

StructureFunctionGrid: StructureFunctionGrid.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

SumRules: SumRules.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la PDFGrid StructureFunctionGrid SumRules *~
