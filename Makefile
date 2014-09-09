# Makefile for tlODE_Solver with Transmission Line ODE System that writes to Python file

# Compilers, linkers, compile options, link options

CPP = g++ -O2

# Lists of Files

CPP_FILES = main.cpp tlODEsolver.cpp pyWriter.cpp transLine.cpp

HEADERS = ODEsolver.h pyWriter.h transLine.h

OUT_EX = TransLineSim

ALL_SOURCES = $(CPP_FILES) $(HEADERS) Makefile 

$(OUT_EX): $(CPP_FILES) Makefile
	$(CPP) -o $(OUT_EX) $(CPP_FILES)

transLinePlot.py: $(OUT_EX)
	./$(OUT_EX)

run: transLinePlot.py
	python transLinePlot.py

tarball: $(ALL_SOURCES)  
	tar -cvf Trans_Line_Files.tar $(ALL_SOURCES) 