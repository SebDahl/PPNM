CXX = g++
CXXFLAGS += -std=c++23 -O2 # -Wall

TARGETS = Out.txt part_a.svg
EXE = main.exe

all: $(TARGETS)

# Run main.exe to generate Out.txt
Out.txt: $(EXE)
	./$(EXE)

# Compile main.exe
$(EXE): main.cpp QR.o matrix.o EVD.o sfuns.o ode.o
	$(CXX) $(CXXFLAGS) main.cpp QR.o matrix.o EVD.o -o $(EXE)

# Compile QR.o
QR.o: QR.cpp QR.h matrix.h
	$(CXX) $(CXXFLAGS) -c QR.cpp -o QR.o

# Compile EVD.o
EVD.o: EVD.cpp EVD.h matrix.h
	$(CXX) $(CXXFLAGS) -c EVD.cpp -o EVD.o

# Compile matrix.o
matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

# Compile ode.o
ode.o: ode.cpp ode.h matrix.h
	$(CXX) $(CXXFLAGS) -c ode.cpp -o ode.o

# Compile sfuns.o
sfuns.o: sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c sfuns.cpp -o sfuns.o

part_a.svg : main.exe result_list.txt deviation_list.txt N_s.txt part_a.gpi
	gnuplot part_a.gpi


clean:
	$(RM) $(EXE) Out.txt *.o *.txt *.svg
