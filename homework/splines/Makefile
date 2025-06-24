CXX = g++
CXXFLAGS += -std=c++23 -O2 # -Wall

TARGETS = Out.txt part_a_int.svg
EXE = main.exe

all: $(TARGETS)

# Run main.exe to generate Out.txt
Out.txt: $(EXE)
	./$(EXE)

# Compile main.exe
$(EXE): main.cpp QR.o matrix.o EVD.o
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

part_a_int.svg: main.exe part_a_int.gpi
	gnuplot part_a_int.gpi


clean:
	$(RM) $(EXE) Out.txt *.o *.txt *.svg
