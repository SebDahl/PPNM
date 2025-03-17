CXX = g++
CXXFLAGS += -std=c++23 -O2 # -Wall

TARGETS = Out.txt E0rmax.svg
EXE = main.exe
HYDROGEN = hydrogen.exe

all: $(TARGETS)

# Run main.exe to generate Out.txt
Out.txt: $(EXE)
	./$(EXE)

# Compile main.exe
$(EXE): main.cpp QR.o matrix.o EVD.o
	$(CXX) $(CXXFLAGS) main.cpp QR.o matrix.o EVD.o -o $(EXE)

# Compile hydrogen.exe
$(HYDROGEN): hydrogen.cpp QR.o matrix.o EVD.o
	$(CXX) $(CXXFLAGS) hydrogen.cpp QR.o matrix.o EVD.o -o $(HYDROGEN)

# Compile QR.o
QR.o: QR.cpp QR.h matrix.h
	$(CXX) $(CXXFLAGS) -c QR.cpp -o QR.o

# Compile EVD.o
EVD.o: EVD.cpp EVD.h matrix.h
	$(CXX) $(CXXFLAGS) -c EVD.cpp -o EVD.o

# Compile matrix.o
matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

# Generate E0rmax.svg using hydrogen.exe
E0rmax.svg: hydrogen.exe E0s_rmax.txt rmaxs.txt E0rmax.gpi
	gnuplot E0rmax.gpi

# Ensure these files exist before plotting
E0s_rmax.txt rmaxs.txt: hydrogen.exe
	@./hydrogen.exe -rmax 10 -dr 0.3

# Generate out.times.data by running hydrogen.exe
# out.times.data: hydrogen.exe
# 	./hydrogen.exe > out.times.data

clean:
	$(RM) $(EXE) $(HYDROGEN) Out.txt *.o *.txt *.data E0rmax.svg
