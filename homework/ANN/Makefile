CXX = g++
CXXFLAGS += -std=c++23 -O2

TARGETS = main.exe main_b.exe part_a.svg
OBJS = QR.o matrix.o EVD.o sfuns.o ode.o minimizer.o

all: $(TARGETS)

# Main target 1
main.exe: main.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Main target 2
main_b.exe: main_b.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Run main.exe to generate Out.txt
Out.txt: main.exe
	./main.exe > Out.txt

# Generate part_a.svg if needed
part_a.svg: main.exe values.txt part_a.gpi
	gnuplot part_a.gpi

# Object file rules
QR.o: QR.cpp QR.h matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

EVD.o: EVD.cpp EVD.h matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

ode.o: ode.cpp ode.h matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

sfuns.o: sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

minimizer.o: minimizer.cpp minimizer.h matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	$(RM) main.exe main_b.exe Out.txt *.o *.txt *.svg
