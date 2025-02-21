CXX = g++
CXXFLAGS += -std=c++23 -Wall -O2
LDLIBS += -lm

# Executable and data files
TARGET = plot.test.svg
DATAFILE = data.txt
EXE = main.exe
OBJS = main.o sfuns.o

# Default rule
all: $(TARGET)

# Rule to generate the plot using Gnuplot
plot.test.svg: test.gpi $(DATAFILE)
	gnuplot $<

# Rule to generate the data file from the C++ program
$(DATAFILE): $(EXE)
	./$(EXE) > $@

# Rule to compile the C++ program
$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Compilation rules for object files
main.o: main.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

sfuns.o: sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	$(RM) *.exe *.svg $(DATAFILE) *.o
