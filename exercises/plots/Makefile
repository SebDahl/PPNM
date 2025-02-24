CXX = g++
CXXFLAGS += -std=c++23 -Wall -O2
LDLIBS += -lm

TARGETS = erf_plot.svg gamma_plot.svg lngamma_plot.svg
EXE = main.exe

all: $(TARGETS)

# Generate both data files in one run
erf_data.txt gamma_data.txt lngamma_data.txt: $(EXE)
	./$(EXE)

# Plot the error function
erf_plot.svg: erf.gpi erf_data.txt
	gnuplot $<

# Plot the gamma function
gamma_plot.svg: gamma.gpi gamma_data.txt
	gnuplot $<
#Plot lngamma function
lngamma_plot.svg: lngamma.gpi lngamma_data.txt
	gnuplot $<

# Compile main.cpp
$(EXE): main.cpp sfuns.o
	$(CXX) $(CXXFLAGS) main.cpp sfuns.o -o $(EXE) $(LDLIBS)

# Compile sfuns.cpp
sfuns.o: sfuns.cpp sfuns.h
	$(CXX) $(CXXFLAGS) -c sfuns.cpp -o sfuns.o

clean:
	$(RM) *.exe *.svg *.txt *.o
