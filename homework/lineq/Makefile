# CXX = g++
# CXXFLAGS += -std=c++23 -O2 # -Wall

# TARGETS = Out.txt
# EXE = main.exe

# all: $(TARGETS)

# # Run the executable to generate Out.txt
# Out.txt: $(EXE)
# 	./$(EXE)

# # Compile main.exe
# $(EXE): main.cpp QR.o matrix.o
# 	$(CXX) $(CXXFLAGS) main.cpp QR.o matrix.o -o $(EXE)

# # Compile QR.o
# QR.o: QR.cpp QR.h matrix.h
# 	$(CXX) $(CXXFLAGS) -c QR.cpp -o QR.o

# # Compile matrix.o
# matrix.o: matrix.cpp matrix.h
# 	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

# clean:
# 	$(RM) $(EXE) Out.txt *.o



CXX = g++
CXXFLAGS += -std=c++23 -O3

TARGETS = out.times.data main.exe count.svg
EXE = main.exe
COUNT = count.exe

all: $(TARGETS)

# Rule to generate timing data using POSIX time
out.times.data: $(COUNT)
	@rm -f out.times.data
	@for N in $$(seq 100 20 1000); do \
		time --format "$$N %e" --output out.times.data --append ./$(COUNT) $$N; \
	done

# Compile main.exe if needed
$(EXE): main.cpp QR.o matrix.o
	$(CXX) $(CXXFLAGS) main.cpp QR.o matrix.o -o $(EXE)

# Compile count.exe
$(COUNT): count.cpp QR.o matrix.o
	$(CXX) $(CXXFLAGS) count.cpp QR.o matrix.o -o $(COUNT)

# Compile QR.o
QR.o: QR.cpp QR.h matrix.h
	$(CXX) $(CXXFLAGS) -c QR.cpp -o QR.o

# Compile matrix.o
matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

count.svg: count.gpi out.times.data
	gnuplot $<

clean:
	$(RM) $(EXE) $(COUNT) out.times.data *.o *.svg
