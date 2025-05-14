# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -ffast-math
LDFLAGS =
LDLIBS = -lm

# Automatically find all .cpp files
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
EXE = complex  # Use this instead of 'complex.exe'

# Target executable
TARGET = complex

# Default target
all: $(TARGET) Out.txt

# Linking step
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Compilation step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Generate Out.txt
Out.txt: $(EXE)
	./$(EXE) | tee Out.txt

# Cleanup
.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) Out.txt

# Test compiler
test:
	@echo $(CXX)
