# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -O2 -pthread

# Automatically find all .cpp files
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# Target executable
TARGET = main

# Default target
all: $(TARGET) Out.txt

# Linking step
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Compilation step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Measure performance using POSIX time
N = 100000000
TIME = time --portability --append --output $@

Out.txt: $(TARGET) Makefile
	@echo "Running performance tests..." > $@
	$(TIME) ./$(TARGET) -terms:$(N) -threads:1 >> $@
	$(TIME) ./$(TARGET) -terms:$(N) -threads:2 >> $@
	$(TIME) ./$(TARGET) -terms:$(N) -threads:3 >> $@
	$(TIME) ./$(TARGET) -terms:$(N) -threads:4 >> $@

# Cleanup
.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) Out.txt

# Test compiler
test:
	@echo $(CXX)
