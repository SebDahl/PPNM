# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23 -ffast-math
LDFLAGS =
LDLIBS = -lm

# Automatically find all .cpp files
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# Target executable
TARGET = genlist  # Change this if you want a different name

# Default target
all: $(TARGET)
	@echo ""
	@echo "Build complete. To run the program:"
	@echo "  ./genlist < inputfile > Out.txt" (Out.txt optional)
	@echo ""

# Linking step (uses g++)
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Compilation step for each .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Create Out.txt by running the executable
# Out.txt: $(TARGET)  # Out.txt depends on the executable
# 	./$(TARGET) > Out.txt  # Run the executable and redirect output to Out.txt

# Cleanup
.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) Out.txt

# Test which compiler is used
test:
	@echo $(CXX)



