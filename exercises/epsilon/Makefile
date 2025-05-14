# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++23


# Automatically find all .cpp files
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# Target executable
TARGET = epsilon  # Change this if you want a different name

# Default target
all: $(TARGET) Out.txt  # Added Out.txt as a dependency to 'all'

# Linking step (uses g++)
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Compilation step for each .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Create Out.txt by running the executable
Out.txt: $(TARGET)  # Out.txt depends on the executable
	./$(TARGET) > Out.txt  # Run the executable and redirect output to Out.txt

# Cleanup
.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) Out.txt

# Test which compiler is used
test:
	@echo $(CXX)