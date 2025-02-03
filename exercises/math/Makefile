CXXFLAGS += -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm

.PHONY: all clean

# Default target (optional, builds main if present)
all: main Out.txt

# Rule to build an executable from a .cc file
%: %.cc
	$(CXX) $(CXXFLAGS) $(LDLIBS) -o $@ $< 
	./$@ > Out.txt   # Run the executable and save output to Out.txt

# Clean rule
clean:
	$(RM) sfuns.dll *.exe Out.txt *.o main math
