CXXFLAGS += -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm

.PHONY: all clean

all: math Out.txt

# Compile both math.cc and sfuns.cc
math: math.cc sfuns.cc
	$(CXX) $(CXXFLAGS) -o math math.cc sfuns.cc $(LDLIBS)
	./math > Out.txt

clean:
	$(RM) sfuns.dll *.exe Out.txt *.o main math
