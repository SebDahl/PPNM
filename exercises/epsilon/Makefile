CXXFLAGS += -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm

.PHONY: all clean

all: epsilon Out.txt

# Compile both math.cc and sfuns.cc
math: epsilon.cc sfuns.c
	$(CXX) $(CXXFLAGS) -o epsilon epsilon.cc $(LDLIBS)
	./epsilon > Out.txt

clean:
	$(RM) sfuns.dll *.exe Out.txt *.o main epsilon
