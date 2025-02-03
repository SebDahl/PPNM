# CXXFLAGS += -Wall -Wextra -Werror -std=c++23
LDLIBS += -lstdc++ -lm
.PHONEY: all
all: Out.txt

Out.txt:main
	./main >Out.txt

main: main.o
	$(CXX) $(LDFLAGS) $(LDLIBS) -o main main.o

main.o: week_2_main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o week_2_main.cc

.PHONEY: clean
clean:
	$(RM) sfuns.dll main.exe Out.txt        
