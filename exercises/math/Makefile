# CXXFLAGS += -Wall -Wextra -Werror -std=c++23
# LDLIBS += -lstdc++ -lm
# .PHONEY: all
# all: Out.txt

# Out.txt:main.
# 	./main >Out.txt

# main: main.o
# 	$(CXX) $(LDFLAGS) $(LDLIBS) -o main main.o

# main.o: math.cc
# 	$(CXX) $(CXXFLAGS) -c -o main.o math.cc

# .PHONEY: clean
# clean:
# 	$(RM) sfuns.dll main.exe Out.txt        

CXXFLAGS += -Wall -Wextra -Werror -std=c++23 -O2
LDLIBS += -lstdc++ -lm
.PHONEY : all
all: Out.txt Out++.txt

Out.txt : main.exe
	mono $< > $@
	#mono main.exe > Out.txt

Out++.txt : main
	./main > Out++.txt

main : main.o sfuns.o
	$(CXX) $(LDFLAGS) $(LDLIBS) -o main $^

main.o : main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o main.cc


sfuns.dll : sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs

.PHONEY: clean
clean:
	$(RM) sfuns.dll main.exe [Oo]ut.txt *.o main

test:
	@echo $(RM) 
	@echo $(CXX) 