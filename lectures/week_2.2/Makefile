CXXFLAGS += -Wall -Wextra -std=c++23 -ferror-limit=1
LDLIBS += -lstdc++ -lm
CODE = $(filter %.cc,$^)
DLLs =$(addprefix -reference ,$(filter %.dll,$^))
MKLIB = mcs -out:$@ -target:library $(DLLS) $(CODE)
MKEXE =mcs -out:$@ -target:exe $(DLLS) $(CODE)

Out.txt: main.exe main
	mono main.exe > $@
	./main > $@

main.exe: main.cs ; $(MKEXE)
