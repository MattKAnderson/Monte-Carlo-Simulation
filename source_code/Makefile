CPP = g++
CPPFLAGS = -O3 -fopenmp -std=c++11
OFLAGS = -c -o
FLAGS = -o
OBJECT = SQclockP.o dataFunctions.o MCsim.o  
HEADER = dataFunctions.h MCsim.h

%.o: %.cpp 
	$(CPP) $(CPPFLAGS) $(OFLAGS) $@ $<

main: $(OBJECT) 
	$(CPP) $(CPPFLAGS) $(FLAGS) $@ $^

clean:
	rm *.o
