OBJECT=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
INCLUDE=$(wildcard *.h)
CFLAGS=-I. -Wall -g -O3 -pthread -std=c++0x

tiLDA : $(OBJECT) $(INCLUDE)
	g++ $(CFLAGS) -o tiLDA $(OBJECT)

%.o : %.cpp $(INCLUDE)
	g++ $(CFLAGS) -o $@ -c $<

clean :
	rm -f *.o tiLDA
