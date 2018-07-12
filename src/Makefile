masscomp: masscomp.o tinyxml2.o
	g++ masscomp.o tinyxml2.o -o masscomp

masscomp.o: masscomp.cpp head.h
	g++ -c masscomp.cpp -o masscomp.o

tinyxml2.o: tinyxml2.cpp tinyxml2.h
	g++ -c tinyxml2.cpp -o tinyxml2.o

clean:
	rm *.o masscomp

