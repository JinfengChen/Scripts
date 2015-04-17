oppclass: main.o rectangle.o
	g++ main.o rectangle.o -o oppclass

main.o: main.cpp rectangle.h
	g++ -c main.cpp

rectangle.o: rectangle.cpp rectangle.h
	g++ -c rectangle.cpp

clean: 
	@echo "cleaning project"
	-rm oppclass *.o
	@echo "clean completed"

.PHONY: clean

