CC=pgc++
OPT=-fast -Minfo=accel

all :	main.o input.o
	${CC} -std=c++11 *.o -o ./wvesim $(OPT)

main.o : main.cpp
	${CC} -std=c++11 -c $(OPT) main.cpp

input.o : input.cpp
	${CC} -std=c++11 -c $(OPT) input.cpp

clean :
	rm *.o
