CC=pgc++
OPT=-fast -Minfo=accel

all :	main.o input.o interp1d.o interp2d.o error.o
	${CC} -std=c++11 *.o -o ../input-data/wvesim $(OPT)

main.o : main.cpp
	${CC} -std=c++11 -c $(OPT) main.cpp

input.o : input.cpp
	${CC} -std=c++11 -c $(OPT) input.cpp

interp1d.o : interp1d.cpp
	${CC} -std=c++11 -c $(OPT) interp1d.cpp

interp2d.o : interp2d.cpp
	${CC} -std=c++11 -c $(OPT) interp2d.cpp

error.o : error.cpp
	${CC} -std=c++11 -c $(OPT) error.cpp

clean :
	rm *.o
