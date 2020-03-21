#objects = main.o coordinate.o segment.o tracer.o
magneto: main.cpp coordinate.cpp event.cpp obstacle.cpp segment.cpp simbox.cpp tracer.cpp event.cpp fluid.cpp measurement.cpp
	g++  -O2 -std=c++11 -o magneto.x main.cpp obstacle.cpp coordinate.cpp  segment.cpp simbox.cpp tracer.cpp event.cpp fluid.cpp measurement.cpp
clean:
	rm magneto.x
