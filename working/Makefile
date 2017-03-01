CPP = g++
FLAGS = -ggdb -Wall -I../include/ -O0
LIBDIR = ../lib

particle-tracker-ncams: particle-tracker-ncams.cpp
	$(CPP) $(FLAGS) ../lib/WesleyanCPV.o ../lib/GDF.o ../lib/Calibration.o ../lib/Camera.o ../lib/Frame.o ../lib/Matrix.o ../lib/ParticleFinder.o ../lib/Position.o ../lib/Track.o ../lib/Tracker.o -o $@ particle-tracker-ncams

clean: 
	rm -f particle-tracker-ncams
    rm -f *.cpp~ *.txt~
	rm -f Makefile~
