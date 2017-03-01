CPP = g++
FLAGS = -g -Wall -I../include/ -O0

all: \
	WesleyanCPV \
	GDF \
	ParticleFinder \
	Position \
	Frame \
	Track \
	Tracker \
	Camera \
	Calibration \
	Matrix \
	Trackfile

WesleyanCPV: WesleyanCPV.cpp ../include/WesleyanCPV.h
	$(CPP) $(FLAGS) -c WesleyanCPV.cpp
	
GDF: GDF.cpp ../include/GDF.h
	$(CPP) $(FLAGS) -c GDF.cpp

ParticleFinder: ParticleFinder.cpp ../include/ParticleFinder.h
	$(CPP) $(FLAGS) -c ParticleFinder.cpp

Position: Position.cpp ../include/Position.h
	$(CPP) $(FLAGS) -c Position.cpp

Frame: Frame.cpp ../include/Frame.h
	$(CPP) $(FLAGS) -c Frame.cpp

Track: Track.cpp ../include/Track.h
	$(CPP) $(FLAGS) -c Track.cpp

Tracker: Tracker.cpp ../include/Tracker.h
	$(CPP) $(FLAGS) -c Tracker.cpp

Camera: Camera.cpp ../include/Camera.h
	$(CPP) $(FLAGS) -c Camera.cpp

Calibration: Calibration.cpp ../include/Calibration.h
	$(CPP) $(FLAGS) -c Calibration.cpp

Matrix: Matrix.cpp ../include/Matrix.h
	$(CPP) $(FLAGS) -c Matrix.cpp

Trackfile: Trackfile.cpp ../include/Trackfile.h
	$(CPP) $(FLAGS) -c Trackfile.cpp

clean:
	rm -f *.o
	rm -f *.cpp~
