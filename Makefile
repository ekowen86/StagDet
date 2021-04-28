CXX=g++
CXXFLAGS=-std=c++11 -O3 -Wall

all: stag_det

stag_det: stag_det.cc
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:
	$(RM) stag_det
