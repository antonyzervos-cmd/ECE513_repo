CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

INCLUDES = -I. -I/usr/include/eigen3 -I/usr/include/suitesparse
LIBS = -lcxsparse

TARGET = spice
SRCS = main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET) 
#dc_op.txt dc_sweep_*.txt

.PHONY: all clean