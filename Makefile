CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall

INCLUDES = -I/usr/include/eigen3

TARGET = spice # onoma ektelesimou

SRCS = main.cpp parser.cpp circuit_equations.cpp options.cpp solver.cpp dc_sweep.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) main.cpp -o $(TARGET)

clean:
	rm -f $(TARGET) 
# dc_op.txt dc_sweep_*.txt

.PHONY: all clean