CXX = g++
CXXFLAGS = -O3 -march=native -ffast-math -funroll-loops -DNDEBUG -std=c++11 -I/usr/include/eigen3
TARGET = main_simulation_3d.out
OBJS = main_simulation_3d.o Cilindro3D.o MallaBurbujas3D.o Burbuja3D.o

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main_simulation_3d.o: main_simulation_3d.cpp Cilindro3D.h MallaBurbujas3D.h
	$(CXX) $(CXXFLAGS) -c main_simulation_3d.cpp

Cilindro3D.o: Cilindro3D.cpp Cilindro3D.h
	$(CXX) $(CXXFLAGS) -c Cilindro3D.cpp

MallaBurbujas3D.o: MallaBurbujas3D.cpp MallaBurbujas3D.h Cilindro3D.h Burbuja3D.h
	$(CXX) $(CXXFLAGS) -c MallaBurbujas3D.cpp

Burbuja3D.o: Burbuja3D.cpp Burbuja3D.h
	$(CXX) $(CXXFLAGS) -c Burbuja3D.cpp

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean