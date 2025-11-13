compilar:
	g++ -O3 -march=native -ffast-math -funroll-loops -DNDEBUG -I/usr/include/eigen3 Cilindro3D.h -o Cilindro3D.out

clean:
	rm -f *.out *.pdf *.xyz