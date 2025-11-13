#ifndef BURBUJA3D_H
#define BURBUJA3D_H

#include <vector>
#include <cmath>

class Burbuja3D {
private:
    // Posición en coordenadas cartesianas
    double pos_x, pos_y, pos_z;
    double tamano;
    int id;
    bool activa;
    
    // Velocidad de la burbuja
    double vel_x, vel_y, vel_z;
    
    // Historial de posiciones
    std::vector<std::vector<double>> trayectoria;

public:
    Burbuja3D(double x, double y, double z, double size, int identificador);
    
    // Movimiento en 3D
    void mover(double nueva_x, double nueva_y, double nueva_z);
    void moverConVelocidad(double dt);
    
    // Coalescencia
    void coalescer(const Burbuja3D& otra);
    
    // Getters
    double getPosX() const { return pos_x; }
    double getPosY() const { return pos_y; }
    double getPosZ() const { return pos_z; }
    double getTamano() const { return tamano; }
    int getId() const { return id; }
    bool isActiva() const { return activa; }
    double getVelX() const { return vel_x; }
    double getVelY() const { return vel_y; }
    double getVelZ() const { return vel_z; }
    const std::vector<std::vector<double>>& getTrayectoria() const { return trayectoria; }
    
    // Setters
    void setActiva(bool estado) { activa = estado; }
    void setTamano(double size) { tamano = size; }
    void setVelocidad(double vx, double vy, double vz);
    
    // Cálculo de distancia entre burbujas
    double distancia(const Burbuja3D& otra) const;
};

#endif