#include "Burbuja3D.h"

Burbuja3D::Burbuja3D(double x, double y, double z, double size, int identificador) 
    : pos_x(x), pos_y(y), pos_z(z), tamano(size), id(identificador), 
      activa(true), vel_x(0), vel_y(0), vel_z(0) {
    trayectoria.push_back({x, y, z});
}

// Agregar en la trayectoria de la burbuja su nueva posicion
void Burbuja3D::mover(double nueva_x, double nueva_y, double nueva_z) {
    pos_x = nueva_x;
    pos_y = nueva_y;
    pos_z = nueva_z;
    trayectoria.push_back({nueva_x, nueva_y, nueva_z});
}

// La nueva posicion de la burbuja dada una velocidad y guardarla en su trayectoria
void Burbuja3D::moverConVelocidad(double dt) {
    pos_x += vel_x * dt;
    pos_y += vel_y * dt;
    pos_z += vel_z * dt;
    trayectoria.push_back({pos_x, pos_y, pos_z});
}

// Setter para la velocidad
void Burbuja3D::setVelocidad(double vx, double vy, double vz) {
    vel_x = vx;
    vel_y = vy;
    vel_z = vz;
}

// Coalescencia con conservacion del volumen y momento de las burbujas
void Burbuja3D::coalescer(const Burbuja3D& otra) {
    // Conservación de masa (volumen proporcional al cubo del radio)
    double vol1 = tamano * tamano * tamano;
    double vol2 = otra.tamano * otra.tamano * otra.tamano;
    double vol_total = vol1 + vol2;
    
    tamano = std::cbrt(vol_total);
    
    // Conservación de momento (promedio ponderado de velocidades)
    double masa_total = vol1 + vol2;
    if (masa_total > 0) {
        vel_x = (vol1 * vel_x + vol2 * otra.vel_x) / masa_total;
        vel_y = (vol1 * vel_y + vol2 * otra.vel_y) / masa_total;
        vel_z = (vol1 * vel_z + vol2 * otra.vel_z) / masa_total;
    }
}

// Calcular la distancia entre dos burbujas
double Burbuja3D::distancia(const Burbuja3D& otra) const {
    double dx = pos_x - otra.pos_x;
    double dy = pos_y - otra.pos_y;
    double dz = pos_z - otra.pos_z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}