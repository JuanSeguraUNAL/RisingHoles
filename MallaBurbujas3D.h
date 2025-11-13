#ifndef MALLA_BURBUJAS_3D_H
#define MALLA_BURBUJAS_3D_H

#include "Cilindro3D.h"
#include "Burbuja3D.h"
#include <vector>
#include <random>
#include <functional>

class MallaBurbujas3D {
private:
    Cilindro3D& cilindro;
    int max_burbujas_por_paso;
    int next_burbuja_id;
    
    std::vector<Burbuja3D> burbujas;
    
    // Mapa 3D de influencia de trayectorias (discretizado)
    std::vector<std::vector<std::vector<double>>> influencia_trayectorias;
    
    // Generadores aleatorios
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
    std::uniform_real_distribution<double> dist_angular;

public:
    MallaBurbujas3D(Cilindro3D& cil, int max_burbujas);
    
    // Generación de burbujas en superficies
    void generarBurbujas();
    
    // Movimiento con física mejorada
    void moverBurbujas();
    
    // Coalescencia con radio de influencia
    void verificarCoalescencia();
    
    // Limpiar burbujas que salieron
    void limpiarBurbujas();
    
    // Getters
    const std::vector<Burbuja3D>& getBurbujas() const { return burbujas; }
    int getNumBurbujasActivas() const;
    
    // Visualización
    void imprimirEstado() const;

private:
    // Métodos auxiliares
    double calcularProbabilidadGeneracion(double h, double r, double theta) const;
    void calcularVelocidadBurbuja(Burbuja3D& burbuja);
    void aplicarConveccion(Burbuja3D& burbuja);
    void aplicarInfluenciaTrayectorias(Burbuja3D& burbuja);
    bool estaDentroCilindro(const Burbuja3D& burbuja) const;
    
    // Generar posición aleatoria en superficies
    std::vector<double> generarPosicionBase();
    std::vector<double> generarPosicionPared();
};

#endif