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

    // Para guardar resultados
    std::string carpeta_resultados;
    int paso_actual;

    // Parámetros de física mejorados
    double factor_conveccion_metal;
    double factor_conveccion_agua;

public:
    MallaBurbujas3D(Cilindro3D& cil, int max_burbujas, const std::string& carpeta = "resultados");
    
    // Generación de burbujas en superficies
    void generarBurbujas();
    
    // Movimiento con física mejorada que considera coeficientes diferentes
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

    // Guardar datos
    void guardarInfluenciaTrayectorias(int paso) const;
    void guardarGeometriaOlla(const Cilindro3D& cilindro) const;
    void guardarEstadoBurbujas(int paso) const;
    void guardarTemperatura(int paso, const Cilindro3D& cilindro) const;
    void guardarCoeficientesDifusividad(const Cilindro3D& cilindro) const; // NUEVO MÉTODO
    
    // Crear la carpeta de resultados
    void crearCarpetaResultados() const;
    
    // Getter para el paso actual
    int getPasoActual() const { return paso_actual; }
    void incrementarPaso() { paso_actual++; }

    // Configurar parámetros de física
    void setFactoresConveccion(double factor_metal, double factor_agua);

private:
    // Métodos auxiliares
    double calcularProbabilidadGeneracion(double h, double r, double theta) const;
    void calcularVelocidadBurbuja(Burbuja3D& burbuja);
    void aplicarConveccion(Burbuja3D& burbuja);
    void aplicarInfluenciaTrayectorias(Burbuja3D& burbuja);
    void aplicarEfectoCoeficientes(Burbuja3D& burbuja); // NUEVO MÉTODO
    bool estaDentroCilindro(const Burbuja3D& burbuja) const;
    
    // Generar posición aleatoria en superficies
    std::vector<double> generarPosicionBase();
    std::vector<double> generarPosicionPared();
    
    // Método para obtener tipo de material en posición de burbuja
    std::string getTipoMaterial(const Burbuja3D& burbuja) const;
};

#endif