#ifndef CILINDRO3D_H
#define CILINDRO3D_H

#include <vector>
#include <functional>
#include <cmath>

class Cilindro3D {
private:
    int radio;          // Radio del cilindro
    int altura;         // Altura del cilindro
    double dt;          // Paso temporal
    double alpha;       // Difusividad térmica
    
    // Malla de temperatura en coordenadas cilíndricas [altura][radio][angulo]
    // Por simetría axial, solo necesitamos 2D [altura][radio]
    std::vector<std::vector<double>> temperatura;
    std::vector<std::vector<double>> temperatura_nueva;
    
    // Fuente de calor en la base
    std::vector<std::vector<double>> fuente_calor;

public:
    Cilindro3D(int r, int h, double delta_t, double alpha_term);
    
    // Configurar distribución inicial de temperatura
    void setDistribucionTemperaturaInicial(const std::vector<std::vector<double>>& temp_init);
    
    // Configurar fuente de calor en la base
    void setFuenteCalorBase(const std::vector<std::vector<double>>& fuente);
    
    // Método para evolucionar la temperatura en el tiempo
    void evolucionarTemperatura();
    
    // Obtener temperatura en un punto específico (coordenadas cilíndricas)
    double getTemperatura(double h, double r, double theta = 0) const;
    
    // Obtener gradiente de temperatura en 3D
    std::vector<double> getGradienteTemperatura(double h, double r, double theta = 0) const;
    
    // Método para diferentes distribuciones de calor
    void configurarFuenteUniforme(double potencia);
    void configurarFuenteGaussiana(double potencia_centro, double sigma);
    void configurarFuenteAnular(double potencia, double radio_int, double radio_ext);
    
    // Conversión entre coordenadas
    std::vector<double> cartesianToCilindricas(double x, double y, double z) const;
    std::vector<double> cilindricasToCartesian(double h, double r, double theta) const;
    
    // Verificar si un punto está dentro del cilindro
    bool estaDentro(double x, double y, double z) const;
    
    // Getters
    int getRadio() const { return radio; }
    int getAltura() const { return altura; }
    double getDT() const { return dt; }
    
    // Visualización
    void imprimirTemperaturaBase() const;
};

#endif