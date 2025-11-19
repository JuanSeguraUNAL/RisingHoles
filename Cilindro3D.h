#ifndef CILINDRO3D_H
#define CILINDRO3D_H

#include <vector>
#include <functional>
#include <cmath>

class Cilindro3D {
private:
    int radio_interior;     // Radio del interior de la olla
    int grosor_pared;       // Grosor de las paredes en puntos
    int radio_total;        // radio_interior + grosor_pared
    int altura;
    double dt;
    double alpha_olla;      // Difusividad térmica del material de la olla (metal)
    double alpha_interior;  // Difusividad térmica del interior (agua/arroz)
    
    // Malla de temperatura en coordenadas cilíndricas [altura][radio_total]
    std::vector<std::vector<double>> temperatura;
    std::vector<std::vector<double>> temperatura_nueva;
    
    // Malla de coeficientes de difusividad [altura][radio_total]
    std::vector<std::vector<double>> alpha_malla;
    
    // Fuente de calor en la base
    std::vector<std::vector<double>> fuente_calor;

public:
    // Constructor modificado para dos coeficientes
    Cilindro3D(int radio_int, int grosor_pared, int h, double delta_t, 
               double alpha_metal, double alpha_agua);

    // Método para verificar si un punto está en la pared
    bool esPared(int r) const { return r >= radio_interior; }
    
    // Método para verificar si un punto está en la base metálica
    bool esBaseMetalica(int h, int r) const { return h == 0 && r >= radio_interior; }
    
    // Configurar distribución inicial de temperatura
    void setDistribucionTemperaturaInicial(const std::vector<std::vector<double>>& temp_init);
    
    // Configurar fuente de calor en la base
    void setFuenteCalorBase(const std::vector<std::vector<double>>& fuente);
    
    // Método para evolucionar la temperatura en el tiempo con coeficientes diferentes
    void evolucionarTemperatura();
    
    // Obtener temperatura en un punto específico (coordenadas cilíndricas)
    double getTemperatura(double h, double r, double theta = 0) const;
    
    // Obtener gradiente de temperatura en 3D
    std::vector<double> getGradienteTemperatura(double h, double r, double theta = 0) const;
    
    // Obtener coeficiente de difusividad en un punto
    double getAlpha(double h, double r) const;
    
    // Método para diferentes distribuciones de calor
    void configurarFuenteUniforme(double potencia);
    void configurarFuenteGaussiana(double potencia_centro, double sigma);
    void configurarFuenteAnular(double potencia, double radio_int, double radio_ext);
    
    // Conversión entre coordenadas
    std::vector<double> cartesianToCilindricas(double x, double y, double z) const;
    std::vector<double> cilindricasToCartesian(double h, double r, double theta) const;
    
    // Verificar si un punto está dentro del cilindro (solo interior)
    bool estaDentro(double x, double y, double z) const;
    
    // Verificar si un punto está en la pared
    bool estaEnPared(double x, double y, double z) const;
    
    // Inicializar malla de coeficientes
    void inicializarMallaCoeficientes();
    
    // Getters
    int getRadioInterior() const { return radio_interior; }
    int getRadioTotal() const { return radio_total; }
    int getGrosorPared() const { return grosor_pared; }
    int getAltura() const { return altura; }
    double getDT() const { return dt; }
    double getAlphaOlla() const { return alpha_olla; }
    double getAlphaInterior() const { return alpha_interior; }
    
    // Diagnóstico
    void diagnosticoParedesCompleto() const;
    void diagnosticoCoeficientes() const;
    
    // Visualización
    void imprimirTemperaturaBase() const;
};

#endif