#include "Cilindro3D.h"
#include <iostream>
#include <algorithm>
#include <cmath>

Cilindro3D::Cilindro3D(int radio_int, int grosor, int h, double delta_t, 
                       double alpha_metal, double alpha_agua) 
    : radio_interior(radio_int), grosor_pared(grosor), altura(h), 
      dt(delta_t), alpha_olla(alpha_metal), alpha_interior(alpha_agua) {
    
    radio_total = radio_interior + grosor_pared;
    
    std::cout << "Inicializando olla - Metal alpha=" << alpha_olla 
              << ", Agua alpha=" << alpha_interior << std::endl;
    
    // Inicializar temperaturas
    temperatura.resize(altura, std::vector<double>(radio_total, 20.0));
    temperatura_nueva = temperatura;
    fuente_calor.resize(1, std::vector<double>(radio_total, 0.0));
    
    // Inicializar coeficientes
    inicializarMallaCoeficientes();
}

void Cilindro3D::inicializarMallaCoeficientes() {
    alpha_malla.resize(altura, std::vector<double>(radio_total, alpha_interior));
    
    // TODO el metal (paredes laterales y base metálica)
    for (int h = 0; h < altura; ++h) {
        for (int r = radio_interior; r < radio_total; ++r) {
            alpha_malla[h][r] = alpha_olla;
        }
    }
}

void Cilindro3D::evolucionarTemperatura() {
    static int paso_count = 0;
    paso_count++;

    // ========== MODELO FÍSICO SIMPLIFICADO ==========
    
    // 1. APLICAR CALOR EN TODA LA BASE (agua Y metal)
    for (int r = 0; r < radio_total; ++r) {
        if (r < radio_interior) {
            // Base interior - fuente principal
            temperatura_nueva[0][r] = temperatura[0][r] + fuente_calor[0][r] * dt;
        } else {
            // Base metálica - también recibe calor directo (pero menos)
            double potencia_metal = fuente_calor[0][radio_interior-1] * 0.7 * dt;
            temperatura_nueva[0][r] = temperatura[0][r] + potencia_metal;
        }
        
        // Limitar temperatura base
        if (temperatura_nueva[0][r] > 120.0) temperatura_nueva[0][r] = 120.0;
    }

    // 2. PROPAGACIÓN VERTICAL EN PAREDES METÁLICAS (MUY RÁPIDA)
    for (int r = radio_interior; r < radio_total; ++r) {
        for (int h = 1; h < altura; ++h) {
            // El metal conduce el calor casi instantáneamente hacia arriba
            double temp_base = temperatura_nueva[0][r];
            double factor_atenuacion = 0.85; // Muy poco atenuación
            
            // La temperatura en la pared es casi igual a la base, con pequeña atenuación
            temperatura_nueva[h][r] = 20.0 + (temp_base - 20.0) * pow(factor_atenuacion, h);
            
            // Mínimo 40°C si la base está caliente
            if (temp_base > 80.0 && temperatura_nueva[h][r] < 40.0) {
                temperatura_nueva[h][r] = 40.0 + h * 2.0;
            }
        }
    }

    // 3. PROPAGACIÓN EN AGUA (más lenta, solo convección/difusión)
    for (int r = 0; r < radio_interior; ++r) {
        for (int h = 1; h < altura; ++h) {
            double temp_inferior = temperatura_nueva[h-1][r];
            double temp_actual = temperatura[h][r];
            
            // El agua se calienta más lentamente desde abajo
            double factor_agua = 0.3; // Mucho más lento que el metal
            
            temperatura_nueva[h][r] = temp_actual + factor_agua * (temp_inferior - temp_actual);
            
            // Transferencia desde paredes calientes
            if (r == radio_interior - 1) {
                double temp_pared = temperatura_nueva[h][radio_interior];
                double transferencia_pared = 0.4 * (temp_pared - temp_actual);
                temperatura_nueva[h][r] += transferencia_pared;
            }
            
            // Pérdidas en superficie
            if (h == altura - 1) {
                temperatura_nueva[h][r] -= 0.05 * (temp_actual - 20.0);
            }
        }
    }

    // 4. TRANSFERENCIA LATERAL AGUA ↔ PARED
    for (int h = 1; h < altura; ++h) {
        double temp_agua_borde = temperatura_nueva[h][radio_interior - 1];
        double temp_pared_interior = temperatura_nueva[h][radio_interior];
        
        // Equilibrar temperaturas en la interfase
        double promedio = (temp_agua_borde + temp_pared_interior) / 2.0;
        temperatura_nueva[h][radio_interior - 1] = promedio * 0.9;
        temperatura_nueva[h][radio_interior] = promedio * 1.1;
    }

    // 5. PÉRDIDAS EN PARED EXTERIOR
    for (int h = 0; h < altura; ++h) {
        double temp_pared = temperatura_nueva[h][radio_total - 1];
        double perdida = 0.08 * (temp_pared - 20.0);
        temperatura_nueva[h][radio_total - 1] -= perdida;
    }

    // 6. ESTABILIDAD
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            if (temperatura_nueva[h][r] < 20.0) temperatura_nueva[h][r] = 20.0;
            if (temperatura_nueva[h][r] > 150.0) temperatura_nueva[h][r] = 150.0;
        }
    }

    temperatura = temperatura_nueva;

    // Diagnóstico cada 20 pasos
    if (paso_count % 20 == 0) {
        std::cout << "Paso " << paso_count << " - Temperaturas clave:" << std::endl;
        std::cout << "  Base centro: " << temperatura[0][0] << "°C" << std::endl;
        std::cout << "  Base borde agua: " << temperatura[0][radio_interior-1] << "°C" << std::endl;
        std::cout << "  Base metal: " << temperatura[0][radio_interior] << "°C" << std::endl;
        std::cout << "  Pared h=10: " << temperatura[10][radio_interior] << "°C" << std::endl;
        std::cout << "  Pared h=20: " << temperatura[20][radio_interior] << "°C" << std::endl;
    }
}

// Los demás métodos se mantienen igual...
double Cilindro3D::getTemperatura(double h, double r, double theta) const {
    int h_idx = std::max(0, std::min(altura-1, static_cast<int>(h)));
    int r_idx = std::max(0, std::min(radio_total-1, static_cast<int>(r)));
    return temperatura[h_idx][r_idx];
}

std::vector<double> Cilindro3D::getGradienteTemperatura(double h, double r, double theta) const {
    std::vector<double> gradiente(3, 0.0);
    
    int h_idx = static_cast<int>(h);
    int r_idx = static_cast<int>(r);
    
    if (h_idx > 0 && h_idx < altura - 1) {
        gradiente[0] = (getTemperatura(h_idx+1, r_idx) - getTemperatura(h_idx-1, r_idx)) / 2.0;
    }
    if (r_idx > 0 && r_idx < radio_total - 1) {
        gradiente[1] = (getTemperatura(h_idx, r_idx+1) - getTemperatura(h_idx, r_idx-1)) / 2.0;
    }
    
    return gradiente;
}

double Cilindro3D::getAlpha(double h, double r) const {
    int h_idx = std::max(0, std::min(altura-1, static_cast<int>(h)));
    int r_idx = std::max(0, std::min(radio_total-1, static_cast<int>(r)));
    return alpha_malla[h_idx][r_idx];
}

std::vector<double> Cilindro3D::cartesianToCilindricas(double x, double y, double z) const {
    std::vector<double> cilindricas(3);
    cilindricas[0] = z;
    cilindricas[1] = std::sqrt(x*x + y*y);
    double theta = std::atan2(y, x);
    if (theta < 0) theta += 2.0 * M_PI;
    cilindricas[2] = theta;
    return cilindricas;
}

std::vector<double> Cilindro3D::cilindricasToCartesian(double h, double r, double theta) const {
    std::vector<double> cartesian(3);
    cartesian[0] = r * std::cos(theta);
    cartesian[1] = r * std::sin(theta);
    cartesian[2] = h;
    return cartesian;
}

bool Cilindro3D::estaDentro(double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y);
    return (z >= 0 && z <= altura && r <= radio_interior);
}

bool Cilindro3D::estaEnPared(double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y);
    return (z >= 0 && z <= altura && r > radio_interior && r <= radio_total);
}

void Cilindro3D::configurarFuenteUniforme(double potencia) {
    for (int r = 0; r < radio_total; ++r) {
        if (r < radio_interior) {
            fuente_calor[0][r] = potencia;
        }
    }
}

void Cilindro3D::configurarFuenteGaussiana(double potencia_centro, double sigma) {
    int centro = radio_interior / 2;
    
    for (int r = 0; r < radio_total; ++r) {
        if (r < radio_interior) {
            double distancia = std::abs(r - centro);
            fuente_calor[0][r] = potencia_centro * std::exp(-distancia*distancia/(2*sigma*sigma));
        }
    }
}

void Cilindro3D::configurarFuenteAnular(double potencia, double radio_int, double radio_ext) {
    for (int r = 0; r < radio_total; ++r) {
        if (r < radio_interior && r >= radio_int && r <= radio_ext) {
            fuente_calor[0][r] = potencia;
        }
    }
}

void Cilindro3D::diagnosticoParedesCompleto() const {
    std::cout << "=== DIAGNÓSTICO PAREDES ===" << std::endl;
    
    std::cout << "BASE (h=0):" << std::endl;
    for (int r = 0; r < radio_total; ++r) {
        std::string tipo = (r < radio_interior) ? "AGUA" : "METAL";
        std::cout << "r=" << r << "(" << tipo << "): " << temperatura[0][r] << "°C";
        if (r == radio_interior-1 || r == radio_interior) std::cout << " <-- INTERFASE";
        std::cout << std::endl;
    }
    
    std::cout << "PARED r=" << radio_interior << ":" << std::endl;
    for (int h = 0; h < altura; h += 5) {
        std::cout << "  h=" << h << ": " << temperatura[h][radio_interior] << "°C";
        if (h > 0) {
            double diff = temperatura[h][radio_interior] - temperatura[0][radio_interior];
            std::cout << " (diff_base=" << diff << "°C)";
        }
        std::cout << std::endl;
    }
}

void Cilindro3D::diagnosticoCoeficientes() const {
    std::cout << "=== DIAGNÓSTICO COEFICIENTES ===" << std::endl;
    
    int puntos_calientes = 0;
    for (int h = 0; h < altura; ++h) {
        for (int r = radio_interior; r < radio_total; ++r) {
            if (temperatura[h][r] >= 50.0) puntos_calientes++;
        }
    }
    
    std::cout << "Paredes >= 50°C: " << puntos_calientes << "/" << (altura * grosor_pared) << std::endl;
    std::cout << "Base metal: " << temperatura[0][radio_interior] << "°C" << std::endl;
    std::cout << "Pared media: " << temperatura[altura/2][radio_interior] << "°C" << std::endl;
    std::cout << "Pared superior: " << temperatura[altura-1][radio_interior] << "°C" << std::endl;
}

void Cilindro3D::imprimirTemperaturaBase() const {
    std::cout << "Temperatura base: ";
    for (int r = 0; r < radio_total; ++r) {
        std::string tipo = (r < radio_interior) ? "A" : "M";
        std::cout << "r" << tipo << r << "=" << temperatura[0][r] << " ";
    }
    std::cout << std::endl;
}

void Cilindro3D::setDistribucionTemperaturaInicial(const std::vector<std::vector<double>>& temp_init) {
    if (temp_init.size() == altura && temp_init[0].size() == radio_total) {
        temperatura = temp_init;
        temperatura_nueva = temp_init;
    }
}

void Cilindro3D::setFuenteCalorBase(const std::vector<std::vector<double>>& fuente) {
    if (fuente.size() == 1 && fuente[0].size() == radio_total) {
        fuente_calor = fuente;
    }
}