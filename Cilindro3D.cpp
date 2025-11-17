#include "Cilindro3D.h"
#include <iostream>
#include <algorithm>

Cilindro3D::Cilindro3D(int r, int h, double delta_t, double alpha_term) 
    : radio(r), altura(h), dt(delta_t), alpha(alpha_term) {
    
    temperatura.resize(altura, std::vector<double>(radio, 20.0));
    temperatura_nueva = temperatura;
    fuente_calor.resize(1, std::vector<double>(radio, 0.0));
}

void Cilindro3D::evolucionarTemperatura() {
    static int paso_count = 0;
    if (paso_count % 100 == 0) {
        std::cout << "Evolucionando temperatura - Paso " << paso_count << std::endl;
    }
    paso_count++;

    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio; ++r) {
            double laplaciano = 0.0;
            double termino_fuente = 0.0;
            double perdida = 0.0;

            // Fuente de calor solo en la base
            if (h == 0) {
                termino_fuente = fuente_calor[0][r] * dt;
            }

            // ========== TÉRMINO RADIAL ==========
            if (r > 0 && r < radio - 1) {
                // Puntos internos - conducción radial completa
                double deriv_radial = (temperatura[h][r+1] - temperatura[h][r-1]) / 2.0;
                laplaciano += (temperatura[h][r-1] - 2*temperatura[h][r] + temperatura[h][r+1])
                            + (1.0/(r+1e-6)) * deriv_radial;
            } 
            else if (r == 0 && radio > 1) {
                // Centro - aproximación especial por simetría
                laplaciano += 2.0 * (temperatura[h][r+1] - temperatura[h][r]);
            } 
            else if (r == radio - 1 && radio > 1) {
                // PARED LATERAL - MEJORADO: Permitir conducción desde el interior
                // Usar diferencia hacia atrás para la derivada radial
                laplaciano += (temperatura[h][r-1] - temperatura[h][r]); // Conducción desde el interior
                
                // Añadir un pequeño término de conducción desde el interior más profundo
                if (r - 2 >= 0) {
                    laplaciano += 0.3 * (temperatura[h][r-2] - temperatura[h][r-1]); // Conducción adicional
                }
            }

            // ========== TÉRMINO VERTICAL ==========
            if (h > 0 && h < altura - 1) {
                // Puntos internos - conducción vertical completa
                laplaciano += (temperatura[h-1][r] - 2*temperatura[h][r] + temperatura[h+1][r]);
            } 
            else if (h == 0 && altura > 1) {
                // Base - solo vecino superior
                laplaciano += (temperatura[h+1][r] - temperatura[h][r]);
            } 
            else if (h == altura - 1 && altura > 1) {
                // Tope superior - solo vecino inferior + pérdida
                laplaciano += (temperatura[h-1][r] - temperatura[h][r]);
                perdida += 0.02 * (temperatura[h][r] - 20.0) * dt; // Pérdida moderada
            }

            // ========== PÉRDIDAS TÉRMICAS ==========
            // Reducir significativamente las pérdidas en paredes laterales
            if (r == radio - 1) {
                // Pérdida MUY suave en paredes para permitir que se calienten
                perdida += 0.001 * (temperatura[h][r] - 20.0) * dt;
            }

            // ========== APLICAR ECUACIÓN ==========
            temperatura_nueva[h][r] = temperatura[h][r] + alpha * dt * laplaciano + termino_fuente - perdida;
            
            // Asegurar que la temperatura no sea negativa
            if (temperatura_nueva[h][r] < 20.0) {
                temperatura_nueva[h][r] = 20.0;
            }
        }
    }
    
    temperatura = temperatura_nueva;
}

double Cilindro3D::getTemperatura(double h, double r, double theta) const {
    int h_idx = std::max(0, std::min(altura-1, static_cast<int>(h)));
    int r_idx = std::max(0, std::min(radio-1, static_cast<int>(r)));
    return temperatura[h_idx][r_idx];
}

std::vector<double> Cilindro3D::getGradienteTemperatura(double h, double r, double theta) const {
    std::vector<double> gradiente(3, 0.0); // [dh, dr, dtheta]
    
    int h_idx = static_cast<int>(h);
    int r_idx = static_cast<int>(r);
    
    // Gradiente vertical
    if (h_idx > 0 && h_idx < altura - 1) {
        gradiente[0] = (getTemperatura(h_idx+1, r_idx) - getTemperatura(h_idx-1, r_idx)) / 2.0;
    }
    
    // Gradiente radial
    if (r_idx > 0 && r_idx < radio - 1) {
        gradiente[1] = (getTemperatura(h_idx, r_idx+1) - getTemperatura(h_idx, r_idx-1)) / 2.0;
    }
    
    // Por simetría axial, gradiente angular es cero
    gradiente[2] = 0.0;
    
    return gradiente;
}

std::vector<double> Cilindro3D::cartesianToCilindricas(double x, double y, double z) const {
    std::vector<double> cilindricas(3);
    cilindricas[0] = z; // altura
    cilindricas[1] = std::sqrt(x*x + y*y); // radio
    
    // Normalizar ángulo a [0, 2*pi]
    double theta = std::atan2(y, x);
    if (theta < 0) {
        theta += 2.0 * M_PI;
    }
    cilindricas[2] = theta;
    
    return cilindricas;
}

std::vector<double> Cilindro3D::cilindricasToCartesian(double h, double r, double theta) const {
    std::vector<double> cartesian(3);
    cartesian[0] = r * std::cos(theta); // x
    cartesian[1] = r * std::sin(theta); // y
    cartesian[2] = h;                   // z
    return cartesian;
}

bool Cilindro3D::estaDentro(double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y);
    return (z >= 0 && z <= altura && r <= radio);
}

void Cilindro3D::configurarFuenteUniforme(double potencia) {
    for (int r = 0; r < radio; ++r) {
        fuente_calor[0][r] = potencia;
    }
}

void Cilindro3D::configurarFuenteGaussiana(double potencia_centro, double sigma) {
    int centro = radio / 2;
    for (int r = 0; r < radio; ++r) {
        double distancia = std::abs(r - centro);
        fuente_calor[0][r] = potencia_centro * std::exp(-distancia * distancia / (2 * sigma * sigma));
    }
}

void Cilindro3D::imprimirTemperaturaBase() const {
    for (int r = 0; r < radio; ++r) {
        std::cout << temperatura[0][r] << " ";
    }
    std::cout << std::endl;
}