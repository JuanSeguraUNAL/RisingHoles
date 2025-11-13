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
    // Ecuación de calor en coordenadas cilíndricas (simplificada por simetría axial)
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio; ++r) {
            if (h == 0) {
                // Base con fuente de calor
                double laplaciano = 0.0;
                if (r > 0 && r < radio - 1) {
                    // Término radial: (1/r) * d/dr(r * dT/dr)
                    double deriv_radial = (temperatura[h][r+1] - temperatura[h][r-1]) / 2.0;
                    if (r > 0) {
                        laplaciano = (temperatura[h][r-1] - 2*temperatura[h][r] + temperatura[h][r+1]) +
                                   (1.0/(r+1e-6)) * deriv_radial;
                    } else {
                        laplaciano = (temperatura[h][r+1] - temperatura[h][r]);
                    }
                    // Término vertical
                    laplaciano += (temperatura[h+1][r] - temperatura[h][r]);
                }
                temperatura_nueva[h][r] = temperatura[h][r] + alpha * dt * laplaciano + fuente_calor[0][r] * dt;
            }
            else if (r == 0 || r == radio - 1 || h == altura - 1) {
                // Condiciones de frontera
                temperatura_nueva[h][r] = temperatura[h][r] - 0.01 * (temperatura[h][r] - 20.0) * dt;
            }
            else {
                // Puntos internos
                double deriv_radial = (temperatura[h][r+1] - temperatura[h][r-1]) / 2.0;
                double laplaciano = (temperatura[h][r-1] - 2*temperatura[h][r] + temperatura[h][r+1]) +
                                  (1.0/(r+1e-6)) * deriv_radial +
                                  (temperatura[h-1][r] - 2*temperatura[h][r] + temperatura[h+1][r]);
                
                temperatura_nueva[h][r] = temperatura[h][r] + alpha * dt * laplaciano;
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
    cilindricas[2] = std::atan2(y, x); // ángulo
    return cilindricas;
}

std::vector<double> Cilindro3D::cilindricasToCartesian(double h, double r, double theta) const {
    std::vector<double> cartesian(3);
    cartesian[0] = r * std::cos(theta); // x
    cartesian[1] = r * std::sin(theta); // y
    cartesian[2] = h; // z
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