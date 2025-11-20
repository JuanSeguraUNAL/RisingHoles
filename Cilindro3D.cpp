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
    
    // Inicializar temperaturas a 20°C (temperatura ambiente)
    temperatura.resize(altura, std::vector<double>(radio_total, 20.0));
    temperatura_nueva = temperatura;
    
    // Fuente de calor en toda la base (agua y metal)
    fuente_calor.resize(altura, std::vector<double>(radio_total, 0.0));
    
    // Inicializar coeficientes
    inicializarMallaCoeficientes();
}

void Cilindro3D::inicializarMallaCoeficientes() {
    alpha_malla.resize(altura, std::vector<double>(radio_total, alpha_interior));
    
    // Configurar coeficientes para las paredes metálicas
    for (int h = 0; h < altura; ++h) {
        for (int r = radio_interior; r < radio_total; ++r) {
            alpha_malla[h][r] = alpha_olla;
        }
    }
    
    // Base metálica también es metal
    for (int r = radio_interior; r < radio_total; ++r) {
        alpha_malla[0][r] = alpha_olla;
    }
}

void Cilindro3D::evolucionarTemperatura() {
    static int paso_count = 0;
    paso_count++;

    // ========== ECUACIÓN DE CALOR MEJORADA ==========
    
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            double alpha = alpha_malla[h][r];
            double laplaciano = 0.0;
            double termino_fuente = 0.0;

            // ========== FUENTE DE CALOR SOLO EN BASE ==========
            if (h == 0) {
                if (r < radio_interior) {
                    termino_fuente = fuente_calor[0][r];
                } else {
                    // Base metálica recibe calor de la zona de agua adyacente
                    termino_fuente = fuente_calor[0][radio_interior - 1] * 0.6;
                }
            }

            // ========== CORRECCIÓN CRÍTICA: LAPLACIANO EN COORDENADAS CILÍNDRICAS ==========
            // ∇²T = ∂²T/∂r² + (1/r)∂T/∂r + ∂²T/∂z²
            
            // TÉRMINO RADIAL: ∂²T/∂r² + (1/r)∂T/∂r
            if (r > 0 && r < radio_total - 1) {
                // Derivada segunda radial estándar
                double d2T_dr2 = temperatura[h][r-1] - 2.0 * temperatura[h][r] + temperatura[h][r+1];
                
                // Derivada primera radial (CORREGIDA: evitar división por cero)
                double dT_dr = (temperatura[h][r+1] - temperatura[h][r-1]) / 2.0;
                double termino_1sobre_r = (r > 0) ? dT_dr / (r + 1e-12) : 0.0; // +1e-12 evita división por cero
                
                laplaciano += d2T_dr2 + termino_1sobre_r;
            }
            else if (r == 0) {
                // EN EL EJE: usar desarrollo de Taylor para coordenadas cilíndricas
                // En r=0, por simetría: ∂T/∂r = 0 y ∇²T = 2 * ∂²T/∂r²
                laplaciano += 4.0 * (temperatura[h][1] - temperatura[h][0]);
            }
            else if (r == radio_total - 1) {
                // PARED EXTERIOR: condición de Neumann (pérdidas)
                double perdida = 0.05 * (20.0 - temperatura[h][r]);
                laplaciano += 2.0 * (temperatura[h][r-1] - temperatura[h][r] - perdida);
            }

            // ========== TÉRMINO VERTICAL MEJORADO ==========
            if (h > 0 && h < altura - 1) {
                // Conducción vertical estándar
                laplaciano += temperatura[h-1][r] - 2.0 * temperatura[h][r] + temperatura[h+1][r];
            }
            else if (h == 0 && altura > 1) {
                // BASE: incluir conducción hacia arriba + fuente
                laplaciano += temperatura[1][r] - temperatura[0][r];
            }
            else if (h == altura - 1 && altura > 1) {
                // SUPERFICIE: pérdidas por convección
                double perdida_superior = 0.02 * (20.0 - temperatura[h][r]);
                laplaciano += temperatura[h-1][r] - temperatura[h][r] - perdida_superior;
            }

            // ========== APLICAR ECUACIÓN DE CALOR ==========
            temperatura_nueva[h][r] = temperatura[h][r] + alpha * dt * laplaciano + termino_fuente * dt;
        }
    }

    // ========== CORRECCIÓN ESPECIAL PARA PAREDES METÁLICAS ==========
    // ¡ESTA ES LA CLAVE PARA EVITAR EL PATRÓN DE "PUNTOS"!
    
    // 1. CONDUCCIÓN MUY EFICIENTE EN PAREDES VERTICALES
    for (int r = radio_interior; r < radio_total; ++r) {
        // La base metálica calienta las paredes
        double temp_base_pared = temperatura_nueva[0][r];
        
        for (int h = 1; h < altura; ++h) {
            // Conducción vertical MUY eficiente en metal
            double factor_conduccion_metal = 0.7; // 70% de transferencia por paso
            
            // Suavizado para evitar patrones de puntos
            double temp_objetivo = temp_base_pared * (1.0 - 0.1 * h/altura) + 20.0 * (0.1 * h/altura);
            
            temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - factor_conduccion_metal) 
                                    + temp_objetivo * factor_conduccion_metal;
        }
    }
    
    // 2. CONDUCCIÓN HORIZONTAL EN PAREDES (evita patrones)
    for (int h = 0; h < altura; ++h) {
        for (int r = radio_interior + 1; r < radio_total - 1; ++r) {
            // Suavizado horizontal en paredes
            double promedio_horizontal = (temperatura_nueva[h][r-1] + temperatura_nueva[h][r+1]) / 2.0;
            double factor_suavizado = 0.3;
            temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - factor_suavizado) 
                                    + promedio_horizontal * factor_suavizado;
        }
    }

    // 3. ACOPLAMIENTO FUERTE AGUA-PARED
    for (int h = 1; h < altura; ++h) {
        if (radio_interior > 0) {
            double T_agua = temperatura_nueva[h][radio_interior - 1];
            double T_pared = temperatura_nueva[h][radio_interior];
            
            // Transferencia eficiente en la interfase
            double k_acoplamiento = 0.4;
            double transferencia = k_acoplamiento * (T_pared - T_agua);
            
            temperatura_nueva[h][radio_interior - 1] += transferencia;
            temperatura_nueva[h][radio_interior] -= transferencia * 0.5;
        }
    }

    // ========== MEZCLA CONVECTIVA SUAVIZADA ==========
    for (int r = 0; r < radio_interior; ++r) {
        for (int h = 1; h < altura - 1; ++h) {
            // Mezcla más suave para evitar patrones artificiales
            double mezcla = 0.08; // Reducido de 0.15 a 0.08
            
            // Promedio con vecinos inmediatos
            double temp_promedio = (temperatura_nueva[h-1][r] + temperatura_nueva[h+1][r] +
                                  temperatura_nueva[h][std::max(0, r-1)] + 
                                  temperatura_nueva[h][std::min(radio_interior-1, r+1)]) / 4.0;
            
            temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - mezcla) 
                                    + temp_promedio * mezcla;
        }
    }

    // ========== TRANSFERENCIA VERTICAL EN EL EJE (MEJORADA) ==========
    for (int h = 1; h < altura; ++h) {
        // El centro se calienta desde la base, pero de forma más física
        double temp_base_centro = temperatura_nueva[0][0];
        double factor_transferencia_vertical = 0.25 * (1.0 - (double)h/altura); // Disminuye con la altura
        
        temperatura_nueva[h][0] = temperatura_nueva[h][0] * (1.0 - factor_transferencia_vertical) 
                                + temp_base_centro * factor_transferencia_vertical;
    }

    // ========== ESTABILIDAD Y LÍMITES ==========
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            // Límites físicos
            temperatura_nueva[h][r] = std::max(20.0, std::min(120.0, temperatura_nueva[h][r]));
            
            // Suavizado final adicional para paredes (elimina patrones residuales)
            if (r >= radio_interior && h > 0 && h < altura - 1) {
                double suavizado_final = 0.1 * (temperatura_nueva[h-1][r] + temperatura_nueva[h+1][r] 
                                            - 2.0 * temperatura_nueva[h][r]);
                temperatura_nueva[h][r] += suavizado_final * 0.1;
            }
        }
    }

    // ========== ACTUALIZAR ==========
    temperatura = temperatura_nueva;

    // ========== DIAGNÓSTICO MEJORADO ==========
    if (paso_count % 25 == 0) {
        std::cout << "=== DIAGNÓSTICO PASO " << paso_count << " ===" << std::endl;
        
        // Temperaturas en el agua
        std::cout << "AGUA - Centro (r=0):" << std::endl;
        for (int h = 0; h < altura; h += altura/4) {
            std::cout << "  h=" << h << ": " << getTemperatura(h, 0) << "°C";
            if (h > 0) {
                double diff_base = getTemperatura(h, 0) - getTemperatura(0, 0);
                std::cout << " (diff_base=" << diff_base << "°C)";
            }
            std::cout << std::endl;
        }
        
        // Temperaturas en paredes
        std::cout << "PAREDES (r=" << radio_interior << "):" << std::endl;
        for (int h = 0; h < altura; h += altura/4) {
            std::cout << "  h=" << h << ": " << getTemperatura(h, radio_interior) << "°C";
            if (h > 0) {
                double diff_base = getTemperatura(h, radio_interior) - getTemperatura(0, radio_interior);
                std::cout << " (diff_base=" << diff_base << "°C)";
            }
            std::cout << std::endl;
        }
        
        // Verificar uniformidad de paredes
        std::cout << "UNIFORMIDAD PAREDES:" << std::endl;
        for (int h = altura/2; h <= altura/2; h += altura/4) {
            double min_temp = 1000.0, max_temp = -1000.0;
            for (int r = radio_interior; r < radio_total; ++r) {
                min_temp = std::min(min_temp, getTemperatura(h, r));
                max_temp = std::max(max_temp, getTemperatura(h, r));
            }
            std::cout << "  h=" << h << ": min=" << min_temp << "°C, max=" << max_temp 
                      << "°C, variación=" << (max_temp - min_temp) << "°C" << std::endl;
        }
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
    for (int r = 0; r < radio_interior; ++r) {
        fuente_calor[0][r] = potencia;
    }
    std::cout << "Fuente uniforme configurada: " << potencia << " K/s" << std::endl;
}

void Cilindro3D::configurarFuenteGaussiana(double potencia_centro, double sigma) {
    int centro = radio_interior / 2;
    
    for (int r = 0; r < radio_interior; ++r) {
        double distancia = std::abs(r - centro);
        fuente_calor[0][r] = potencia_centro * std::exp(-distancia*distancia/(2*sigma*sigma));
    }
    std::cout << "Fuente gaussiana configurada: centro=" << potencia_centro 
              << " K/s, sigma=" << sigma << std::endl;
}

void Cilindro3D::configurarFuenteAnular(double potencia, double radio_int, double radio_ext) {
    for (int r = 0; r < radio_interior; ++r) {
        if (r >= radio_int && r <= radio_ext) {
            fuente_calor[0][r] = potencia;
        }
    }
    std::cout << "Fuente anular configurada: " << potencia << " K/s, r=[" 
              << radio_int << "," << radio_ext << "]" << std::endl;
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
    
    // Verificar estabilidad numérica
    double Fo_agua = alpha_interior * dt;
    double Fo_metal = alpha_olla * dt;
    std::cout << "Números de Fourier efectivos:" << std::endl;
    std::cout << "  Agua: " << Fo_agua << std::endl;
    std::cout << "  Metal: " << Fo_metal << std::endl;
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
    if (fuente.size() == altura && fuente[0].size() == radio_total) {
        fuente_calor = fuente;
    }
}