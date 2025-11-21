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
    
    // 1. CONDUCCIÓN MUY EFICIENTE EN PAREDES VERTICALES
    for (int r = radio_interior; r < radio_total; ++r) {
        // La base metálica calienta las paredes
        double temp_base_pared = temperatura_nueva[0][r];
        
        for (int h = 1; h < altura; ++h) {
            // Conducción vertical MUY eficiente en metal - MÁS FUERTE
            double factor_conduccion_metal = 0.8; // Aumentado de 0.7 a 0.8 (80% de transferencia)
            
            // Suavizado para evitar patrones de puntos
            double temp_objetivo = temp_base_pared * (1.0 - 0.08 * h/altura) + 20.0 * (0.08 * h/altura); // Menor atenuación
            
            temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - factor_conduccion_metal) 
                                    + temp_objetivo * factor_conduccion_metal;
            
            // GARANTIZAR que las paredes estén más calientes que el agua adyacente
            if (r == radio_interior && h > 0) {
                double temp_agua_adyacente = temperatura_nueva[h][radio_interior - 1];
                if (temperatura_nueva[h][r] < temp_agua_adyacente + 3.0) {
                    temperatura_nueva[h][r] = temp_agua_adyacente + 3.0; // Mínimo 3°C más caliente
                }
            }
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

    // 3. ACOPLAMIENTO AGUA-PARED BALANCEADO (las paredes calientan el agua, no al revés)
    for (int h = 1; h < altura; ++h) {
        if (radio_interior > 0) {
            double T_agua = temperatura_nueva[h][radio_interior - 1];
            double T_pared = temperatura_nueva[h][radio_interior];
            
            // Solo transferir calor si la pared está más caliente
            if (T_pared > T_agua) {
                double k_acoplamiento = 0.3; // Reducido de 0.4 a 0.3
                double transferencia = k_acoplamiento * (T_pared - T_agua);
                
                temperatura_nueva[h][radio_interior - 1] += transferencia;
                // La pared se enfría menos que antes
                temperatura_nueva[h][radio_interior] -= transferencia * 0.3; // Reducido de 0.5 a 0.3
            }
        }
    }

    // ========== CORRECCIÓN PARA EL AGUA (MÁS CONSERVADORA) ==========
    
    // 1. CONDUCCIÓN VERTICAL DESDE LA BASE (MÁS SUAVE)
    for (int r = 0; r < radio_interior; ++r) {
        double temp_base = temperatura_nueva[0][r];
        
        for (int h = 1; h < altura; ++h) {
            // Factor más conservador para evitar que el agua se caliente demasiado
            double factor_vertical = 0.25 * (1.0 - 0.8 * (double)h/altura); // Reducido de 0.35 a 0.25
            
            // Solo aplicar si la base está significativamente más caliente
            if (temp_base > temperatura_nueva[h][r] + 8.0) { // Aumentado el umbral de 5.0 a 8.0
                temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - factor_vertical) 
                                        + temp_base * factor_vertical;
            }
        }
    }

    // 2. CONDUCCIÓN RADIAL MEJORADA EN AGUA (MÁS SUAVE)
    for (int h = 1; h < altura - 1; ++h) {
        for (int r = 1; r < radio_interior - 1; ++r) {
            // Calcular gradientes radiales reales
            double gradiente_izq = temperatura_nueva[h][r-1] - temperatura_nueva[h][r];
            double gradiente_der = temperatura_nueva[h][r+1] - temperatura_nueva[h][r];
            
            // Aplicar conducción radial más suave
            double factor_conduccion_radial = 0.1; // Reducido de 0.15 a 0.1
            
            if (std::abs(gradiente_izq) > 3.0) { // Aumentado el umbral de 2.0 a 3.0
                temperatura_nueva[h][r] += factor_conduccion_radial * gradiente_izq;
                temperatura_nueva[h][r-1] -= factor_conduccion_radial * gradiente_izq * 0.5;
            }
            
            if (std::abs(gradiente_der) > 3.0) { // Aumentado el umbral de 2.0 a 3.0
                temperatura_nueva[h][r] += factor_conduccion_radial * gradiente_der;
                temperatura_nueva[h][r+1] -= factor_conduccion_radial * gradiente_der * 0.5;
            }
        }
    }

    // 3. MEZCLA CONVECTIVA MÁS SUAVE
    for (int h = 1; h < altura - 1; ++h) {
        for (int r = 1; r < radio_interior - 1; ++r) {
            // La convección es más suave
            double diferencia_max = 0.0;
            diferencia_max = std::max(diferencia_max, std::abs(temperatura_nueva[h][r] - temperatura_nueva[h-1][r]));
            diferencia_max = std::max(diferencia_max, std::abs(temperatura_nueva[h][r] - temperatura_nueva[h+1][r]));
            diferencia_max = std::max(diferencia_max, std::abs(temperatura_nueva[h][r] - temperatura_nueva[h][r-1]));
            diferencia_max = std::max(diferencia_max, std::abs(temperatura_nueva[h][r] - temperatura_nueva[h][r+1]));
            
            // Convección más conservadora
            double fuerza_conveccion = 0.03 + 0.07 * (diferencia_max / 40.0); // Reducida
            fuerza_conveccion = std::min(0.15, fuerza_conveccion); // Reducido el máximo
            
            // Mezcla con vecinos ponderada
            double temp_promedio = (
                temperatura_nueva[h-1][r] + temperatura_nueva[h+1][r] +
                temperatura_nueva[h][r-1] + temperatura_nueva[h][r+1]
            ) / 4.0;
            
            temperatura_nueva[h][r] = temperatura_nueva[h][r] * (1.0 - fuerza_conveccion) 
                                    + temp_promedio * fuerza_conveccion;
        }
    }

    // 4. TRANSFERENCIA DESDE PAREDES CALIENTES HACIA EL AGUA (MÁS CONTROLADA)
    for (int h = 1; h < altura; ++h) {
        for (int r = radio_interior - 3; r < radio_interior; ++r) {
            if (r >= 0) {
                double distancia_a_pared = radio_interior - r;
                double temp_pared = temperatura_nueva[h][radio_interior];
                double temp_agua = temperatura_nueva[h][r];
                
                // Solo transferir si la pared está más caliente
                if (temp_pared > temp_agua) {
                    // Transferencia más controlada
                    double factor_penetracion = 0.3 / (distancia_a_pared + 0.5); // Reducido de 0.5 a 0.3
                    double transferencia = factor_penetracion * (temp_pared - temp_agua) * dt;
                    
                    temperatura_nueva[h][r] += transferencia;
                }
            }
        }
    }

    // ========== GARANTIZAR JERARQUÍA FÍSICA: PAREDES > AGUA ==========
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_interior; ++r) {
            // Para puntos de agua cerca de la pared, asegurar que estén más fríos
            if (r >= radio_interior - 2) {
                double temp_pared_adyacente = temperatura_nueva[h][radio_interior];
                if (temperatura_nueva[h][r] > temp_pared_adyacente - 2.0) {
                    temperatura_nueva[h][r] = temp_pared_adyacente - 2.0; // Mínimo 2°C más frío
                }
            }
        }
    }

    // ========== ESTABILIDAD Y LÍMITES ==========
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            // Límites físicos
            temperatura_nueva[h][r] = std::max(20.0, std::min(120.0, temperatura_nueva[h][r]));
        }
    }

    // ========== ACTUALIZAR ==========
    temperatura = temperatura_nueva;

    // ========== DIAGNÓSTICO MEJORADO ==========
    if (paso_count % 25 == 0) {
        std::cout << "=== DIAGNÓSTICO BALANCEADO PASO " << paso_count << " ===" << std::endl;
        
        // Comparación directa agua vs pared
        std::cout << "COMPARACIÓN AGUA vs PARED:" << std::endl;
        for (int h = 0; h < altura; h += altura/3) {
            double temp_agua_borde = getTemperatura(h, radio_interior - 1);
            double temp_pared = getTemperatura(h, radio_interior);
            double diferencia = temp_pared - temp_agua_borde;
            std::cout << "  h=" << h << ": Agua=" << temp_agua_borde << "°C, Pared=" << temp_pared 
                      << "°C, Diferencia=" << diferencia << "°C";
            if (diferencia < 0) {
                std::cout << " ⚠️ PARED MÁS FRÍA!"; // Esto no debería pasar
            }
            std::cout << std::endl;
        }
        
        // Distribución en agua
        std::cout << "DISTRIBUCIÓN EN AGUA (h=" << altura/2 << "):" << std::endl;
        for (int r = 0; r < radio_interior; r += 2) {
            std::cout << "  r=" << r << ": " << getTemperatura(altura/2, r) << "°C" << std::endl;
        }
        
        // Eficiencias
        double temp_base_agua = getTemperatura(0, radio_interior - 1);
        double temp_base_pared = getTemperatura(0, radio_interior);
        double temp_medio_agua = getTemperatura(altura/2, radio_interior - 1);
        double temp_medio_pared = getTemperatura(altura/2, radio_interior);
        
        std::cout << "EFICIENCIAS:" << std::endl;
        std::cout << "  Conducción vertical agua: " << (temp_medio_agua / temp_base_agua * 100) << "%" << std::endl;
        std::cout << "  Conducción vertical pared: " << (temp_medio_pared / temp_base_pared * 100) << "%" << std::endl;
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