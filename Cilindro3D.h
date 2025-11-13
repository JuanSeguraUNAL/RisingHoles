#ifndef CILINDRO_3D_EVOLUTIVO_H
#define CILINDRO_3D_EVOLUTIVO_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <cmath>

class Cilindro3DEvolutivo {
private:
    int radio, altura;
    double temperatura_ambiente;
    
    Eigen::MatrixXd temperatura;
    Eigen::MatrixXi es_frontera;
    
    // Parámetros térmicos del agua
    double alpha;          // Difusividad térmica [m²/s]
    double k;              // Conductividad térmica [W/m·K]
    double rho;            // Densidad [kg/m³]
    double cp;             // Calor específico [J/kg·K]
    
    // Fuentes de calor
    Eigen::MatrixXd fuente_calor;
    double potencia_total;
    
    // Parámetros para condiciones de frontera
    double dr_param, dz_param, dt_param;
    
public:
    Cilindro3DEvolutivo(int radio, int altura, double temp_ambiente = 25.0)
        : radio(radio), altura(altura), temperatura_ambiente(temp_ambiente),
          temperatura(Eigen::MatrixXd::Constant(altura, radio + 1, temp_ambiente)),
          es_frontera(Eigen::MatrixXi::Zero(altura, radio + 1)),
          alpha(1.43e-7), k(0.6), rho(1000.0), cp(4186.0),
          fuente_calor(Eigen::MatrixXd::Zero(altura, radio + 1)),
          potencia_total(0.0), dr_param(0.001), dz_param(0.001), dt_param(1.0) {
        
        inicializarFronteras();
    }
    
private:
    void inicializarFronteras() {
        // Fondo del cilindro
        es_frontera.row(0).setConstant(1);
        // Pared lateral
        es_frontera.col(radio).setConstant(1);
        // Eje central
        es_frontera.col(0).setConstant(2);
    }

public:
    // ===== CONFIGURACIÓN DE FUENTES DE CALOR =====
    
    void setFuenteCalorUniformeFondo(double potencia) {
        fuente_calor.row(0).setConstant(potencia / (radio + 1));
        potencia_total = potencia;
    }
    
    void setCalentamientoParedes(double potencia_pared) {
        for (int z = 0; z < altura; ++z) {
            fuente_calor(z, radio) = potencia_pared;
        }
    }

    // ===== RESOLUCIÓN DE LA ECUACIÓN DE CALOR =====
    
    void resolverEcuacionCalor(double dt, double dr, double dz) {
        // Guardar parámetros para condiciones de frontera
        dt_param = dt;
        dr_param = dr;
        dz_param = dz;
        
        Eigen::MatrixXd temperatura_nueva = temperatura;
        
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                // Saltar fronteras (condiciones Dirichlet)
                if (es_frontera(z, r) == 1) continue;
                
                // Calcular Laplaciano en coordenadas cilíndricas
                double laplaciano = calcularLaplacianoCilindrico(r, z, dr, dz);
                
                // Ecuación de calor: ∂T/∂t = α∇²T + Q/(ρ*cp)
                double termino_fuente = fuente_calor(z, r) / (rho * cp);
                double derivada_temporal = alpha * laplaciano + termino_fuente;
                
                // Método de Euler explícito
                temperatura_nueva(z, r) = temperatura(z, r) + dt * derivada_temporal;
                
                // Estabilidad: limitar cambios bruscos
                double cambio_max = 5.0; // °C por paso
                double cambio = temperatura_nueva(z, r) - temperatura(z, r);
                if (std::abs(cambio) > cambio_max) {
                    temperatura_nueva(z, r) = temperatura(z, r) + std::copysign(cambio_max, cambio);
                }
            }
        }
        
        // Aplicar condiciones de frontera
        aplicarCondicionesFrontera(temperatura_nueva, dr, dz, dt);
        
        temperatura = temperatura_nueva;
    }
    
private:
    double calcularLaplacianoCilindrico(int r, int z, double dr, double dz) const {
        double termino_r = 0.0, termino_z = 0.0;
        
        // Término radial: (1/r)∂/∂r(r ∂T/∂r)
        if (r == 0) {
            // En el eje, usar expansión de Taylor
            termino_r = 2.0 * (temperatura(z, 1) - temperatura(z, 0)) / (dr * dr);
        } else {
            double r_val = static_cast<double>(r) * dr;
            
            // Derivada segunda radial
            double d2T_dr2 = 0.0;
            if (r > 0 && r < radio) {
                d2T_dr2 = (temperatura(z, r+1) - 2.0 * temperatura(z, r) + temperatura(z, r-1)) / (dr * dr);
            }
            
            // Derivada primera radial (diferencias centradas)
            double dT_dr = 0.0;
            if (r > 0 && r < radio) {
                dT_dr = (temperatura(z, r+1) - temperatura(z, r-1)) / (2.0 * dr);
            } else if (r == radio) {
                dT_dr = (temperatura(z, r) - temperatura(z, r-1)) / dr;
            }
            
            termino_r = d2T_dr2 + (1.0 / r_val) * dT_dr;
        }
        
        // Término vertical: ∂²T/∂z²
        if (z > 0 && z < altura - 1) {
            termino_z = (temperatura(z+1, r) - 2.0 * temperatura(z, r) + temperatura(z-1, r)) / (dz * dz);
        }
        
        return termino_r + termino_z;
    }
    
    void aplicarCondicionesFrontera(Eigen::MatrixXd& T, double dr, double dz, double dt) const {
        // Fondo: temperatura fija (fuente de calor)
        for (int r = 0; r <= radio; ++r) {
            if (fuente_calor(0, r) > 0) {
                T(0, r) = std::max(temperatura_ambiente, T(0, r));
            }
        }
        
        // Paredes: condición mixta
        for (int z = 0; z < altura; ++z) {
            if (fuente_calor(z, radio) > 0) {
                T(z, radio) = std::max(temperatura_ambiente, T(z, radio));
            } else {
                // Aislada: derivada cero (usar valor del interior)
                T(z, radio) = T(z, radio-1);
            }
        }
        
        // Superficie superior: pérdida de calor al ambiente (simplificado)
        for (int r = 0; r <= radio; ++r) {
            double h = 10.0; // Coeficiente de transferencia de calor [W/m²·K]
            double area_celda = 2.0 * M_PI * r * dr; // Área superficial
            if (area_celda > 0) {
                double perdida = h * (T(altura-1, r) - temperatura_ambiente) * dt / (rho * cp * dz);
                T(altura-1, r) -= perdida;
                T(altura-1, r) = std::max(temperatura_ambiente, T(altura-1, r));
            }
        }
    }

public:
    // ===== ENFRIAMIENTO POR BURBUJAS =====
    
    void aplicarEnfriamientoBurbujas(const std::vector<std::tuple<double, double, double>>& posiciones_burbujas,
                                   double factor_enfriamiento = 0.1) {
        for (const auto& pos : posiciones_burbujas) {
            double x = std::get<0>(pos);
            double y = std::get<1>(pos);
            double z = std::get<2>(pos);
            
            // Convertir a cilíndricas
            double r = std::sqrt(x*x + y*y);
            int r_idx = static_cast<int>(r);
            int z_idx = static_cast<int>(z);
            
            if (r_idx >= 0 && r_idx <= radio && z_idx >= 0 && z_idx < altura) {
                // Enfriar localmente
                double delta_T = -factor_enfriamiento;
                temperatura(z_idx, r_idx) += delta_T;
                temperatura(z_idx, r_idx) = std::max(temperatura_ambiente, temperatura(z_idx, r_idx));
            }
        }
    }
    
    // ===== CÁLCULO DE CONVECCIÓN (simplificado) =====
    void calcularCampoConveccion(double dr, double dz) {
        // Para esta versión, usamos un campo de convección simplificado
        // basado en el gradiente de temperatura
        
        std::cout << "Calculando campo de convección simplificado..." << std::endl;
        
        // En una implementación completa, aquí resolverías las ecuaciones de Navier-Stokes
        // Por ahora, usamos un modelo simplificado donde la velocidad es proporcional
        // al gradiente de temperatura
        
        // Este método debería actualizar los campos v_r y v_z
        // Por simplicidad, asumimos un campo constante por ahora
    }
    
    // ===== MÉTODOS DE ACCESO =====
    double getTemperatura(double r, double z) const {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);
        if (r_idx < 0 || r_idx > radio || z_idx < 0 || z_idx >= altura) {
            return temperatura_ambiente;
        }
        return temperatura(z_idx, r_idx);
    }
    
    double getTemperaturaCartesianas(double x, double y, double z) const {
        double r = std::sqrt(x*x + y*y);
        return getTemperatura(r, z);
    }
    
    void setTemperatura(double r, double z, double temp) {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);
        if (r_idx >= 0 && r_idx <= radio && z_idx >= 0 && z_idx < altura) {
            temperatura(z_idx, r_idx) = temp;
        }
    }
    
    const Eigen::MatrixXd& getTemperaturaRef() const { return temperatura; }
    Eigen::MatrixXd& getTemperaturaRef() { return temperatura; }
    
    // ===== MÉTODOS DE CONVECCIÓN (para compatibilidad) =====
    double getVelocidadRadial(int r, int z) const {
        // Modelo simplificado: velocidad proporcional a gradiente de temperatura
        if (z > 0 && z < altura - 1 && r > 0 && r < radio) {
            double dT_dr = (getTemperatura(r+1, z) - getTemperatura(r-1, z)) / (2.0 * dr_param);
            return 0.01 * dT_dr; // Factor de escala
        }
        return 0.0;
    }
    
    double getVelocidadVertical(int r, int z) const {
        // Flotabilidad básica + efecto térmico
        if (z > 0 && z < altura - 1) {
            double T_local = getTemperatura(r, z);
            double flotabilidad = 0.05 * (T_local - temperatura_ambiente) / 50.0;
            return 0.1 + flotabilidad; // Velocidad base + flotabilidad
        }
        return 0.1;
    }
    
    // ===== INFORMACIÓN =====
    void imprimirEstadisticasTemperatura() const {
        std::cout << "  Temp min: " << temperatura.minCoeff() << "°C, max: " 
                  << temperatura.maxCoeff() << "°C, avg: " << temperatura.mean() << "°C" << std::endl;
    }
    
    void imprimirEstadisticasFlujo() const {
        std::cout << "  Campo de convección calculado (modelo simplificado)" << std::endl;
    }
    
    int getRadio() const { return radio; }
    int getAltura() const { return altura; }
    double getTemperaturaAmbiente() const { return temperatura_ambiente; }
    
    // ===== PARÁMETROS TÉRMICOS =====
    void setDifusividadTermica(double new_alpha) { alpha = new_alpha; }
    void setPotenciaTotal(double potencia) { potencia_total = potencia; }
};

#endif