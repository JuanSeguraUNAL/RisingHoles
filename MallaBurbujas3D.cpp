#include "MallaBurbujas3D.h"
#include <cmath>
#include <iostream>
#include <algorithm>

MallaBurbujas3D::MallaBurbujas3D(Cilindro3D& cil, int max_burbujas) 
    : cilindro(cil), max_burbujas_por_paso(max_burbujas), next_burbuja_id(0),
      generator(std::random_device{}()), 
      distribution(0.0, 1.0),
      dist_angular(0.0, 2 * M_PI) {
    
    // Inicializar mapa de influencia 3D
    int altura = cilindro.getAltura();
    int radio = cilindro.getRadio();
    int sectores_angulares = 36; // 10 grados por sector
    
    influencia_trayectorias.resize(altura, 
        std::vector<std::vector<double>>(radio, 
            std::vector<double>(sectores_angulares, 0.0))); // CORREGIDO: sectores_angulares
}

void MallaBurbujas3D::generarBurbujas() {
    int burbujas_generadas = 0;
    
    // Generar en la base
    for (int i = 0; i < max_burbujas_por_paso / 2; ++i) {
        auto pos = generarPosicionBase();
        double x = pos[0], y = pos[1], z = pos[2];
        auto cil_coords = cilindro.cartesianToCilindricas(x, y, z);
        
        double prob = calcularProbabilidadGeneracion(cil_coords[0], cil_coords[1], cil_coords[2]);
        if (distribution(generator) < prob && burbujas_generadas < max_burbujas_por_paso) {
            burbujas.emplace_back(x, y, z, 1.0, next_burbuja_id++);
            burbujas_generadas++;
        }
    }
    
    // Generar en las paredes
    for (int i = 0; i < max_burbujas_por_paso / 2; ++i) {
        auto pos = generarPosicionPared();
        double x = pos[0], y = pos[1], z = pos[2];
        auto cil_coords = cilindro.cartesianToCilindricas(x, y, z);
        
        double prob = calcularProbabilidadGeneracion(cil_coords[0], cil_coords[1], cil_coords[2]);
        if (distribution(generator) < prob && burbujas_generadas < max_burbujas_por_paso) {
            burbujas.emplace_back(x, y, z, 1.0, next_burbuja_id++);
            burbujas_generadas++;
        }
    }
}

void MallaBurbujas3D::moverBurbujas() {
    double dt = cilindro.getDT();
    
    for (auto& burbuja : burbujas) {
        if (!burbuja.isActiva()) continue;
        
        // Calcular velocidad basada en múltiples factores
        calcularVelocidadBurbuja(burbuja);
        
        // Mover la burbuja
        burbuja.moverConVelocidad(dt);
        
        // Actualizar mapa de influencia
        auto cil_coords = cilindro.cartesianToCilindricas(
            burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
        
        int h_idx = static_cast<int>(cil_coords[0]);
        int r_idx = static_cast<int>(cil_coords[1]);
        int theta_idx = static_cast<int>(cil_coords[2] * influencia_trayectorias[0][0].size() / (2 * M_PI));
        
        if (h_idx >= 0 && h_idx < influencia_trayectorias.size() &&
            r_idx >= 0 && r_idx < influencia_trayectorias[0].size() &&
            theta_idx >= 0 && theta_idx < influencia_trayectorias[0][0].size()) {
            influencia_trayectorias[h_idx][r_idx][theta_idx] += 1.0;
        }
    }
}

void MallaBurbujas3D::calcularVelocidadBurbuja(Burbuja3D& burbuja) {
    double vx = 0, vy = 0, vz = 0;
    
    // 1. Flotabilidad (siempre hacia arriba)
    vz = 2.0 + burbuja.getTamano() * 0.5;
    
    // 2. Convección por gradiente de temperatura
    aplicarConveccion(burbuja);
    
    // 3. Influencia de trayectorias previas
    aplicarInfluenciaTrayectorias(burbuja);
    
    // 4. Movimiento browniano (pequeña componente aleatoria)
    vx += (distribution(generator) - 0.5) * 0.1;
    vy += (distribution(generator) - 0.5) * 0.1;
    
    burbuja.setVelocidad(vx, vy, vz);
}

void MallaBurbujas3D::aplicarConveccion(Burbuja3D& burbuja) {
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    
    auto gradiente = cilindro.getGradienteTemperatura(cil_coords[0], cil_coords[1], cil_coords[2]);
    
    // Convertir gradiente cilíndrico a cartesiano
    double theta = cil_coords[2];
    double conv_x = gradiente[1] * std::cos(theta) - gradiente[2] * std::sin(theta);
    double conv_y = gradiente[1] * std::sin(theta) + gradiente[2] * std::cos(theta);
    double conv_z = gradiente[0];
    
    burbuja.setVelocidad(
        burbuja.getVelX() + conv_x * 0.1,
        burbuja.getVelY() + conv_y * 0.1,
        burbuja.getVelZ() + conv_z * 0.1
    );
}

void MallaBurbujas3D::aplicarInfluenciaTrayectorias(Burbuja3D& burbuja) {
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    
    int h_idx = static_cast<int>(cil_coords[0]);
    int r_idx = static_cast<int>(cil_coords[1]);
    int theta_idx = static_cast<int>(cil_coords[2] * influencia_trayectorias[0][0].size() / (2 * M_PI));
    
    // Buscar trayectorias cercanas en el mapa de influencia
    double influencia_total = 0.0;
    double influencia_x = 0.0;
    double influencia_y = 0.0;
    double influencia_z = 0.0;
    
    int radio_busqueda = 2;
    int sectores = influencia_trayectorias[0][0].size();
    
    for (int dh = -radio_busqueda; dh <= radio_busqueda; ++dh) {
        for (int dr = -radio_busqueda; dr <= radio_busqueda; ++dr) {
            for (int dtheta = -1; dtheta <= 1; ++dtheta) {
                int h_actual = h_idx + dh;
                int r_actual = r_idx + dr;
                int theta_actual = (theta_idx + dtheta + sectores) % sectores;
                
                if (h_actual >= 0 && h_actual < influencia_trayectorias.size() &&
                    r_actual >= 0 && r_actual < influencia_trayectorias[0].size()) {
                    
                    double influencia = influencia_trayectorias[h_actual][r_actual][theta_actual];
                    if (influencia > 0) {
                        // Convertir a coordenadas cartesianas para la dirección
                        auto pos_influencia = cilindro.cilindricasToCartesian(
                            h_actual, r_actual, theta_actual * 2 * M_PI / sectores);
                        
                        double dx = pos_influencia[0] - burbuja.getPosX();
                        double dy = pos_influencia[1] - burbuja.getPosY();
                        double dz = pos_influencia[2] - burbuja.getPosZ();
                        
                        double dist = std::sqrt(dx*dx + dy*dy + dz*dz + 1e-6);
                        
                        influencia_x += influencia * dx / dist;
                        influencia_y += influencia * dy / dist;
                        influencia_z += influencia * dz / dist;
                        influencia_total += influencia;
                    }
                }
            }
        }
    }
    
    // Aplicar influencia si hay trayectorias cercanas
    if (influencia_total > 0) {
        double factor_influencia = 0.5; // Factor de seguimiento de trayectorias
        burbuja.setVelocidad(
            burbuja.getVelX() + influencia_x / influencia_total * factor_influencia,
            burbuja.getVelY() + influencia_y / influencia_total * factor_influencia,
            burbuja.getVelZ() + influencia_z / influencia_total * factor_influencia
        );
    }
}

void MallaBurbujas3D::verificarCoalescencia() {
    double radio_coalescencia = 0.5;
    
    for (size_t i = 0; i < burbujas.size(); ++i) {
        if (!burbujas[i].isActiva()) continue;
        
        for (size_t j = i + 1; j < burbujas.size(); ++j) {
            if (!burbujas[j].isActiva()) continue;
            
            if (burbujas[i].distancia(burbujas[j]) < radio_coalescencia * 
                (burbujas[i].getTamano() + burbujas[j].getTamano())) {
                
                burbujas[i].coalescer(burbujas[j]);
                burbujas[j].setActiva(false);
            }
        }
    }
}

std::vector<double> MallaBurbujas3D::generarPosicionBase() {
    double theta = dist_angular(generator);
    double r = distribution(generator) * (cilindro.getRadio() - 1);
    
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = 0;
    
    return {x, y, z};
}

std::vector<double> MallaBurbujas3D::generarPosicionPared() {
    double theta = dist_angular(generator);
    double z = distribution(generator) * (cilindro.getAltura() - 1);
    double r = cilindro.getRadio() - 0.1;
    
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    
    return {x, y, z};
}

double MallaBurbujas3D::calcularProbabilidadGeneracion(double h, double r, double theta) const {
    double temp = cilindro.getTemperatura(h, r, theta);
    double temp_umbral = 50.0;
    double temp_max = 100.0;
    
    if (temp < temp_umbral) return 0.0;
    return std::min(1.0, (temp - temp_umbral) / (temp_max - temp_umbral));
}

bool MallaBurbujas3D::estaDentroCilindro(const Burbuja3D& burbuja) const {
    return cilindro.estaDentro(burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
}

void MallaBurbujas3D::limpiarBurbujas() {
    for (auto& burbuja : burbujas) {
        if (!estaDentroCilindro(burbuja) || burbuja.getPosZ() >= cilindro.getAltura()) {
            burbuja.setActiva(false);
        }
    }
}

int MallaBurbujas3D::getNumBurbujasActivas() const {
    return std::count_if(burbujas.begin(), burbujas.end(),
                        [](const Burbuja3D& b) { return b.isActiva(); });
}

void MallaBurbujas3D::imprimirEstado() const {
    std::cout << "Burbujas activas: " << getNumBurbujasActivas() 
              << " / Total: " << burbujas.size() << std::endl;
}