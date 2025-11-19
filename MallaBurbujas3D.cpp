#include "MallaBurbujas3D.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sys/stat.h>

// Constructor modificado
MallaBurbujas3D::MallaBurbujas3D(Cilindro3D& cil, int max_burbujas, const std::string& carpeta) 
    : cilindro(cil), max_burbujas_por_paso(max_burbujas), next_burbuja_id(0),
      carpeta_resultados(carpeta), paso_actual(0),
      factor_conveccion_metal(0.15),  // Mayor convección en metal
      factor_conveccion_agua(0.08),   // Menor convección en agua
      generator(std::random_device{}()), 
      distribution(0.0, 1.0),
      dist_angular(0.0, 2 * M_PI) {
    
    // Crear carpeta de resultados
    crearCarpetaResultados();
    
    // Inicializar malla de influencia
    int altura = cilindro.getAltura();
    int radio_total = cilindro.getRadioTotal();
    int sectores_angulares = 36;
    
    influencia_trayectorias.resize(altura, 
        std::vector<std::vector<double>>(radio_total, 
            std::vector<double>(sectores_angulares, 0.0)));
    
    std::cout << "Malla de burbujas inicializada con coeficientes diferentes" << std::endl;
    std::cout << "  - Factor convección metal: " << factor_conveccion_metal << std::endl;
    std::cout << "  - Factor convección agua: " << factor_conveccion_agua << std::endl;
}

// Método para crear la carpeta de resultados
void MallaBurbujas3D::crearCarpetaResultados() const {
    // Crear directorio principal
    system(("mkdir -p " + carpeta_resultados).c_str());
    system(("mkdir -p " + carpeta_resultados + "/influencia").c_str());
    system(("mkdir -p " + carpeta_resultados + "/burbujas").c_str());
    system(("mkdir -p " + carpeta_resultados + "/temperatura").c_str());
    system(("mkdir -p " + carpeta_resultados + "/geometria").c_str());
    system(("mkdir -p " + carpeta_resultados + "/coeficientes").c_str()); // Nueva carpeta
    std::cout << "Carpetas de resultados creadas en: " << carpeta_resultados << std::endl;
}

// NUEVO MÉTODO: Guardar coeficientes de difusividad
void MallaBurbujas3D::guardarCoeficientesDifusividad(const Cilindro3D& cilindro) const {
    std::string filename = carpeta_resultados + "/coeficientes/coeficientes_olla.txt";
    
    std::ofstream archivo(filename);
    if (!archivo.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return;
    }
    
    int altura = cilindro.getAltura();
    int radio_total = cilindro.getRadioTotal();
    int sectores = 36;

    archivo << "# Coeficientes de Difusividad Térmica" << std::endl;
    archivo << "# Dimensiones: Altura=" << altura << " RadioTotal=" << radio_total << std::endl;
    archivo << "# Formato: altura radio angulo coeficiente x y z tipo_material" << std::endl;
    archivo << "# tipo_material: 0=agua, 1=metal" << std::endl;
    
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            for (int theta_idx = 0; theta_idx < sectores; ++theta_idx) {
                double alpha = cilindro.getAlpha(h, r);
                double angulo = theta_idx * 2.0 * M_PI / sectores;
                int tipo_material = (r < cilindro.getRadioInterior()) ? 0 : 1;
                
                auto pos_cartesian = cilindro.cilindricasToCartesian(h, r, angulo);
                double x = pos_cartesian[0];
                double y = pos_cartesian[1];
                double z = pos_cartesian[2];
                
                archivo << h << " " << r << " " << angulo << " " << alpha
                        << " " << x << " " << y << " " << z << " " << tipo_material << std::endl;
            }
        }
    }
    
    archivo.close();
    std::cout << "Coeficientes de difusividad guardados: " << filename << std::endl;
}

// Método para guardar el mapa de influencia
void MallaBurbujas3D::guardarInfluenciaTrayectorias(int paso) const {
    std::string filename = carpeta_resultados + "/influencia/influencia_paso_" + 
                          std::to_string(paso) + ".txt";
    
    std::ofstream archivo(filename);
    if (!archivo.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return;
    }
    
    int altura = influencia_trayectorias.size();
    int radio = influencia_trayectorias[0].size();
    int sectores = influencia_trayectorias[0][0].size();
    
    // Encabezado con información de dimensiones
    archivo << "# Mapa de Influencia de Trayectorias - Paso: " << paso << std::endl;
    archivo << "# Dimensiones: Altura=" << altura << " Radio=" << radio << " Sectores=" << sectores << std::endl;
    archivo << "# Formato: altura radio angulo valor_influencia x y z tipo_material" << std::endl;
    
    // Guardar TODOS los datos de la malla
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio; ++r) {
            for (int theta_idx = 0; theta_idx < sectores; ++theta_idx) {
                double valor = influencia_trayectorias[h][r][theta_idx];
                double angulo = theta_idx * 2.0 * M_PI / sectores;
                int tipo_material = (r < cilindro.getRadioInterior()) ? 0 : 1;
                
                // Convertir a coordenadas cartesianas para visualización
                auto pos_cartesian = cilindro.cilindricasToCartesian(h, r, angulo);
                double x = pos_cartesian[0];
                double y = pos_cartesian[1];
                double z = pos_cartesian[2];
                
                archivo << h << " " << r << " " << angulo << " " << valor 
                        << " " << x << " " << y << " " << z << " " << tipo_material << std::endl;
            }
        }
    }
    
    archivo.close();
}

// Método adicional para guardar la geometría de la olla (solo una vez)
void MallaBurbujas3D::guardarGeometriaOlla(const Cilindro3D& cilindro) const {
    std::string filename = carpeta_resultados + "/geometria/geometria_olla.txt";
    
    std::ofstream archivo(filename);
    if (!archivo.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return;
    }
    
    int altura = cilindro.getAltura();
    int radio_interior = cilindro.getRadioInterior();
    int radio_total = cilindro.getRadioTotal();
    
    // Encabezado actualizado
    archivo << "# Geometría de la Olla" << std::endl;
    archivo << "# Dimensiones: Altura=" << altura << " RadioInterior=" << radio_interior 
            << " RadioTotal=" << radio_total << std::endl;
    archivo << "# Coeficientes: AlphaMetal=" << cilindro.getAlphaOlla() 
            << " AlphaAgua=" << cilindro.getAlphaInterior() << std::endl;
    archivo << "# Formato: tipo altura radio x y z" << std::endl;
    archivo << "# tipo: 0=base_interior, 1=base_pared, 2=pared_lateral, 3=interior" << std::endl;
    
    // Guardar puntos de la base INTERIOR
    for (int r = 0; r < radio_interior; ++r) {
        for (int theta_idx = 0; theta_idx < 36; ++theta_idx) {
            double angulo = theta_idx * 2.0 * M_PI / 36;
            auto pos = cilindro.cilindricasToCartesian(0, r, angulo);
            archivo << "0 " << 0 << " " << r << " " 
                    << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
    }
    
    // Guardar puntos de la base PARED
    for (int r = radio_interior; r < radio_total; ++r) {
        for (int theta_idx = 0; theta_idx < 36; ++theta_idx) {
            double angulo = theta_idx * 2.0 * M_PI / 36;
            auto pos = cilindro.cilindricasToCartesian(0, r, angulo);
            archivo << "1 " << 0 << " " << r << " " 
                    << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
    }
    
    // Guardar puntos de las paredes laterales
    for (int h = 0; h < altura; ++h) {
        for (int theta_idx = 0; theta_idx < 36; ++theta_idx) {
            double angulo = theta_idx * 2.0 * M_PI / 36;
            auto pos = cilindro.cilindricasToCartesian(h, radio_total-1, angulo);
            archivo << "2 " << h << " " << radio_total-1 << " " 
                    << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
    }
    
    // Guardar algunos puntos del interior para visualización
    for (int h = 1; h < altura-1; h += 2) {
        for (int r = 1; r < radio_interior-1; r += 2) {
            for (int theta_idx = 0; theta_idx < 12; ++theta_idx) {
                double angulo = theta_idx * 2.0 * M_PI / 12;
                auto pos = cilindro.cilindricasToCartesian(h, r, angulo);
                archivo << "3 " << h << " " << r << " " 
                        << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
            }
        }
    }
    
    archivo.close();
    std::cout << "Geometría de la olla guardada: " << filename << std::endl;
    
    // Guardar también los coeficientes
    guardarCoeficientesDifusividad(cilindro);
}

// Método para guardar estado de las burbujas (MODIFICADO)
void MallaBurbujas3D::guardarEstadoBurbujas(int paso) const {
    std::string filename = carpeta_resultados + "/burbujas/burbujas_paso_" + 
                          std::to_string(paso) + ".txt";
    
    std::ofstream archivo(filename);
    if (!archivo.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return;
    }
    
    // Encabezado mejorado
    archivo << "# Estado de Burbujas - Paso: " << paso << std::endl;
    archivo << "# Formato: id x y z tamano activa tipo_material alpha_local" << std::endl;
    
    // Guardar datos de cada burbuja
    for (const auto& burbuja : burbujas) {
        if (burbuja.isActiva()) {
            auto cil_coords = cilindro.cartesianToCilindricas(
                burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
            double alpha_local = cilindro.getAlpha(cil_coords[0], cil_coords[1]);
            std::string tipo_material = getTipoMaterial(burbuja);
            
            archivo << burbuja.getId() << " "
                    << burbuja.getPosX() << " "
                    << burbuja.getPosY() << " "
                    << burbuja.getPosZ() << " "
                    << burbuja.getTamano() << " "
                    << (burbuja.isActiva() ? 1 : 0) << " "
                    << tipo_material << " "
                    << alpha_local << std::endl;
        }
    }
    
    archivo.close();
}

// Método para guardar distribución de temperatura
void MallaBurbujas3D::guardarTemperatura(int paso, const Cilindro3D& cilindro) const {
    std::string filename = carpeta_resultados + "/temperatura/temperatura_paso_" + 
                          std::to_string(paso) + ".txt";
    
    std::ofstream archivo(filename);
    if (!archivo.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << std::endl;
        return;
    }
    
    int altura = cilindro.getAltura();
    int radio_total = cilindro.getRadioTotal();
    int sectores = 36;

    // Encabezado actualizado
    archivo << "# Distribución de Temperatura - Paso: " << paso << std::endl;
    archivo << "# Dimensiones: Altura=" << altura << " RadioTotal=" << radio_total << " Sectores=" << sectores << std::endl;
    archivo << "# Coeficientes: AlphaMetal=" << cilindro.getAlphaOlla() << " AlphaAgua=" << cilindro.getAlphaInterior() << std::endl;
    archivo << "# Formato: altura radio angulo temperatura alpha x y z tipo_material" << std::endl;
    
    // Guardar temperatura en TODOS los puntos de la malla 3D
    for (int h = 0; h < altura; ++h) {
        for (int r = 0; r < radio_total; ++r) {
            for (int theta_idx = 0; theta_idx < sectores; ++theta_idx) {
                double temp = cilindro.getTemperatura(h, r);
                double alpha = cilindro.getAlpha(h, r);
                double angulo = theta_idx * 2.0 * M_PI / sectores;
                int tipo_material = (r < cilindro.getRadioInterior()) ? 0 : 1;
                
                // Convertir a coordenadas cartesianas
                auto pos_cartesian = cilindro.cilindricasToCartesian(h, r, angulo);
                double x = pos_cartesian[0];
                double y = pos_cartesian[1];
                double z = pos_cartesian[2];
                
                archivo << h << " " << r << " " << angulo << " " << temp << " " << alpha
                        << " " << x << " " << y << " " << z << " " << tipo_material << std::endl;
            }
        }
    }
    
    archivo.close();
}

// MÉTODOS QUE FALTABAN - IMPLEMENTACIONES

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
    
    // Incrementar paso después de mover todas las burbujas
    incrementarPaso();
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

void MallaBurbujas3D::limpiarBurbujas() {
    for (auto& burbuja : burbujas) {
        if (!estaDentroCilindro(burbuja) || burbuja.getPosZ() >= cilindro.getAltura()) {
            burbuja.setActiva(false);
        }
    }
}

void MallaBurbujas3D::imprimirEstado() const {
    std::cout << "Burbujas activas: " << getNumBurbujasActivas() 
              << " / Total: " << burbujas.size() << std::endl;
}

int MallaBurbujas3D::getNumBurbujasActivas() const {
    return std::count_if(burbujas.begin(), burbujas.end(),
                        [](const Burbuja3D& b) { return b.isActiva(); });
}

// MÉTODOS AUXILIARES

double MallaBurbujas3D::calcularProbabilidadGeneracion(double h, double r, double theta) const {
    double temp = cilindro.getTemperatura(h, r, theta);
    double temp_umbral = 50.0;
    double temp_max = 100.0;
    
    if (temp < temp_umbral) return 0.0;
    return std::min(1.0, (temp - temp_umbral) / (temp_max - temp_umbral));
}

std::vector<double> MallaBurbujas3D::generarPosicionBase() {
    double theta = dist_angular(generator);
    // Solo generar en el radio interior (donde hay agua)
    double r = distribution(generator) * (cilindro.getRadioInterior() - 1);
    
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    double z = 0;
    
    return {x, y, z};
}

std::vector<double> MallaBurbujas3D::generarPosicionPared() {
    double theta = dist_angular(generator);
    double z = distribution(generator) * (cilindro.getAltura() - 1);
    // Generar en la pared interior (límite entre agua y metal)
    double r = cilindro.getRadioInterior() - 0.1;
    
    double x = r * std::cos(theta);
    double y = r * std::sin(theta);
    
    return {x, y, z};
}

bool MallaBurbujas3D::estaDentroCilindro(const Burbuja3D& burbuja) const {
    return cilindro.estaDentro(burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
}

// NUEVO MÉTODO: Aplicar efecto de coeficientes diferentes en el movimiento
void MallaBurbujas3D::aplicarEfectoCoeficientes(Burbuja3D& burbuja) {
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    
    double alpha_local = cilindro.getAlpha(cil_coords[0], cil_coords[1]);
    double alpha_agua = cilindro.getAlphaInterior();
    
    // Factor que modifica el movimiento basado en la difusividad local
    // Mayores coeficientes (metal) permiten movimientos más rápidos/definidos
    double factor_alpha = alpha_local / alpha_agua;
    
    // Aplicar efecto en la velocidad
    double vx = burbuja.getVelX();
    double vy = burbuja.getVelY();
    double vz = burbuja.getVelZ();
    
    // En metal, las burbujas tienden a moverse más rápido y con menos turbulencia
    if (factor_alpha > 1.5) { // Es metal
        // Movimiento más direccional y menos aleatorio
        vx *= (1.0 + 0.1 * (factor_alpha - 1.0));
        vy *= (1.0 + 0.1 * (factor_alpha - 1.0));
        vz *= (1.0 + 0.05 * (factor_alpha - 1.0)); // Menor efecto en dirección vertical
    }
    
    burbuja.setVelocidad(vx, vy, vz);
}

// Método modificado para calcular velocidad con efectos de coeficientes
void MallaBurbujas3D::calcularVelocidadBurbuja(Burbuja3D& burbuja) {
    double vx = 0, vy = 0, vz = 0;
    
    // 1. Flotabilidad (siempre hacia arriba)
    vz = 2.0 + burbuja.getTamano() * 0.5;
    
    // 2. Convección por gradiente de temperatura
    aplicarConveccion(burbuja);
    
    // 3. Influencia de trayectorias previas
    aplicarInfluenciaTrayectorias(burbuja);
    
    // 4. Efecto de coeficientes de difusividad (NUEVO)
    aplicarEfectoCoeficientes(burbuja);
    
    // 5. Movimiento browniano (modulado por el material)
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    double alpha_local = cilindro.getAlpha(cil_coords[0], cil_coords[1]);
    double factor_ruido = (alpha_local > cilindro.getAlphaInterior()) ? 0.05 : 0.1;
    
    vx += (distribution(generator) - 0.5) * factor_ruido;
    vy += (distribution(generator) - 0.5) * factor_ruido;
    
    burbuja.setVelocidad(vx, vy, vz);
}

// Método modificado para aplicar convección considerando coeficientes
void MallaBurbujas3D::aplicarConveccion(Burbuja3D& burbuja) {
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    
    auto gradiente = cilindro.getGradienteTemperatura(cil_coords[0], cil_coords[1], cil_coords[2]);
    double alpha_local = cilindro.getAlpha(cil_coords[0], cil_coords[1]);
    
    // Convertir gradiente cilíndrico a cartesiano
    double theta = cil_coords[2];
    double conv_x = gradiente[1] * std::cos(theta) - gradiente[2] * std::sin(theta);
    double conv_y = gradiente[1] * std::sin(theta) + gradiente[2] * std::cos(theta);
    double conv_z = gradiente[0];
    
    // Factor de convección dependiente del material
    double factor_conv = (alpha_local > cilindro.getAlphaInterior()) ? 
                        factor_conveccion_metal : factor_conveccion_agua;
    
    burbuja.setVelocidad(
        burbuja.getVelX() + conv_x * factor_conv,
        burbuja.getVelY() + conv_y * factor_conv,
        burbuja.getVelZ() + conv_z * factor_conv
    );
}

// MÉTODO QUE FALTABA: aplicarInfluenciaTrayectorias
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

// NUEVO MÉTODO: Obtener tipo de material para una burbuja
std::string MallaBurbujas3D::getTipoMaterial(const Burbuja3D& burbuja) const {
    auto cil_coords = cilindro.cartesianToCilindricas(
        burbuja.getPosX(), burbuja.getPosY(), burbuja.getPosZ());
    
    if (cil_coords[1] < cilindro.getRadioInterior()) {
        return "agua";
    } else {
        return "metal";
    }
}

// NUEVO MÉTODO: Configurar factores de convección
void MallaBurbujas3D::setFactoresConveccion(double factor_metal, double factor_agua) {
    factor_conveccion_metal = factor_metal;
    factor_conveccion_agua = factor_agua;
    std::cout << "Factores de convección actualizados:" << std::endl;
    std::cout << "  - Metal: " << factor_conveccion_metal << std::endl;
    std::cout << "  - Agua: " << factor_conveccion_agua << std::endl;
}