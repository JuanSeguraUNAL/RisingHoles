#ifndef MALLA_ARROZ_BURBUJAS_H
#define MALLA_ARROZ_BURBUJAS_H

#include "cilindro_3d.h"
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

class Burbuja {
public:
    double r, z;        // Posición en coordenadas cilíndricas
    double radio;       // Tamaño de la burbuja
    bool activa;        // Si la burbuja sigue activa
    int id;             // Identificador único
    
    Burbuja(double r, double z, double radio = 1.0, int id = 0) 
        : r(r), z(z), radio(radio), activa(true), id(id) {}
};

class MallaArrozBurbujas{
private:
    Cilindro3D& cilindro;  // Referencia al cilindro con convección
    Eigen::MatrixXd masa_arroz;  // Masa de arroz en cada celda (r,z)
    Eigen::MatrixXi ocupado_por_burbuja;  // 1 si hay burbuja, 0 si no
    
    std::vector<Burbuja> burbujas;
    int siguiente_id_burbuja;
    
    // Generador de números aleatorios
    std::mt19937 generador;
    std::uniform_real_distribution<double> distribucion;
    
    // Parámetros de simulación
    double flotabilidad_base;
    double factor_conveccion;
    double masa_inicial_arroz;
    double temperatura_umbral_burbujas;

public:
    MallaArrozBurbujas(Cilindro3D& cil, double masa_inicial = 1.0) 
        : cilindro(cil), 
          masa_arroz(Eigen::MatrixXd::Constant(cil.getAltura(), cil.getRadio() + 1, masa_inicial)),
          ocupado_por_burbuja(Eigen::MatrixXi::Zero(cil.getAltura(), cil.getRadio() + 1)),
          siguiente_id_burbuja(0),
          generador(std::random_device{}()),
          distribucion(0.0, 1.0),
          flotabilidad_base(0.1),
          factor_conveccion(1.0),
          masa_inicial_arroz(masa_inicial),
          temperatura_umbral_burbujas(50.0){}

    // Obtener masa en posición (r,z)
    double getMasaArroz(double r, double z) const {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);
        
        if (r_idx < 0 || r_idx > cilindro.getRadio() || z_idx < 0 || z_idx >= cilindro.getAltura()) {
            return 0.0;
        }
        return masa_arroz(z_idx, r_idx);
    }

    // Establecer masa en posición (r,z)
    void setMasaArroz(double r, double z, double masa) {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);
        
        if (r_idx >= 0 && r_idx <= cilindro.getRadio() && z_idx >= 0 && z_idx < cilindro.getAltura()) {
            masa_arroz(z_idx, r_idx) = masa;
        }
    }

    // Redistribuir masa cuando una burbuja se mueve
    void redistribuirMasa(double r_vieja, double z_vieja, double r_nueva, double z_nueva) {
        int r_old = static_cast<int>(r_vieja);
        int z_old = static_cast<int>(z_vieja);
        int r_new = static_cast<int>(r_nueva);
        int z_new = static_cast<int>(z_nueva);
        
        // Liberar la posición vieja (masa = 0 donde pasó la burbuja)
        if (r_old >= 0 && r_old <= cilindro.getRadio() && z_old >= 0 && z_old < cilindro.getAltura()) {
            masa_arroz(z_old, r_old) = 0.0;
            ocupado_por_burbuja(z_old, r_old) = 0;
        }
        
        // Marcar nueva posición como ocupada
        if (r_new >= 0 && r_new <= cilindro.getRadio() && z_new >= 0 && z_new < cilindro.getAltura()) {
            ocupado_por_burbuja(z_new, r_new) = 1;
        }
        
        // Redistribuir masa a celdas vecinas (conservación)
        redistribuirMasaVecinos(r_old, z_old);
    }

private:
    void redistribuirMasaVecinos(int r_centro, int z_centro) {
        double masa_total_redistribuir = masa_inicial_arroz;
        int vecinos_validos = 0;
        
        // Contar vecinos válidos
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dr = -1; dr <= 1; ++dr) {
                if (dr == 0 && dz == 0) continue; // Saltar el centro
                
                int r_vecino = r_centro + dr;
                int z_vecino = z_centro + dz;
                
                if (esPosicionValida(r_vecino, z_vecino) && ocupado_por_burbuja(z_vecino, r_vecino) == 0) {
                    vecinos_validos++;
                }
            }
        }
        
        if (vecinos_validos == 0) return;
        
        double masa_por_vecino = masa_total_redistribuir / vecinos_validos;
        
        // Distribuir masa a vecinos
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dr = -1; dr <= 1; ++dr) {
                if (dr == 0 && dz == 0) continue;
                
                int r_vecino = r_centro + dr;
                int z_vecino = z_centro + dz;
                
                if (esPosicionValida(r_vecino, z_vecino) && ocupado_por_burbuja(z_vecino, r_vecino) == 0) {
                    masa_arroz(z_vecino, r_vecino) += masa_por_vecino;
                }
            }
        }
    }

    bool esPosicionValida(int r, int z) const {
        return (r >= 0 && r <= cilindro.getRadio() && z >= 0 && z < cilindro.getAltura());
    }

public:
     void generarBurbujasSuperficie(int max_burbujas_por_paso = 5){
        if(max_burbujas_por_paso < = 0) return;
        int burbujas_generadas = 0;

        // 1. Construir distribucion de probabilidad sobre toda la superficie
        std::vector<std::pair<int, int>> posiciones_validas,
        std::vector<double> probabilidades;
        double suma_total = 0.0;

        // Considerar TODAS las posiciones de la superficie interna
        for (int r = 0; r <= cilindro.getRadio(); ++r) {
            for (int z = 0; z < cilindro.getAltura(); ++z) {
                // Solo considerar posiciones de superficie (fondo o paredes)
                if (!esPosicionSuperficie(r, z)) continue;
                
                // Verificar que no esté ocupada
                if (estaOcupadoPorBurbuja(r, z)) continue;
                
                double temperatura = cilindro.getTemperatura(r, z);
                double probabilidad = calcularProbabilidadTemperatura(temperatura);
                
                if (probabilidad > 0.0) {
                    posiciones_validas.emplace_back(r, z);
                    probabilidades.push_back(probabilidad);
                    suma_total += probabilidad;
                }
            }
        }

        // Si no hay posiciones validas, no genere burbujas
        if (posiciones_validas.empty() || suma_total <= 0.0) {
            return; // No hay posiciones válidas
        }

        // 2. Normalizar probabilidades
        for (double& prob : probabilidades) {
            prob /= suma_total;
        }

        // 3. Crear distribucion discreta para selccion aleatoria
        std::discrete_distribution<int> distribucion_burbujas(probabilidades.begin(), probabilidades.end());

        // 4. Generar burbujas usando la distribucion
        int burbujas_generadas = 0;
        int intentos_maximos = max_burbujas_por_paso * 10;
        int intentos = 0;
        while (burbujas_generadas < max_burbujas_por_paso && intentos < intentos_maximos){
            intentos++;

            // Seleccionar posicion aleatoria segun distribucion de probabilidad
            int indice = distribucion_burbujas(generador);
            int r = posiciones_validas[indice].first;
            int z = posiciones_validas[indice].second;

            // Verificar que todavía esté disponible (podría haber cambiado)
            if (estaOcupadoPorBurbuja(r, z)) {
                // Recalcular distribución si muchas posiciones se ocuparon
                if (intentos % 5 == 0) {
                    return;
                }
                continue;
            }

            // Generar la burbuja
            double temperatura = cilindro.getTemperatura(r, z);
            std::string tipo = obtenerTipoSuperficie(r, z);
            generarBurbujaEnPosicion(r, z, temperatura, tipo);
            burbujas_generadas++;
        }
    }

private:
    // Verificar si una posicion es de superficie (fondo o paredes de la olla)
    bool esPosicionSuperficie(int r, int z) const {
        // Fondo de la olla
        if (z == 0) return true;

        // Paredes laterales
        if (r == cilindro.getRadio()) return true;

        return false;
    }

    // Calcular probabilidad basada en temperatura
    double calcularProbabilidadTemperatura(double temperatura) const {
        if (temperatura < temperatura_umbral_burbujas){
            return 0.0;
        }

        // Funcion no lineal sensible a cambios en temperaturas altas
        double normalizada = (temperatura - temperatura_umbral_burbujas) / (120.0 - temperatura_umbral_burbujas); // 120°C como maximo

        normalizada = std::min(1.0, std::max(0.0, normalizada));

        // Funcion exponencial donde pequeñas diferencias en temperatura alta generan grandes diferencias en probabilidad
        double probabilidad = std::pow(normalizada, 2.0)

        // Suavizar con funcion sigmoide
        probabilidad = 1.0 / (1 + std::exp(-4.0 * (probabilidad - 0.5)));

        return probabilidad;
    }

    // Determinar tipo de superficie
    std::string obtenerTipoSuperficie(int, r, int z) const {
        if (z == 0) return "FONDO";
        if (r == cilindro.getRadio()) return "PARED";
        if (r == 0) return "EJE";
        return "INTERIOR";
    }

public:
    void moverBurbujas(double dt = 1.0){
        // Verificar coalescencia antes de mover burbujas
        verificarCoalescencia();

        for (auto& burbuja : burbujas){
            if(!burbuja.activa) continue;

            moverBurbuja(burbuja, dt);
        }

        // Eliminar burbujas inactivas
        burbujas.erase(
            std::remove_if(burbujas.begin(), burbujas.end(), 
                         [](const Burbuja& b) { return !b.activa; }),
            burbujas.end()
        );
    }

private:
    // Compara cada par de burbujas activas
    void verificarCoalescencia() {
        for (size_t i = 0; i < burbujas.size(); ++i) {
            if (!burbujas[i].activa) continue;
            
            for (size_t j = i + 1; j < burbujas.size(); ++j) {
                if (!burbujas[j].activa) continue;
                
                if (debenCoalescer(burbujas[i], burbujas[j])) {
                    coalescerBurbujas(burbujas[i], burbujas[j]);
                }
            }
        }
    }

    bool debenCoalescer(const Burbuja& b1, const Burbuja& b2) const {
        // Coalescer si están en la misma celda o muy cercanas
        double distancia_r = std::abs(b1.r - b2.r);
        double distancia_z = std::abs(b1.z - b2.z);
        
        // Misma celda o celdas adyacentes
        return (distancia_r <= 1.0 && distancia_z <= 1.0);
    }

    void coalescerBurbujas(Burbuja& b1, Burbuja& b2) {
        // Fusionar b2 en b1
        b1.radio = std::sqrt(b1.radio * b1.radio + b2.radio * b2.radio); // Conservar área
        // b1.r y b1.z se mantienen (o puedes promediar)
        
        // Marcar b2 como inactiva
        b2.activa = false;
        
        // Liberar la posición de b2
        int r2 = static_cast<int>(b2.r);
        int z2 = static_cast<int>(b2.z);
        if (esPosicionValida(r2, z2)) {
            ocupado_por_burbuja(z2, r2) = 0;
            // Redistribuir masa de la posición liberada
            redistribuirMasaVecinos(r2, z2);
        }
    }

    void moverBurbuja(Burbuja& burbuja, double dt){
        double r_old = burbuja.r;
        double z_old = burbuja.z;

        // Calcular probabilidades de movimiento
        std::vector<std::pair<double, double>> direcciones; // (dr, dz)
        std::vector<double> probabilidades;

        calcularProbabilidadesMovimiento(burbuja.r, burbuja.z, direcciones, probabilidades);

        if (probabilidades.empty()){
            burbuja.activa = false;
            return;
        }

        // Seleccionar movimiento basado en probabilidades
        double aleatorio = distribucion(generador);
        double acumulado = 0.0;
        int direccion_elegida = -1;
        
        for (size_t i = 0; i < probabilidades.size(); ++i) {
            acumulado += probabilidades[i];
            if (aleatorio <= acumulado) {
                direccion_elegida = i;
                break;
            }
        }

        if (direccion_elegida == -1) {
            direccion_elegida = probabilidades.size() - 1;
        }

        // Aplicar movimiento
        double dr = direcciones[direccion_elegida].first;
        double dz = direcciones[direccion_elegida].second;

        double r_new = burbuja.r + dr;
        double z_new = burbuja.z + dz;

        // Verificar que no baje de altura y esté dentro de límites
        if (z_new < burbuja.z) {
            z_new = burbuja.z; // No puede bajar
        }

        if (r_new < 0) r_new = 0;
        if (r_new > cilindro.getRadio()) r_new = cilindro.getRadio();
        if (z_new >= cilindro.getAltura()) {
            // Burbuja llegó a la superficie - desaparece
            redistribuirMasa(r_old, z_old, r_new, z_new);
            burbuja.activa = false;
            std::cout << "Burbuja " << burbuja.id << " llegó a la superficie" << std::endl;
            return;
        }

        // Verificar que la nueva posicion no este ocupada
        int r_new_idx = static_cast<int>(r_new);
        int z_new_idx = static_cast<int>(z_new);

        if (ocupado_por_burbuja(z_new_idx, r_new_idx) == 1) {
            // Posición ocupada - no se mueve
            return;
        }

        // Actualizar posición y redistribuir masa
        burbuja.r = r_new;
        burbuja.z = z_new;
        redistribuirMasa(r_old, z_old, r_new, z_new);
    }

    void calcularProbabilidadesMovimiento(double r, double z, 
                                        std::vector<std::pair<double, double>>& direcciones,
                                        std::vector<double>& probabilidades){
        direcciones.clear();
        probabilidades.clear();

        // Direcciones posibles: arriba, arriba-izquierda, arriba-derecha
        std::vector<std::pair<double, double>> movimientos_posibles = {
            {0.0, 1.0},    // Arriba
            {-1.0, 1.0},   // Arriba-izquierda  
            {1.0, 1.0}     // Arriba-derecha
        };

        // Obtener velocidades de convección
        double v_r = cilindro.getVelocidadRadial(static_cast<int>(r), static_cast<int>(z));
        double v_z = cilindro.getVelocidadVertical(static_cast<int>(r), static_cast<int>(z));

        for (const auto& movimiento : movimientos_posibles){
            double dr = movimiento.first;
            double dz = movimiento.second;
            
            double r_new = r + dr;
            double z_new = z + dz;

            // Verificar límites
            if (r_new < 0 || r_new > cilindro.getRadio() || z_new >= cilindro.getAltura()) {
                continue;
            }

            // Verificar que no esté ocupado
            int r_idx = static_cast<int>(r_new);
            int z_idx = static_cast<int>(z_new);
            if (ocupado_por_burbuja(z_idx, r_idx) == 1) {
                continue;
            }

            // Calcular probabilidad base (inversamente proporcional a masa)
            double masa_destino = getMasaArroz(r_new, z_new);
            double prob_base = 1.0 / (1.0 + masa_destino);

            // Efecto de convección (proporcional a velocidad en esa dirección)
            double efecto_conveccion = 1.0 + factor_conveccion * (v_r * dr + v_z * dz);

            // Efecto de flotabilidad (siempre favorece subir)
            double efecto_flotabilidad = (dz > 0) ? (1.0 + flotabilidad_base) : 1.0;

            double probabilidad_total = prob_base * efecto_conveccion * efecto_flotabilidad;

            direcciones.push_back(movimiento);
            probabilidades.push_back(probabilidad_total);
        }

        // Normalizar probabilidades
        double suma = 0.0;
        for (double prob : probabilidades) {
            suma += prob;
        }
        
        if (suma > 0) {
            for (double& prob : probabilidades) {
                prob /= suma;
            }
        }
    }

public:
    // Visualizacion y estadisticas

    void exportarEstadoCompleto(const std::string& nombre_archivo) const {
        std::ofstream archivo(nombre_archivo);
        archivo << "r z temperatura masa_arroz ocupado_burbuja v_r v_z" << std::endl;
        
        for (int z = 0; z < cilindro.getAltura(); ++z) {
            for (int r = 0; r <= cilindro.getRadio(); ++r) {
                archivo << r << " " << z << " " 
                       << cilindro.getTemperatura(r, z) << " "
                       << masa_arroz(z, r) << " "
                       << ocupado_por_burbuja(z, r) << " "
                       << cilindro.getVelocidadRadial(r, z) << " "
                       << cilindro.getVelocidadVertical(r, z) << std::endl;
            }
        }
        archivo.close();
    }

    void imprimirEstadisticas() const {
        std::cout << "=== ESTADÍSTICAS MALLA ARROZ ===" << std::endl;
        std::cout << "Burbujas activas: " << burbujas.size() << std::endl;
        std::cout << "Masa total arroz: " << masa_arroz.sum() << std::endl;
        std::cout << "Masa promedio por celda: " << masa_arroz.mean() << std::endl;
        std::cout << "Celdas ocupadas por burbujas: " << ocupado_por_burbuja.sum() << std::endl;
    }

    // Getters
    const std::vector<Burbuja>& getBurbujas() const { return burbujas; }
    const Eigen::MatrixXd& getMasaArrozRef() const { return masa_arroz; }
    const Eigen::MatrixXi& getOcupadoBurbujaRef() const { return ocupado_por_burbuja; }

    // Setters
    void setFlotabilidadBase(double flotabilidad) { flotabilidad_base = flotabilidad; }
    void setFactorConveccion(double factor) { factor_conveccion = factor; }
    void setTemperaturaUmbral(double temp) { temperatura_umbral_burbujas = temp; }

};

#endif