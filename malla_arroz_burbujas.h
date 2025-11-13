#ifndef MALLA_ARROZ_BURBUJAS_3D_SIMPLE_H
#define MALLA_ARROZ_BURBUJAS_3D_SIMPLE_H

#include "Cilindro3D.h"
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <algorithm>

class Burbuja3D {
public:
    double x, y, z;
    double radio;
    bool activa;
    int id;
    
    Burbuja3D(double x, double y, double z, double radio = 1.0, int id = 0) 
        : x(x), y(y), z(z), radio(radio), activa(true), id(id) {}
    
    void aCilindricas(double& r, double& theta, double& z_out) const {
        r = std::sqrt(x*x + y*y);
        theta = std::atan2(y, x);
        z_out = z;
    }
    
    double distancia(const Burbuja3D& otra) const {
        double dx = x - otra.x;
        double dy = y - otra.y;
        double dz = z - otra.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

class MallaArrozBurbujas3DSimple {
private:
    Cilindro3DEvolutivo& cilindro;
    
    // Malla 3D simple usando vectores anidados (sin Eigen::Tensor)
    int nx, ny, nz;
    double dx, dy, dz;
    std::vector<std::vector<std::vector<double>>> masa_arroz;
    std::vector<std::vector<std::vector<int>>> ocupado_por_burbuja;
    
    std::vector<Burbuja3D> burbujas;
    int siguiente_id_burbuja;
    
    std::mt19937 generador;
    std::uniform_real_distribution<double> distribucion;
    
    double flotabilidad_base;
    double factor_conveccion;
    double masa_inicial_arroz;
    double temperatura_umbral_burbujas;
    
    double x_min, x_max, y_min, y_max, z_min, z_max;

public:
    MallaArrozBurbujas3DSimple(Cilindro3DEvolutivo& cil, int res_x = 30, int res_y = 30, int res_z = 20, 
                              double masa_inicial = 1.0) 
        : cilindro(cil), nx(res_x), ny(res_y), nz(res_z),
          siguiente_id_burbuja(0),
          generador(std::random_device{}()),
          distribucion(0.0, 1.0),
          flotabilidad_base(0.1),
          factor_conveccion(1.0),
          masa_inicial_arroz(masa_inicial),
          temperatura_umbral_burbujas(50.0) {
        
        // Calcular límites
        double radio = cilindro.getRadio();
        double altura = cilindro.getAltura();
        
        x_min = -radio; x_max = radio;
        y_min = -radio; y_max = radio; 
        z_min = 0; z_max = altura;
        
        dx = (x_max - x_min) / (nx - 1);
        dy = (y_max - y_min) / (ny - 1);
        dz = (z_max - z_min) / (nz - 1);
        
        // Inicializar malla 3D
        inicializarMalla3D(masa_inicial);
    }
    
private:
    void inicializarMalla3D(double masa_inicial) {
        masa_arroz.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, masa_inicial)));
        ocupado_por_burbuja.resize(nx, std::vector<std::vector<int>>(ny, std::vector<int>(nz, 0)));
    }
    
    void aIndices(double x, double y, double z, int& i, int& j, int& k) const {
        i = static_cast<int>((x - x_min) / dx);
        j = static_cast<int>((y - y_min) / dy);
        k = static_cast<int>((z - z_min) / dz);
        
        i = std::max(0, std::min(nx - 1, i));
        j = std::max(0, std::min(ny - 1, j));
        k = std::max(0, std::min(nz - 1, k));
    }
    
    bool estaDentroCilindro(double x, double y, double z) const {
        double r = std::sqrt(x*x + y*y);
        return (r <= cilindro.getRadio() && z >= 0 && z <= cilindro.getAltura());
    }
    
    bool esPosicionValida(double x, double y, double z) const {
        return estaDentroCilindro(x, y, z);
    }

public:
    double getMasaArroz(double x, double y, double z) const {
        int i, j, k;
        aIndices(x, y, z, i, j, k);
        return masa_arroz[i][j][k];
    }
    
    void setMasaArroz(double x, double y, double z, double masa) {
        int i, j, k;
        aIndices(x, y, z, i, j, k);
        if (estaDentroCilindro(x, y, z)) {
            masa_arroz[i][j][k] = masa;
        }
    }
    
    bool estaOcupadoPorBurbuja(double x, double y, double z) const {
        int i, j, k;
        aIndices(x, y, z, i, j, k);
        return ocupado_por_burbuja[i][j][k] == 1;
    }

    // ===== GENERACIÓN DE BURBUJAS =====
    void generarBurbujasSuperficie3D(int max_burbujas_por_paso = 3) {
        if (max_burbujas_por_paso <= 0) return;
        
        int burbujas_generadas = 0;
        
        // Generar en puntos aleatorios de la superficie
        for (int intento = 0; intento < max_burbujas_por_paso * 10 && burbujas_generadas < max_burbujas_por_paso; ++intento) {
            // Elegir aleatoriamente entre fondo y paredes
            if (distribucion(generador) < 0.7) {
                // Fondo
                double theta = 2.0 * M_PI * distribucion(generador);
                double r = cilindro.getRadio() * 0.8 * distribucion(generador); // No todo el radio
                double x = r * std::cos(theta);
                double y = r * std::sin(theta);
                double z = 0.0;
                
                if (!estaOcupadoPorBurbuja(x, y, z)) {
                    double temp = cilindro.getTemperaturaCartesianas(x, y, z);
                    if (temp > temperatura_umbral_burbujas && distribucion(generador) < 0.3) {
                        generarBurbujaEnPosicion3D(x, y, z, temp);
                        burbujas_generadas++;
                    }
                }
            } else {
                // Paredes
                double theta = 2.0 * M_PI * distribucion(generador);
                double z = cilindro.getAltura() * distribucion(generador) * 0.5; // Mitad inferior
                double x = cilindro.getRadio() * std::cos(theta);
                double y = cilindro.getRadio() * std::sin(theta);
                
                if (!estaOcupadoPorBurbuja(x, y, z)) {
                    double temp = cilindro.getTemperaturaCartesianas(x, y, z);
                    if (temp > temperatura_umbral_burbujas && distribucion(generador) < 0.2) {
                        generarBurbujaEnPosicion3D(x, y, z, temp);
                        burbujas_generadas++;
                    }
                }
            }
        }
        
        if (burbujas_generadas > 0) {
            std::cout << "  Generadas " << burbujas_generadas << " burbujas" << std::endl;
        }
    }

private:
    void generarBurbujaEnPosicion3D(double x, double y, double z, double temperatura) {
        Burbuja3D nueva_burbuja(x, y, z, 0.5, siguiente_id_burbuja++);
        burbujas.push_back(nueva_burbuja);
        
        int i, j, k;
        aIndices(x, y, z, i, j, k);
        ocupado_por_burbuja[i][j][k] = 1;
    }

public:
    // ===== MOVIMIENTO 3D =====
    void moverBurbujas3D(double dt = 1.0) {
        verificarCoalescencia3D();
        
        for (auto& burbuja : burbujas) {
            if (!burbuja.activa) continue;
            moverBurbuja3D(burbuja, dt);
        }
        
        // Limpiar inactivas
        burbujas.erase(
            std::remove_if(burbujas.begin(), burbujas.end(),
                         [](const Burbuja3D& b) { return !b.activa; }),
            burbujas.end()
        );
    }

private:
    void moverBurbuja3D(Burbuja3D& burbuja, double dt) {
        double x_old = burbuja.x, y_old = burbuja.y, z_old = burbuja.z;
        
        // Movimiento simplificado
        double dx = (distribucion(generador) - 0.5) * 0.5;
        double dy = (distribucion(generador) - 0.5) * 0.5;
        double dz = 1.0;
        
        double x_new = burbuja.x + dx * dt;
        double y_new = burbuja.y + dy * dt;
        double z_new = burbuja.z + dz * dt;
        
        // Verificar límites más estrictamente
        if (!esPosicionValida(x_new, y_new, z_new)) {
            return;
        }
        
        if (z_new >= cilindro.getAltura()) {
            // CORRECCIÓN: Cuando llega a la superficie, mantener la masa
            // No llamar a redistribuirMasa3D que podría perder masa
            burbuja.activa = false;
            
            // Liberar la posición pero mantener la masa
            int i_old, j_old, k_old;
            aIndices(x_old, y_old, z_old, i_old, j_old, k_old);
            ocupado_por_burbuja[i_old][j_old][k_old] = 0;
            return;
        }
        
        if (estaOcupadoPorBurbuja(x_new, y_new, z_new)) {
            return;
        }
        
        // Actualizar posición y redistribuir masa
        burbuja.x = x_new;
        burbuja.y = y_new;
        burbuja.z = z_new;
        redistribuirMasa3D(x_old, y_old, z_old, x_new, y_new, z_new);
    }

    void verificarCoalescencia3D() {
        for (size_t i = 0; i < burbujas.size(); ++i) {
            if (!burbujas[i].activa) continue;
            
            for (size_t j = i + 1; j < burbujas.size(); ++j) {
                if (!burbujas[j].activa) continue;
                
                if (burbujas[i].distancia(burbujas[j]) < 1.0) {
                    coalescerBurbujas3D(burbujas[i], burbujas[j]);
                }
            }
        }
    }
    
    void coalescerBurbujas3D(Burbuja3D& b1, Burbuja3D& b2) {
        // Fusionar conservando volumen
        double vol1 = std::pow(b1.radio, 3);
        double vol2 = std::pow(b2.radio, 3);
        double vol_total = vol1 + vol2;
        b1.radio = std::cbrt(vol_total);
        
        // Promediar posición
        b1.x = (b1.x + b2.x) / 2.0;
        b1.y = (b1.y + b2.y) / 2.0;
        b1.z = (b1.z + b2.z) / 2.0;
        
        // CORRECCIÓN: Manejar la masa de b2 al coalescer
        double masa_b2 = getMasaArroz(b2.x, b2.y, b2.z);
        double masa_b1 = getMasaArroz(b1.x, b1.y, b1.z);
        
        // Transferir masa de b2 a b1
        setMasaArroz(b1.x, b1.y, b1.z, masa_b1 + masa_b2);
        
        // Liberar posición de b2
        int i, j, k;
        aIndices(b2.x, b2.y, b2.z, i, j, k);
        masa_arroz[i][j][k] = 0.0;  // Asegurar que se pone a cero
        ocupado_por_burbuja[i][j][k] = 0;
        
        b2.activa = false;
    }

    void redistribuirMasa3D(double x_old, double y_old, double z_old,
                           double x_new, double y_new, double z_new) {
        int i_old, j_old, k_old;
        int i_new, j_new, k_new;
        
        aIndices(x_old, y_old, z_old, i_old, j_old, k_old);
        aIndices(x_new, y_new, z_new, i_new, j_new, k_new);
        
        // Solo redistribuir si es una posición válida dentro del cilindro
        if (!estaDentroCilindro(x_new, y_new, z_new)) {
            // Si la nueva posición está fuera, mantener la masa en la posición vieja
            return;
        }
        
        double masa_vieja = masa_arroz[i_old][j_old][k_old];
        
        // SOLUCIÓN CORREGIDA: SUMAR la masa en lugar de reemplazarla
        if (i_old != i_new || j_old != j_new || k_old != k_new) {
            // Sumar la masa a la nueva posición (no reemplazar)
            masa_arroz[i_new][j_new][k_new] += masa_vieja;
            // Solo poner a cero la posición vieja si no es la misma
            masa_arroz[i_old][j_old][k_old] = 0.0;
        }
        
        // Actualizar ocupación
        ocupado_por_burbuja[i_old][j_old][k_old] = 0;
        ocupado_por_burbuja[i_new][j_new][k_new] = 1;
    }


public:
    // ===== VISUALIZACIÓN =====
    void exportarEstadoCompleto3D(const std::string& nombre_archivo) const {
        std::ofstream archivo(nombre_archivo);
        archivo << "x y z temperatura masa_arroz ocupado_burbuja" << std::endl;
        
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    double x = x_min + i * dx;
                    double y = y_min + j * dy;
                    double z_val = z_min + k * dz;
                    
                    if (estaDentroCilindro(x, y, z_val)) {
                        double temp = cilindro.getTemperaturaCartesianas(x, y, z_val);
                        archivo << x << " " << y << " " << z_val << " "
                               << temp << " "
                               << masa_arroz[i][j][k] << " "
                               << ocupado_por_burbuja[i][j][k] << std::endl;
                    }
                }
            }
        }
        archivo.close();
    }
    
    void imprimirEstadisticas3D() const {
        std::cout << "  Burbujas activas: " << burbujas.size() << std::endl;
        
        // Calcular masa total
        double masa_total = 0.0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    double x = x_min + i * dx;
                    double y = y_min + j * dy;
                    double z = z_min + k * dz;
                    if (estaDentroCilindro(x, y, z)) {
                        masa_total += masa_arroz[i][j][k];
                    }
                }
            }
        }
        std::cout << "  Masa total arroz: " << masa_total << std::endl;
    }

    // Getters
    const std::vector<Burbuja3D>& getBurbujas() const { return burbujas; }

    double getTemperaturaUmbral() const { 
        return temperatura_umbral_burbujas; 
    }
    
    double getFlotabilidadBase() const { 
        return flotabilidad_base; 
    }
    
    double getFactorConveccion() const { 
        return factor_conveccion; 
    }
    
    // Setters
    void setFlotabilidadBase(double flotabilidad) { flotabilidad_base = flotabilidad; }
    void setFactorConveccion(double factor) { factor_conveccion = factor; }
    void setTemperaturaUmbral(double temp) { temperatura_umbral_burbujas = temp; }
};

#endif