#include "Cilindro3D.h"
#include "malla_arroz_burbujas.h"
#include <iostream>
#include <chrono>

class SimuladorDebug {
private:
    Cilindro3DEvolutivo cilindro;
    MallaArrozBurbujas3DSimple malla;
    double tiempo_actual;

public:
    SimuladorDebug() 
        : cilindro(10, 8, 25.0),
          malla(cilindro, 20, 20, 16, 1.0),
          tiempo_actual(0.0) {
        
        // Configurar parámetros más sensibles
        malla.setTemperaturaUmbral(40.0);  // Más bajo para debug
    }

    void ejecutarDebug() {
        std::cout << "=== MODO DEBUG - GENERACIÓN DE BURBUJAS ===" << std::endl;
        
        // Probar diferentes temperaturas
        for (int paso = 0; paso < 10; ++paso) {
            std::cout << "\n--- Paso " << paso << " ---" << std::endl;
            
            // Simular aumento de temperatura
            double temperatura_simulada = 25.0 + paso * 10.0;
            std::cout << "Temperatura simulada: " << temperatura_simulada << "°C" << std::endl;
            
            // Forzar temperatura en el cilindro para testing
            for (int r = 0; r <= cilindro.getRadio(); ++r) {
                cilindro.setTemperatura(r, 0, temperatura_simulada); // Fondo caliente
            }
            
            // Intentar generar burbujas
            std::cout << "Intentando generar burbujas..." << std::endl;
            malla.generarBurbujasSuperficie3D(5);
            
            // Verificar estado
            const auto& burbujas = malla.getBurbujas();
            std::cout << "Burbujas activas: " << burbujas.size() << std::endl;
            
            // Mostrar detalles de las burbujas
            for (const auto& burbuja : burbujas) {
                std::cout << "  Burbuja " << burbuja.id << " en (" 
                          << burbuja.x << ", " << burbuja.y << ", " << burbuja.z << ")" << std::endl;
            }
            
            // Verificar condiciones de generación
            verificarCondicionesGeneracion();
        }
    }

private:
    void verificarCondicionesGeneracion() {
        std::cout << "--- VERIFICANDO CONDICIONES ---" << std::endl;
        
        // Verificar algunas posiciones del fondo
        int puntos_verificados = 0;
        int puntos_validos = 0;
        
        for (int r = 0; r <= 5; r += 1) { // Solo verificar algunas posiciones
            double x = r * 0.5;
            double y = 0.0;
            double z = 0.0;
            
            double temp = cilindro.getTemperaturaCartesianas(x, y, z);
            bool ocupado = malla.estaOcupadoPorBurbuja(x, y, z);
            
            std::cout << "Posición (" << x << ", " << y << ", " << z << "): " 
                      << "temp=" << temp << "°C, "
                      << "ocupado=" << (ocupado ? "SI" : "NO") << std::endl;
            
            puntos_verificados++;
            if (temp > 40.0 && !ocupado) {
                puntos_validos++;
            }
        }
        
        std::cout << "Puntos válidos para generación: " << puntos_validos << "/" << puntos_verificados << std::endl;
    }
};

int main() {
    std::cout << "DEBUG - GENERACIÓN DE BURBUJAS" << std::endl;
    
    SimuladorDebug debug;
    debug.ejecutarDebug();
    
    return 0;
}