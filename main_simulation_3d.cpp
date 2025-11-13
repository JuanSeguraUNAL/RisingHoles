#include "Cilindro3D.h"
#include "MallaBurbujas3D.h"
#include <iostream>
#include <chrono>
#include <thread>

int main() {
    // Parámetros de simulación 3D
    const int RADIO_OLLA = 15;
    const int ALTURA_OLLA = 25;
    const double DT = 0.05;
    const double ALPHA = 0.01;
    const int MAX_BURBUJAS_PASO = 15;
    const int PASOS_TOTALES = 500;
    
    // Crear cilindro 3D
    Cilindro3D olla(RADIO_OLLA, ALTURA_OLLA, DT, ALPHA);
    olla.configurarFuenteGaussiana(15.0, 2.0);
    
    // Crear malla de burbujas 3D
    MallaBurbujas3D malla_burbujas(olla, MAX_BURBUJAS_PASO);
    
    std::cout << "Iniciando simulación 3D de patrones en cocción de arroz..." << std::endl;
    
    for (int paso = 0; paso < PASOS_TOTALES; ++paso) {
        std::cout << "Paso " << paso << ":" << std::endl;
        
        // Evolución del sistema
        olla.evolucionarTemperatura();
        malla_burbujas.generarBurbujas();
        malla_burbujas.moverBurbujas();
        malla_burbujas.verificarCoalescencia();
        malla_burbujas.limpiarBurbujas();
        
        // Mostrar estado
        malla_burbujas.imprimirEstado();
        
        if (paso % 50 == 0) {
            std::cout << "--- Temperatura en centro de la base: " 
                      << olla.getTemperatura(0, 0) << "°C ---" << std::endl;
            
            // Mostrar algunas burbujas como ejemplo
            int contador = 0;
            for (const auto& burbuja : malla_burbujas.getBurbujas()) {
                if (burbuja.isActiva() && contador < 3) {
                    std::cout << "Burbuja " << burbuja.getId() 
                              << " en (" << burbuja.getPosX() << ", " 
                              << burbuja.getPosY() << ", " << burbuja.getPosZ() 
                              << ") tamaño: " << burbuja.getTamano() << std::endl;
                    contador++;
                }
            }
        }
        
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    
    std::cout << "Simulación 3D completada." << std::endl;
    return 0;
}