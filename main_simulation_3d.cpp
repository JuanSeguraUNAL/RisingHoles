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
    const double ALPHA = 0.8;
    const int MAX_BURBUJAS_PASO = 15;
    const int PASOS_TOTALES = 500;  // Reducido para prueba
    const int INTERVALO_GUARDADO = 5;  // Guardar cada 5 pasos
    
    // Crear cilindro 3D
    Cilindro3D olla(RADIO_OLLA, ALTURA_OLLA, DT, ALPHA);
    olla.configurarFuenteGaussiana(15.0, 2.0);
    
    // Crear malla de burbujas 3D con carpeta de resultados
    MallaBurbujas3D malla_burbujas(olla, MAX_BURBUJAS_PASO, "resultados_simulacion");
    
    std::cout << "Iniciando simulación 3D de patrones en cocción de arroz..." << std::endl;
    std::cout << "Los resultados se guardarán en: resultados_simulacion/" << std::endl;
    std::cout << "Malla: " << ALTURA_OLLA << "x" << RADIO_OLLA << "x36 = " 
              << ALTURA_OLLA * RADIO_OLLA * 36 << " puntos por archivo" << std::endl;
    
    // Guardar geometría inicial (solo una vez)
    malla_burbujas.guardarGeometriaOlla(olla);
    
    for (int paso = 0; paso < PASOS_TOTALES; ++paso) {
        std::cout << "Paso " << paso << ":" << std::endl;
        
        // Evolución del sistema
        olla.evolucionarTemperatura();
        malla_burbujas.generarBurbujas();
        malla_burbujas.moverBurbujas();
        malla_burbujas.verificarCoalescencia();
        malla_burbujas.limpiarBurbujas();
        
        // Guardar resultados en intervalos
        if (paso % INTERVALO_GUARDADO == 0) {
            malla_burbujas.guardarInfluenciaTrayectorias(paso);
            malla_burbujas.guardarEstadoBurbujas(paso);
            malla_burbujas.guardarTemperatura(paso, olla);
            std::cout << "  -> Datos COMPLETOS guardados para paso " << paso << std::endl;
        }
        
        // Mostrar estado
        malla_burbujas.imprimirEstado();
        
        if (paso % 20 == 0) {
            std::cout << "--- Temperatura en centro de la base: " 
                      << olla.getTemperatura(0, 0) << "°C ---" << std::endl;
        }
        
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    
    // Guardar resultados finales
    malla_burbujas.guardarInfluenciaTrayectorias(PASOS_TOTALES);
    malla_burbujas.guardarEstadoBurbujas(PASOS_TOTALES);
    malla_burbujas.guardarTemperatura(PASOS_TOTALES, olla);
    
    std::cout << "Simulación 3D completada." << std::endl;
    std::cout << "Resultados guardados en: resultados_simulacion/" << std::endl;
    std::cout << "Archivos generados:" << std::endl;
    std::cout << "  - geometria_olla.txt: Estructura de la olla" << std::endl;
    std::cout << "  - influencia_paso_X.txt: Mapas de influencia completos" << std::endl;
    std::cout << "  - temperatura_paso_X.txt: Distribuciones de temperatura" << std::endl;
    std::cout << "  - burbujas_paso_X.txt: Estados de las burbujas" << std::endl;
    
    return 0;
}