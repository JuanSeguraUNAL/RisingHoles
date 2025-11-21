#include "Cilindro3D.h"
#include "MallaBurbujas3D.h"
#include <iostream>
#include <chrono>
#include <thread>

int main() {
    // Parámetros de simulación 3D - CON COEFICIENTES DIFERENTES
    const int RADIO_INTERIOR = 12; 
    const int GROSOR_PARED = 3;   
    const int ALTURA_OLLA = 20;
    const double DT = 0.05;
    
    // COEFICIENTES DE DIFUSIVIDAD TÉRMICA DIFERENTES
    // Valores típicos (puedes ajustarlos según el material específico):
    // - Agua: ~0.15 mm²/s
    // - Acero inoxidable: ~4-5 mm²/s  
    const double ALPHA_METAL = 5;    // Mayor difusividad para metal
    const double ALPHA_AGUA = 0.15;    // Menor difusividad para agua/arroz
    
    const int MAX_BURBUJAS_PASO = 6;
    const int PASOS_TOTALES = 500;
    const int INTERVALO_GUARDADO = 5;
    
    // Crear cilindro 3D con coeficientes diferentes
    Cilindro3D olla(RADIO_INTERIOR, GROSOR_PARED, ALTURA_OLLA, DT, ALPHA_METAL, ALPHA_AGUA);
    olla.configurarFuenteGaussiana(12.0, 4.0);
    
    // Crear malla de burbujas 3D con carpeta de resultados
    MallaBurbujas3D malla_burbujas(olla, MAX_BURBUJAS_PASO, "resultados_simulacion");
    
    std::cout << "==========================================================" << std::endl;
    std::cout << "INICIANDO SIMULACIÓN 3D CON COEFICIENTES DE DIFUSIVIDAD DIFERENTES" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Materiales y coeficientes:" << std::endl;
    std::cout << "  - Metal (olla): alpha = " << ALPHA_METAL << std::endl;
    std::cout << "  - Agua/arroz (interior): alpha = " << ALPHA_AGUA << std::endl;
    std::cout << "  - Relación metal/agua: " << ALPHA_METAL/ALPHA_AGUA << "x" << std::endl;
    std::cout << std::endl;
    std::cout << "Dimensiones de la olla:" << std::endl;
    std::cout << "  - Radio interior: " << olla.getRadioInterior() << std::endl;
    std::cout << "  - Grosor pared: " << olla.getGrosorPared() << std::endl;
    std::cout << "  - Radio total: " << olla.getRadioTotal() << std::endl;
    std::cout << "  - Altura: " << olla.getAltura() << std::endl;
    std::cout << std::endl;
    std::cout << "Parámetros de simulación:" << std::endl;
    std::cout << "  - Pasos totales: " << PASOS_TOTALES << std::endl;
    std::cout << "  - Delta tiempo: " << DT << std::endl;
    std::cout << "  - Máx burbujas/paso: " << MAX_BURBUJAS_PASO << std::endl;
    std::cout << "  - Intervalo guardado: " << INTERVALO_GUARDADO << std::endl;
    std::cout << std::endl;
    std::cout << "Los resultados se guardarán en: resultados_simulacion/" << std::endl;
    std::cout << "==========================================================" << std::endl;
    
    // Información de dimensiones de la malla
    std::cout << "Malla de simulación: " << ALTURA_OLLA << "x" << olla.getRadioTotal() << "x36 = " 
              << ALTURA_OLLA * olla.getRadioTotal() * 36 << " puntos por archivo" << std::endl;
    
    // Guardar geometría inicial (solo una vez)
    std::cout << "Guardando geometría de la olla..." << std::endl;
    malla_burbujas.guardarGeometriaOlla(olla);
    
    // Diagnóstico inicial
    std::cout << std::endl << "=== DIAGNÓSTICO INICIAL ===" << std::endl;
    olla.diagnosticoParedesCompleto();
    olla.diagnosticoCoeficientes();
    
    std::cout << std::endl << "=== INICIANDO SIMULACIÓN ===" << std::endl;
    
    // Bucle principal de simulación
    for (int paso = 0; paso < PASOS_TOTALES; ++paso) {
        std::cout << "Paso " << paso << ":" << std::endl;
        
        // Evolución del sistema
        olla.evolucionarTemperatura();
        malla_burbujas.generarBurbujas();
        malla_burbujas.moverBurbujas();
        malla_burbujas.verificarCoalescencia();
        malla_burbujas.limpiarBurbujas();

        //
        if (paso % 20 == 0) {
            std::cout << "=== ESTADO PAREDES METÁLICAS ===" << std::endl;
            
            // Temperaturas clave en las paredes
            std::cout << "Base metálica (r=" << olla.getRadioInterior() << "): " 
                    << olla.getTemperatura(0, olla.getRadioInterior()) << "°C" << std::endl;
            std::cout << "Pared media (h=" << olla.getAltura()/2 << ", r=" << olla.getRadioInterior() << "): " 
                    << olla.getTemperatura(olla.getAltura()/2, olla.getRadioInterior()) << "°C" << std::endl;
            std::cout << "Pared superior (h=" << olla.getAltura()-1 << ", r=" << olla.getRadioInterior() << "): " 
                    << olla.getTemperatura(olla.getAltura()-1, olla.getRadioInterior()) << "°C" << std::endl;
            
            // Verificar capacidad de generar burbujas
            double temp_umbral = 50.0;
            bool paredes_generan_burbujas = olla.getTemperatura(0, olla.getRadioInterior()) >= temp_umbral;
            std::cout << "¿Paredes pueden generar burbujas? " << (paredes_generan_burbujas ? "SÍ" : "NO") << std::endl;
        }
        
        // Guardar resultados en intervalos
        if (paso % INTERVALO_GUARDADO == 0) {
            std::cout << "  -> Guardando datos completos..." << std::endl;
            malla_burbujas.guardarInfluenciaTrayectorias(paso);
            malla_burbujas.guardarEstadoBurbujas(paso);
            malla_burbujas.guardarTemperatura(paso, olla);
            std::cout << "  -> Datos COMPLETOS guardados para paso " << paso << std::endl;
        }
        
        // Mostrar estado actual
        malla_burbujas.imprimirEstado();
        
        // Diagnóstico detallado cada 50 pasos
        if (paso % 50 == 0) {
            std::cout << "--- Diagnóstico Detallado Paso " << paso << " ---" << std::endl;
            
            // Temperaturas clave
            std::cout << "Temperaturas clave:" << std::endl;
            std::cout << "  - Centro base (agua): " << olla.getTemperatura(0, 0) << "°C" 
                      << " [alpha=" << olla.getAlpha(0, 0) << "]" << std::endl;
            std::cout << "  - Borde base (agua): " << olla.getTemperatura(0, olla.getRadioInterior()-1) << "°C" 
                      << " [alpha=" << olla.getAlpha(0, olla.getRadioInterior()-1) << "]" << std::endl;
            std::cout << "  - Pared exterior base (metal): " << olla.getTemperatura(0, olla.getRadioTotal()-1) << "°C" 
                      << " [alpha=" << olla.getAlpha(0, olla.getRadioTotal()-1) << "]" << std::endl;
            std::cout << "  - Centro volumen (agua): " << olla.getTemperatura(olla.getAltura()/2, 0) << "°C" 
                      << " [alpha=" << olla.getAlpha(olla.getAltura()/2, 0) << "]" << std::endl;
            
            // Gradientes de temperatura
            auto grad_centro = olla.getGradienteTemperatura(0, 0);
            auto grad_pared = olla.getGradienteTemperatura(0, olla.getRadioTotal()-1);
            std::cout << "Gradientes de temperatura:" << std::endl;
            std::cout << "  - Centro base: dh=" << grad_centro[0] << ", dr=" << grad_centro[1] << std::endl;
            std::cout << "  - Pared base: dh=" << grad_pared[0] << ", dr=" << grad_pared[1] << std::endl;
            
            // Diagnóstico completo cada 100 pasos
            if (paso % 100 == 0) {
                std::cout << "--- DIAGNÓSTICO COMPLETO PASO " << paso << " ---" << std::endl;
                olla.diagnosticoParedesCompleto();
                olla.diagnosticoCoeficientes();
                
                // Mostrar distribución de temperatura en la base
                std::cout << "Distribución de temperatura en base (r=0 a " << olla.getRadioTotal()-1 << "):" << std::endl;
                for (int r = 0; r < olla.getRadioTotal(); ++r) {
                    std::string tipo = (r < olla.getRadioInterior()) ? "A" : "M"; // A=Agua, M=Metal
                    std::cout << "r" << tipo << r << "=" << olla.getTemperatura(0, r) << " ";
                    if ((r+1) % 5 == 0) std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        }
        
        // Pequeña pausa para visualización en tiempo real (opcional)
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    
    std::cout << std::endl << "=== FINALIZANDO SIMULACIÓN ===" << std::endl;
    
    // Guardar resultados finales
    std::cout << "Guardando resultados finales..." << std::endl;
    malla_burbujas.guardarInfluenciaTrayectorias(PASOS_TOTALES);
    malla_burbujas.guardarEstadoBurbujas(PASOS_TOTALES);
    malla_burbujas.guardarTemperatura(PASOS_TOTALES, olla);
    
    // Diagnóstico final
    std::cout << std::endl << "=== DIAGNÓSTICO FINAL ===" << std::endl;
    olla.diagnosticoParedesCompleto();
    olla.diagnosticoCoeficientes();
    
    // Resumen final
    std::cout << std::endl << "==========================================================" << std::endl;
    std::cout << "SIMULACIÓN 3D CON COEFICIENTES DIFERENTES COMPLETADA" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Resultados guardados en: resultados_simulacion/" << std::endl;
    std::cout << std::endl;
    std::cout << "Archivos generados:" << std::endl;
    std::cout << "  - geometria_olla.txt: Estructura completa de la olla" << std::endl;
    std::cout << "  - influencia_paso_X.txt: Mapas de influencia de trayectorias" << std::endl;
    std::cout << "  - temperatura_paso_X.txt: Distribuciones de temperatura 3D" << std::endl;
    std::cout << "  - burbujas_paso_X.txt: Estados y posiciones de burbujas" << std::endl;
    std::cout << std::endl;
    std::cout << "Estadísticas finales:" << std::endl;
    std::cout << "  - Total de pasos simulados: " << PASOS_TOTALES << std::endl;
    std::cout << "  - Burbujas activas finales: " << malla_burbujas.getNumBurbujasActivas() << std::endl;
    std::cout << "  - Temperatura máxima en base: " << olla.getTemperatura(0, 0) << "°C" << std::endl;
    std::cout << "  - Diferencial metal/agua en base: " 
              << (olla.getTemperatura(0, olla.getRadioTotal()-1) - olla.getTemperatura(0, 0)) 
              << "°C" << std::endl;
    std::cout << "==========================================================" << std::endl;
    
    return 0;
}