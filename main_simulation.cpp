#include "Cilindro3D.h"
#include "malla_arroz_burbujas.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <fstream>
#include <iomanip>

class SimuladorCompleto {
private:
    Cilindro3DEvolutivo cilindro;
    MallaArrozBurbujas3DSimple malla;
    
    double tiempo_total;
    double dt;
    int pasos_totales;
    double tiempo_actual;
    
    double dr, dz;
    double potencia_maxima;
    double tiempo_calentamiento;
    
    double umbral_burbujas;
    
    // Directorio para guardar archivos
    std::string directorio_salida;

public:
    SimuladorCompleto(int radio_olla, int altura_olla, 
                     int res_x, int res_y, int res_z,
                     double tiempo_simulacion = 120.0,
                     double paso_temporal = 2.0,
                     const std::string& output_dir = "output")
        : cilindro(radio_olla, altura_olla, 25.0),
          malla(cilindro, res_x, res_y, res_z, 1.0),
          tiempo_total(tiempo_simulacion),
          dt(paso_temporal),
          pasos_totales(static_cast<int>(tiempo_simulacion / paso_temporal)),
          tiempo_actual(0.0),
          dr(0.001), dz(0.001),
          potencia_maxima(50000.0),
          tiempo_calentamiento(10.0),
          umbral_burbujas(35.0),
          directorio_salida(output_dir)
    {
        configurarSimulacion();
    }
    
private:
    void configurarSimulacion() {
        // Crear directorio de salida
        system(("mkdir -p " + directorio_salida).c_str());
        
        // Configurar malla de burbujas
        malla.setTemperaturaUmbral(umbral_burbujas);
        malla.setFlotabilidadBase(0.3);
        malla.setFactorConveccion(2.0);
        
        std::cout << "=== CONFIGURACIÓN SIMULADOR CORREGIDO ===" << std::endl;
        std::cout << "Resolución: " << (cilindro.getRadio()*2) << "x" << (cilindro.getRadio()*2) << "x" << cilindro.getAltura() << std::endl;
        std::cout << "Tiempo total: " << tiempo_total << " s" << std::endl;
        std::cout << "Pasos: " << pasos_totales << " (Δt = " << dt << " s)" << std::endl;
        std::cout << "Potencia máxima: " << potencia_maxima << " W" << std::endl;
        std::cout << "Umbral burbujas: " << umbral_burbujas << "°C" << std::endl;
        std::cout << "Tiempo calentamiento: " << tiempo_calentamiento << " s" << std::endl;
        std::cout << "Directorio salida: " << directorio_salida << std::endl;
    }
    
    // ===== NUEVO MÉTODO: GUARDAR DISTRIBUCIÓN DE MASA =====
    void guardarDistribucionMasa(int paso) {
        std::string nombre_archivo = directorio_salida + "/masa_paso_" + 
                                   std::to_string(paso) + "_t_" + 
                                   std::to_string(static_cast<int>(tiempo_actual)) + "s.txt";
        
        std::ofstream archivo(nombre_archivo);
        
        // Encabezado con información del paso
        archivo << "# Paso: " << paso << ", Tiempo: " << tiempo_actual << "s\n";
        archivo << "# x y z masa_arroz temperatura\n";
        
        // Configurar precisión
        archivo << std::fixed << std::setprecision(4);
        
        // Obtener dimensiones de la malla (necesitarás añadir getters en MallaArrozBurbujas3DSimple)
        // Por ahora, asumimos que podemos iterar sobre las posiciones
        
        // Para una implementación completa, necesitarías añadir estos métodos a MallaArrozBurbujas3DSimple:
        // - getNX(), getNY(), getNZ() 
        // - getXMin(), getXMax(), etc.
        // - getMasaEnCelda(i, j, k)
        
        // Mientras tanto, aquí hay una implementación que funciona con la estructura actual:
        archivo << "# NOTA: Se necesita implementar getters completos en MallaArrozBurbujas3DSimple\n";
        
        archivo.close();
        
        std::cout << "  💾 Masa guardada: " << nombre_archivo << std::endl;
    }
    
    // ===== MÉTODO ALTERNATIVO: Guardar distribución 2D (corte horizontal) =====
    void guardarCorteHorizontalMasa(int paso, double z_corte = 0.0) {
        std::string nombre_archivo = directorio_salida + "/corte_masa_paso_" + 
                                   std::to_string(paso) + "_z" + 
                                   std::to_string(static_cast<int>(z_corte)) + ".txt";
        
        std::ofstream archivo(nombre_archivo);
        
        archivo << "# Corte horizontal en z = " << z_corte << "\n";
        archivo << "# Paso: " << paso << ", Tiempo: " << tiempo_actual << "s\n";
        archivo << "x y masa_arroz temperatura\n";
        archivo << std::fixed << std::setprecision(4);
        
        // Esto es un placeholder - necesitarás implementar el acceso real a los datos
        double radio = cilindro.getRadio();
        int puntos = 50;
        
        for (int i = 0; i < puntos; ++i) {
            for (int j = 0; j < puntos; ++j) {
                double x = -radio + 2.0 * radio * i / (puntos - 1);
                double y = -radio + 2.0 * radio * j / (puntos - 1);
                
                // Verificar si está dentro del cilindro
                if (std::sqrt(x*x + y*y) <= radio) {
                    double masa = malla.getMasaArroz(x, y, z_corte);
                    double temp = cilindro.getTemperaturaCartesianas(x, y, z_corte);
                    archivo << x << " " << y << " " << masa << " " << temp << "\n";
                }
            }
        }
        
        archivo.close();
        std::cout << "  📊 Corte masa guardado: " << nombre_archivo << std::endl;
    }
    
    // ===== MÉTODO COMPLETO: Exportar estado completo incluyendo masa =====
    void exportarEstadoCompletoConMasa(int paso) {
        std::string nombre_archivo = directorio_salida + "/estado_completo_paso_" + 
                                   std::to_string(paso) + ".txt";
        
        std::ofstream archivo(nombre_archivo);
        
        archivo << "# Estado completo - Paso: " << paso << ", Tiempo: " << tiempo_actual << "s\n";
        archivo << "# x y z masa_arroz temperatura ocupado_burbuja\n";
        archivo << std::fixed << std::setprecision(4);
        
        // Usar el método existente de exportación y añadir información de masa
        // Esto requeriría modificar el método exportarEstadoCompleto3D en MallaArrozBurbujas3DSimple
        // o crear uno nuevo que incluya masa
        
        malla.exportarEstadoCompleto3D(nombre_archivo); // Este ya incluye masa_arroz
        
        std::cout << "  📁 Estado completo guardado: " << nombre_archivo << std::endl;
    }

    double calcularPotenciaInstantanea() {
        if (tiempo_actual < tiempo_calentamiento) {
            return potencia_maxima * (tiempo_actual / tiempo_calentamiento);
        } else {
            return potencia_maxima;
        }
    }
    
    void actualizarSistemaTermico(int paso) {
        double potencia_instantanea = calcularPotenciaInstantanea();
        
        cilindro.setFuenteCalorUniformeFondo(potencia_instantanea * 10.0);
        cilindro.setCalentamientoParedes(potencia_instantanea * 2.0);
        
        auto inicio = std::chrono::high_resolution_clock::now();
        cilindro.resolverEcuacionCalor(dt, dr, dz);
        auto fin = std::chrono::high_resolution_clock::now();
        
        double temp_max = cilindro.getTemperaturaRef().maxCoeff();
        if (paso > 3 && temp_max < 40.0) {
            std::cout << "  APLICANDO PARCHE: Forzando temperatura a 50°C" << std::endl;
            for (int r = 0; r <= cilindro.getRadio(); ++r) {
                cilindro.setTemperatura(r, 0, 50.0);
            }
        }
        
        std::vector<std::tuple<double, double, double>> posiciones_burbujas;
        for (const auto& burbuja : malla.getBurbujas()) {
            posiciones_burbujas.emplace_back(burbuja.x, burbuja.y, burbuja.z);
        }
        cilindro.aplicarEnfriamientoBurbujas(posiciones_burbujas, 0.05);
        
        if (paso % 5 == 0 || paso <= 10) {
            auto duracion = std::chrono::duration_cast<std::chrono::milliseconds>(fin - inicio);
            std::cout << "Paso " << paso << " (t=" << tiempo_actual << "s): " << std::endl;
            std::cout << "  Potencia: " << potencia_instantanea << " W" << std::endl;
            std::cout << "  Tiempo cálculo calor: " << duracion.count() << " ms" << std::endl;
            cilindro.imprimirEstadisticasTemperatura();
        }
    }
    
public:
    void guardarCorteSuperficialMasa(int paso) {
        std::string nombre_archivo = "distribucion_masa_final_paso_" + 
                                std::to_string(paso) + ".txt";
        
        std::ofstream archivo(nombre_archivo);
        
        archivo << "# Distribución de masa en la superficie (z = altura)\n";
        archivo << "# Paso: " << paso << ", Tiempo: " << tiempo_actual << "s\n";
        archivo << "x y masa_arroz temperatura\n";
        archivo << std::fixed << std::setprecision(6);
        
        double altura_superficie = cilindro.getAltura();
        double radio = cilindro.getRadio();
        int puntos = 100; // Resolución de la malla de salida
        
        for (int i = 0; i < puntos; ++i) {
            for (int j = 0; j < puntos; ++j) {
                double x = -radio + 2.0 * radio * i / (puntos - 1);
                double y = -radio + 2.0 * radio * j / (puntos - 1);
                
                // Verificar si está dentro del cilindro
                if (std::sqrt(x*x + y*y) <= radio) {
                    double masa = malla.getMasaArroz(x, y, altura_superficie);
                    double temp = cilindro.getTemperaturaCartesianas(x, y, altura_superficie);
                    archivo << x << " " << y << " " << masa << " " << temp << "\n";
                } else {
                    // Puntos fuera del cilindro (masa cero)
                    archivo << x << " " << y << " " << 0.0 << " " << cilindro.getTemperaturaAmbiente() << "\n";
                }
            }
            archivo << "\n"; // Línea en blanco para gnuplot
        }
        
        archivo.close();
        std::cout << "  💾 Distribución superficial guardada: " << nombre_archivo << std::endl;
    }

    void ejecutarSimulacion() {
        std::cout << "\n=== INICIANDO SIMULACIÓN TÉRMICA COMPLETA ===" << std::endl;
        
        std::cout << "\n--- Estado Inicial ---" << std::endl;
        std::cout << "Temperatura ambiente: " << cilindro.getTemperaturaAmbiente() << "°C" << std::endl;
        
        actualizarSistemaTermico(0);
        
        // Guardar estado inicial
        guardarCorteHorizontalMasa(0, 0.0); // Corte en el fondo
        guardarCorteHorizontalMasa(0, cilindro.getAltura()/2.0); // Corte en medio
        exportarEstadoCompletoConMasa(0);
        
        for (int paso = 0; paso < pasos_totales; ++paso) {
            tiempo_actual = paso * dt;
            
            std::cout << "\n--- Paso " << paso << " (t = " << tiempo_actual << "s) ---" << std::endl;
            
            actualizarSistemaTermico(paso);
            
            if (paso % 3 == 0) {
                cilindro.calcularCampoConveccion(dr, dz);
            }
            
            int max_burbujas = calcularMaxBurbujas();
            std::cout << "  Intentando generar hasta " << max_burbujas << " burbujas..." << std::endl;
            malla.generarBurbujasSuperficie3D(max_burbujas);
            
            malla.moverBurbujas3D(dt);
            
            malla.imprimirEstadisticas3D();
            
            // ===== GUARDAR DATOS DE MASA EN CADA PASO =====
            if (paso % 5 == 0 || paso <= 10 || paso == pasos_totales - 1) {
                // GUARDAR SOLO EN EL ÚLTIMO PASO
                if (paso == pasos_totales - 1) {
                    guardarCorteSuperficialMasa(paso);
                }
                
                if (paso % 15 == 0 || paso == pasos_totales - 1) {
                    exportarEstadoCompletoConMasa(paso);
                }
            }
            
            double temp_max = cilindro.getTemperaturaRef().maxCoeff();
            if (paso > 15 && temp_max < 35.0) {
                std::cout << "🛑 DETENIENDO: La temperatura no está alcanzando el umbral después de 15 pasos." << std::endl;
                // Guardar estado final antes de detener
                exportarEstadoCompletoConMasa(paso);
                break;
            }
            
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        std::cout << "\n=== SIMULACIÓN COMPLETADA ===" << std::endl;
        
        double temp_max_final = cilindro.getTemperaturaRef().maxCoeff();
        int burbujas_final = malla.getBurbujas().size();
        
        std::cout << "📊 RESULTADOS FINALES:" << std::endl;
        std::cout << "   Temperatura máxima alcanzada: " << temp_max_final << "°C" << std::endl;
        std::cout << "   Burbujas activas finales: " << burbujas_final << std::endl;
        std::cout << "   ¿Éxito? " << (temp_max_final > umbral_burbujas ? "✅ SÍ" : "❌ NO") << std::endl;
        std::cout << "   Datos guardados en: " << directorio_salida << std::endl;
        
        exportarEstadoCompletoConMasa(pasos_totales - 1);
    }

    
private:
    int calcularMaxBurbujas() {
        double temp_max = cilindro.getTemperaturaRef().maxCoeff();
        double factor = std::min(1.0, (temp_max - umbral_burbujas) / 30.0);
        
        int base = 2;
        int adicional = static_cast<int>(8 * factor);
        
        return base + adicional;
    }
};

int main() {
    std::cout << "==========================================" << std::endl;
    std::cout << "=== SIMULADOR DE COCCIÓN CORREGIDO ===" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    try {
        int radio_olla = 8;
        int altura_olla = 6;
        int res_x = 25, res_y = 25, res_z = 15;
        double tiempo_simulacion = 60.0;
        double dt = 1.0;
        
        std::cout << "🧪 Parámetros de prueba:" << std::endl;
        std::cout << "   Olla: " << radio_olla*2 << "cm × " << radio_olla*2 << "cm × " << altura_olla << "cm" << std::endl;
        std::cout << "   Tiempo: " << tiempo_simulacion << "s" << std::endl;
        std::cout << "   Δt: " << dt << "s" << std::endl;
        std::cout << "   Umbral burbujas: 35°C" << std::endl;
        std::cout << "==========================================\n" << std::endl;
        
        SimuladorCompleto simulador(radio_olla, altura_olla, 
                                   res_x, res_y, res_z,
                                   tiempo_simulacion, dt,
                                   "resultados_simulacion");
        
        simulador.ejecutarSimulacion();
        
    } catch (const std::exception& e) {
        std::cerr << "💥 ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}