#ifndef CILINDRO_3D_H
#define CILINDRO_3D_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

class Cilindro3D {
private:
    int radio, altura;
    double temperatura_ambiente;

    Eigen::MatrixXd temperatura;
    Eigen::MatrixXi es_frontera; // 0 si no es frontera y 1 si lo es

    // Parametros fisicos
    double rho;  // Densidad del agua [kg/m3]
    double g;    // Gravedad [m/s2]
    double beta; // Coeficiente de expansion termica [1/K]
    double mu;   // Viscocidad dinamica [Pa * s]

    // Campos de flujo
    Eigen::MatrixXd v_r;    // Velocidad radial
    Eigen::MatrixXd v_z;    // Velocidad vertical
    Eigen::MatrixXd psi;    // Funcion de corriente
    Eigen::MatrixXd omega;  // Vorticidad 

public:
    // Constructor para cilindro con simetria axial
    Cilindro3D(int radio, int altura, double temp_ambiente = 25.0) 
        : radio(radio), altura(altura), temperatura_ambiente(temp_ambiente),
          temperatura(Eigen::MatrixXd::Constant(altura, radio + 1, temp_ambiente)),
          es_frontera(Eigen::MatrixXi::Zero(altura, radio + 1)), 
          rho(1000.0), g(9.81), beta(2.07e-4), mu(1e-3),

          v_r(Eigen::MatrixXd::Zero(altura, radio + 1)),
          v_z(Eigen::MatrixXd::Zero(altura, radio + 1)),
          psi(Eigen::MatrixXd::Zero(altura, radio + 1)),
          omega(Eigen::MatrixXd::Zero(altura, radio + 1)) {
        
        inicializarFronterasCilindricas();
    }

private:
    void inicializarFronterasCilindricas() {
        // Fondo del cilindro (z = 0)
        es_frontera.row(0).setConstant(1);

        // Pared lateral (r = radio)
        es_frontera.col(radio).setConstant(1);

        // Eje central (r = 0) - condiciones de simetria
        es_frontera.col(0).setConstant(2); // 2 = condicion de simetria
    }

    void construirLaplacianoCilindrico(Eigen::SparseMatrix<double>& L, double dr, double dz) {
        int n = altura * (radio + 1);
        L.resize(n, n);
        std::vector<Eigen::Triplet<double>> triplets;

        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                int idx = z * (radio + 1) + r;

                // Condiciones de frontera (Valor fijo de Dirichlet)
                if (es_frontera(z, r) == 1 || es_frontera(z, r) == 2) {
                    triplets.push_back(Eigen::Triplet<double>(idx, idx, 1.0));
                    continue;
                }

                // Laplaciano en coordenadas cilíndricas (simetría axial):
                // ∇²f = (1/r)∂/∂r(r ∂f/∂r) + ∂²f/∂z²

                double coef_centro = 0.0;
                
                // Termino radial
                if (r > 0) {
                    double r_val = static_cast<double>(r);
                    // Derivada en direccion radial positiva
                    triplets.push_back(Eigen::Triplet<double>(idx, idx + 1, 1.0 / (dr*dr) + 1.0 / (2 * r_val * dr)));
                    // Derivada en direccion radial negativa
                    triplets.push_back(Eigen::Triplet<double>(idx, idx - 1, 1.0 / (dr*dr) - 1.0 / (2 * r_val * dr)));
                    coef_centro -= 2.0 / (dr*dr);
                } else {
                    // En el eje usar condicion de simetria
                    triplets.push_back(Eigen::Triplet<double>(idx, idx + 1, 2.0 / (dr*dr)));
                    coef_centro -= 2.0 / (dr*dr);
                }

                // Termino vertical
                if (z > 0) {
                    triplets.push_back(Eigen::Triplet<double>(idx, idx - (radio + 1), 1.0 / (dz*dz)));
                    coef_centro -= 1.0 / (dz*dz);
                }
                if (z < altura - 1) {
                    triplets.push_back(Eigen::Triplet<double>(idx, idx + (radio + 1), 1.0 / (dz*dz)));
                    coef_centro -= 1.0 / (dz*dz);
                }

                // Termino central
                triplets.push_back(Eigen::Triplet<double>(idx, idx, coef_centro));
            }
        }

        L.setFromTriplets(triplets.begin(), triplets.end());
    }

public:
    // Conversion de coordenadas

    // Convertir coordenadas cartesianas a cilindricas
    void cartesianasACilindricas(double x, double y, double z, double& r, double& theta, double& z_out) const {
        r = std::sqrt(x*x + y*y);
        theta = std::atan2(y, x);
        z_out = z;
    }

    // Convertir coordenadas cilindricas a cartesianas - CORRECCIÓN: typo "doble"
    void cilindricasACartesianas(double r, double theta, double z, double& x, double& y, double& z_out) const {
        x = r * std::cos(theta);
        y = r * std::sin(theta);
        z_out = z;
    }

    // Fuentes de calor para cilindro

    // Fuente de calor uniforme en el fondo
    void aplicarFuenteCalorUniformeFondo(double temp_uniforme) {
        temperatura.row(0).setConstant(temp_uniforme);
    }

    // Fuente de calor anular en el fondo
    void aplicarFuenteCalorAnular(double temp_max, double radio_interno, double radio_externo) {
        for (int r = 0; r <= radio; ++r) {
            if (r >= radio_interno && r <= radio_externo) {
                temperatura(0, r) = temp_max;
            }
        }
    }

    // Fuente de calor Gaussiana en el fondo (simetria radial)
    void aplicarFuenteCalorGaussianaRadial(double temperatura_max, double desviacion) {
        for (int r = 0; r <= radio; ++r) {
            double perfil = temperatura_max * std::exp(-(r * r) / (2 * desviacion * desviacion));
            temperatura(0, r) = perfil;
        }
    }

    // Calentamiento desde las paredes laterales
    void aplicarCalentamientoParedes(double temp_pared) {
        for (int z = 0; z < altura; ++z) {
            temperatura(z, radio) = temp_pared;
        }
    }

    // Acceso a temperaturas

    // Obtener temperatura en coordenadas cilindricas (r, z)
    double getTemperatura(double r, double z) const {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);

        if (r_idx < 0 || r_idx > radio || z_idx < 0 || z_idx >= altura) {
            return temperatura_ambiente;
        }
        return temperatura(z_idx, r_idx);  // CORRECCIÓN: faltaba punto y coma
    }

    // Obtener temperatura en coordenadas cartesianas
    double getTemperaturaCartesianas(double x, double y, double z) const {
        double r, theta, z_val;
        cartesianasACilindricas(x, y, z, r, theta, z_val);
        return getTemperatura(r, z_val);
    }

    // Establecer temperatura en coordenadas cilindricas
    void setTemperatura(double r, double z, double temp) {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);

        if (r_idx >= 0 && r_idx <= radio && z_idx >= 0 && z_idx < altura) {
            temperatura(z_idx, r_idx) = temp;
        }
    }

    // Verificacion de fronteras
    bool esPosicionFrontera(double r, double z) const {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);

        if (r_idx < 0 || r_idx > radio || z_idx < 0 || z_idx >= altura) return false;

        return es_frontera(z_idx, r_idx) == 1;
    }

    bool esEjeCentral(double r, double z) const {
        int r_idx = static_cast<int>(r);
        int z_idx = static_cast<int>(z);

        if (r_idx < 0 || r_idx > radio || z_idx < 0 || z_idx >= altura) return false;

        return es_frontera(z_idx, r_idx) == 2;
    }

    // Calcular el gradiente en coordenadas cilindrica
    void calcularGradienteCilindrico(Eigen::MatrixXd& grad_r, Eigen::MatrixXd& grad_z) const {
        grad_r.resize(altura, radio + 1);
        grad_z.resize(altura, radio + 1);

        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                // Manejo espacial en eje central
                if (r == 0) {
                    grad_r(z, r) = 0.0; // Dada simetria radial
                } else {
                    double left = (r > 0) ? temperatura(z, r - 1) : temperatura(z, r);
                    double right = (r < radio) ? temperatura(z, r + 1) : temperatura(z, r);
                    grad_r(z, r) = (right - left) / 2.0;
                }

                // Gradiente vertical
                double down = (z > 0) ? temperatura(z-1, r) : temperatura(z, r);
                double up = (z < altura - 1) ? temperatura(z + 1, r) : temperatura(z, r);
                grad_z(z, r) = (up - down) / 2.0;
            }
        }
    }

    // Calcular campo de velocidad para conveccion
    void calcularCampoConveccion(double dr = 0.001, double dz = 0.001) {
        // 1. Calcular gradiente radial de temperatura (Fuente de vorticidad)
        Eigen::MatrixXd dTdr = Eigen::MatrixXd::Zero(altura, radio + 1);
        
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                if (r == 0) {
                    // En el eje usar diferencia forward por simetria
                    dTdr(z, r) = (temperatura(z, 1) - temperatura(z, 0)) / dr;
                } else if (r == radio) {
                    // En la pared: usar diferencia backward 
                    dTdr(z, r) = (temperatura(z, radio) - temperatura(z, radio - 1)) / dr;
                } else {
                    // Diferencia centrada
                    dTdr(z, r) = (temperatura(z, r + 1) - temperatura(z, r - 1)) / (2 * dr);
                }
            }
        }

        // 2. Calcular el termino fuente para ecuacion de vorticidad
        // En coordenadas cilíndricas: ∇²ω = (ρ·g·β/μ) * ∂T/∂r
        Eigen::VectorXd fuente_vorticidad(altura * (radio + 1));
        for (int z = 0; z < altura; ++z) { 
            for (int r = 0; r <= radio; ++r) {
                int idx = z * (radio + 1) + r;
                if (es_frontera(z, r) == 1 || es_frontera(z, r) == 2) {
                    fuente_vorticidad(idx) = 0.0; // Condicion de Dirichlet en fronteras
                } else {
                    fuente_vorticidad(idx) = (rho * g * beta / mu) * dTdr(z, r);
                }
            }
        }

        // 3. Resolver ecuacion de vorticidad: ∇²ω = fuente
        Eigen::SparseMatrix<double> L_vort;
        construirLaplacianoCilindrico(L_vort, dr, dz);

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver_vort;
        solver_vort.compute(L_vort);

        Eigen::VectorXd omega_flat = solver_vort.solve(fuente_vorticidad);

        // Remodelar a matriz 2D
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                int idx = z * (radio + 1) + r;
                omega(z, r) = omega_flat(idx);
            }
        }

        // 4. Resolver ecuación de función de corriente: ∇²ψ = -ω
        Eigen::VectorXd fuente_corriente(altura * (radio + 1));
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                int idx = z * (radio + 1) + r;
                if (es_frontera(z, r) == 1 || es_frontera(z, r) == 2) {
                    fuente_corriente(idx) = 0.0; // ψ = 0 en fronteras
                } else {
                    fuente_corriente(idx) = -omega(z, r);
                }
            }
        }

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver_psi;
        solver_psi.compute(L_vort);  // Mismo operador Laplaciano

        Eigen::VectorXd psi_flat = solver_psi.solve(fuente_corriente);
    
        // Remodelar a matriz 2D
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                int idx = z * (radio + 1) + r;
                psi(z, r) = psi_flat(idx);
            }
        }

        // 5. Calcular velocidades a partir de la funcion de corriente
        // En coordenadas cilíndricas: v_r = -∂ψ/∂z, v_z = (1/r)∂/∂r(rψ)
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                // Velocidad radial: v_r = -∂ψ/∂z
                if (z == 0) {
                    v_r(z, r) = -(psi(z + 1, r) - psi(z, r)) / dz;
                } else if (z == altura - 1) {
                    v_r(z, r) = -(psi(z, r) - psi(z - 1, r)) / dz;
                } else {
                    v_r(z, r) = -(psi(z + 1, r) - psi(z - 1, r)) / (2 * dz);
                }

                // Velocidad vertical: v_z = (1/r)∂/∂r(rψ)
                if (r == 0) {
                    // En el eje, usar expansion de Taylor para evitar division por cero
                    v_z(z, r) = (2.0 * psi(z, 1) - 2.0 * psi(z, 0)) / (dr * dr);
                } else {
                    double r_val = static_cast<double>(r) * dr;
                    if (r == radio) {
                        v_z(z, r) = ( (r_val * psi(z, r) - (r_val - dr) * psi(z, r - 1)) ) / (r_val * dr);
                    } else {
                        double term_forward = (r_val + dr) * psi(z, r + 1);
                        double term_backward = (r_val - dr) * psi(z, r - 1);
                        v_z(z, r) = (term_forward - term_backward) / (2 * r_val * dr);
                    }
                }

                // Asegurar condiciones de frontera: Velocidad cero en paredes solidas
                if (es_frontera(z, r) == 1) {
                    v_r(z, r) = 0.0;
                    v_z(z, r) = 0.0;
                }
            }
        }
    }

    // Getters
    double getVelocidadRadial(int r, int z) const {
        if (r < 0 || r > radio || z < 0 || z >= altura) return 0.0;
        return v_r(z, r);
    }

    double getVelocidadVertical(int r, int z) const {
        if (r < 0 || r > radio || z < 0 || z >= altura) return 0.0;
        return v_z(z, r);
    }
    
    const Eigen::MatrixXd& getTemperaturaRef() const { return temperatura; }
    const Eigen::MatrixXd& getVelocidadRadialRef() const { return v_r; }
    const Eigen::MatrixXd& getVelocidadVerticalRef() const { return v_z; }
    const Eigen::MatrixXd& getFuncionCorrienteRef() const { return psi; }
    const Eigen::MatrixXd& getVorticidadRef() const { return omega; }

    // Visualizacion
    void exportarCampoFlujo(const std::string& nombre_archivo) const {
        std::ofstream archivo(nombre_archivo);
        archivo << "r z temperatura v_r v_z psi omega" << std::endl;
        
        for (int z = 0; z < altura; ++z) {
            for (int r = 0; r <= radio; ++r) {
                archivo << r << " " << z << " " << temperatura(z, r) << " "
                    << v_r(z, r) << " " << v_z(z, r) << " " 
                    << psi(z, r) << " " << omega(z, r) << std::endl;
            }
        }
        archivo.close();
    }

    void imprimirEstadisticasFlujo() const {
        std::cout << "=== ESTADÍSTICAS DE FLUJO ===" << std::endl;
        std::cout << "Velocidad radial - Max: " << v_r.maxCoeff() 
                << ", Min: " << v_r.minCoeff() << std::endl;
        std::cout << "Velocidad vertical - Max: " << v_z.maxCoeff() 
                << ", Min: " << v_z.minCoeff() << std::endl;
        std::cout << "Función corriente - Max: " << psi.maxCoeff() 
                << ", Min: " << psi.minCoeff() << std::endl;
        std::cout << "Vorticidad - Max: " << omega.maxCoeff() 
                << ", Min: " << omega.minCoeff() << std::endl;
    }
};

#endif