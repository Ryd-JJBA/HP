#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
using namespace std;

float eps0 = 1;
float mu0 = 1;

int N = 1;
float dx = 0.2;  // Aumentado para reducir el número de pasos
float dt = 0.01; // Aumentado para reducir el número de iteraciones
int pasos = N / dx;
float tiempo = 0.1; // Reducido para simulación más rápida
float A = 1;
float omega = 1;

double eps(float x, float y, float z) {
    return eps0 * log(1 + (x * y * z));
}

double sigma(float x, float y, float z) {
    return ((dt / dx) * (dt / dx)) / (mu0 * eps(x, y, z));
}

double fuente(float omega, float t) {
    return A * cos(omega * t);
}

int main() {
    cout << "Iniciando simulación..." << endl;
    cout << "Pasos del grid: " << pasos << "x" << pasos << "x" << pasos << endl;
    cout << "Número de iteraciones temporales: " << int(tiempo / dt) << endl;
    
	//Construcción, pasado
    vector<vector<vector<double> > > Ex_old(pasos, vector<vector<double> >(pasos, vector<double>(pasos, 0.0)));
    vector<vector<vector<double> > > Ey_old(pasos, vector<vector<double> >(pasos, vector<double>(pasos, 0.0)));
    vector<vector<vector<double> > > Ez_old(pasos, vector<vector<double> >(pasos, vector<double>(pasos, 0.0)));

    vector<vector<vector<double> > > Ex_pres = Ex_old;
    vector<vector<vector<double> > > Ey_pres = Ey_old;
    vector<vector<vector<double> > > Ez_pres = Ez_old;

    vector<vector<vector<double> > > Ex_fut = Ex_old;
    vector<vector<vector<double> > > Ey_fut = Ey_old;
    vector<vector<vector<double> > > Ez_fut = Ez_old;
    
    vector<vector<vector<double> > > sigma_grid(pasos, vector<vector<double> >(pasos, vector<double>(pasos, 0.0)));
    
    // Fuente extendida en x, centrada en y y z
    int j_c = pasos / 2;
    int k_c = pasos / 2;
    for (int i = 1; i < pasos / 5; i++) {
        Ez_old[i][j_c][k_c] = fuente(omega, 0);
    }

    cout << "Construyendo primer presente..." << endl;
    // Construcción del primer presente
    for (int i = 1; i < pasos - 1; i++) {
        for (int j = 1; j < pasos - 1; j++) {
            for (int k = 1; k < pasos - 1; k++) {
                Ex_pres[i][j][k] = Ex_old[i][j][k] + 0.5 * sigma(i, j, k) * (Ex_old[i][j + 1][k] + Ex_old[i][j - 1][k] - 2 * Ex_old[i][j][k] + Ex_old[i][j][k + 1] + Ex_old[i][j][k - 1] - 2 * Ex_old[i][j][k] - 0.25 * (Ey_old[i + 1][j + 1][k] + Ey_old[i - 1][j - 1][k] - Ey_old[i + 1][j - 1][k] + Ey_old[i - 1][j + 1][k] + Ez_old[i + 1][j][k + 1] + Ez_old[i - 1][j][k - 1] - Ez_old[i + 1][j][k - 1] - Ez_old[i - 1][j][k + 1]));
                Ey_pres[i][j][k] = Ey_old[i][j][k] + 0.5 * sigma(i, j, k) * (Ey_old[i + 1][j][k] + Ey_old[i - 1][j][k] - 2 * Ey_old[i][j][k] + Ey_old[i][j][k + 1] + Ey_old[i][j][k - 1] - 2 * Ey_old[i][j][k] - 0.25 * (Ex_old[i + 1][j + 1][k] + Ex_old[i - 1][j - 1][k] - Ex_old[i + 1][j - 1][k] + Ex_old[i - 1][j + 1][k] + Ez_old[i][j + 1][k + 1] + Ez_old[i][j - 1][k - 1] - Ez_old[i][j + 1][k - 1] - Ez_old[i][j - 1][k + 1]));
                Ez_pres[i][j][k] = Ez_old[i][j][k] + 0.5 * sigma(i, j, k) * (Ez_old[i + 1][j][k] + Ez_old[i - 1][j][k] - 2 * Ez_old[i][j][k] + Ez_old[i][j + 1][k] + Ez_old[i][j - 1][k] - 2 * Ez_old[i][j][k] - 0.25 * (Ex_old[i + 1][j][k + 1] + Ex_old[i - 1][j][k - 1] - Ex_old[i + 1][j][k - 1] + Ex_old[i - 1][j][k + 1] + Ey_old[i][j + 1][k + 1] + Ey_old[i][j - 1][k - 1] - Ey_old[i][j + 1][k - 1] - Ey_old[i][j - 1][k + 1]));
            }
        }
    }

    cout << "Iniciando evolución temporal..." << endl;
    // Evolución temporal
    int total_steps = int(tiempo / dt);
    for (int t = 0; t < total_steps; t++) {
	    
        if (t % (total_steps / 10 + 1) == 0) {
            cout << "Progreso: " << (t * 100 / total_steps) << "%" << endl;
        }
        
        for (int i = 1; i < pasos / 5; i++) {
            Ez_pres[i][j_c][k_c] = fuente(omega, t * dt);
        }

        for (int i = 1; i < pasos - 1; i++) {
            for (int j = 1; j < pasos - 1; j++) {
                for (int k = 1; k < pasos - 1; k++) {
                    Ex_fut[i][j][k] = 2 * Ex_pres[i][j][k] - Ex_old[i][j][k] + sigma(i, j, k) * (Ex_pres[i][j + 1][k] + Ex_pres[i][j - 1][k] - 4 * Ex_pres[i][j][k] + Ex_pres[i][j][k + 1] + Ex_pres[i][j][k - 1] - 0.25 * (Ey_pres[i + 1][j + 1][k] + Ey_pres[i - 1][j - 1][k] - Ey_pres[i + 1][j - 1][k] - Ey_pres[i - 1][j + 1][k] + Ez_pres[i + 1][j][k + 1] + Ez_pres[i - 1][j][k - 1] - Ez_pres[i + 1][j][k - 1] - Ez_pres[i - 1][j][k + 1]));
                    Ey_fut[i][j][k] = 2 * Ey_pres[i][j][k] - Ey_old[i][j][k] + sigma(i, j, k) * (Ey_pres[i + 1][j][k] + Ey_pres[i - 1][j][k] - 4 * Ey_pres[i][j][k] + Ey_pres[i][j][k + 1] + Ey_pres[i][j][k - 1] - 0.25 * (Ex_pres[i + 1][j + 1][k] + Ex_pres[i - 1][j - 1][k] - Ex_pres[i + 1][j - 1][k] - Ex_pres[i - 1][j + 1][k] + Ez_pres[i][j + 1][k + 1] + Ez_pres[i][j - 1][k - 1] - Ez_pres[i][j + 1][k - 1] - Ez_pres[i][j - 1][k + 1]));
                    Ez_fut[i][j][k] = 2 * Ez_pres[i][j][k] - Ez_old[i][j][k] + sigma(i, j, k) * (Ez_pres[i + 1][j][k] + Ez_pres[i - 1][j][k] - 4 * Ez_pres[i][j][k] + Ez_pres[i][j + 1][k] + Ez_pres[i][j - 1][k] - 0.25 * (Ex_pres[i + 1][j][k + 1] + Ex_pres[i - 1][j][k - 1] - Ex_pres[i + 1][j][k - 1] - Ex_pres[i - 1][j][k + 1] + Ey_pres[i][j + 1][k + 1] + Ey_pres[i][j - 1][k - 1] - Ey_pres[i][j + 1][k - 1] - Ey_pres[i][j - 1][k + 1]));
                }
            }
        }

        // Condiciones de frontera absorbentes (ABC, es un método con el que básicamente conseguimos que la onda se vaya a "infinito" para que no hayan efectos de reflejo en los extremos)
        for (int i = 0; i < pasos; i++){
            for (int j = 0; j < pasos; j++){
                for (int k = 0; k < pasos; k++){
                    if (i == 0 || i == pasos - 1 || j == 0 || j == pasos - 1 || k == 0 || k == pasos - 1){
                        Ex_fut[i][j][k] = 0.0;
                        Ey_fut[i][j][k] = 0.0;
                        Ez_fut[i][j][k] = 0.0;
                    }
                }
            }
        }

        // Avanzar en el tiempo
        Ex_old = Ex_pres;
        Ey_old = Ey_pres;
        Ez_old = Ez_pres;

        Ex_pres = Ex_fut;
        Ey_pres = Ey_fut;
        Ez_pres = Ez_fut;
//Bueno, no me quiere compilar esta porquería pero esta parte serviría para guardar datos en un plano constante
	int k_plano = pasos/2;  // plano z = constante
	ofstream z_cte("frame_t_" + to_string(t) + ".dat");
	for (int i=0; i<pasos; i++){
		for (int j=0; j<pasos; j++){
			z_cte << i*dx << " " << j*dx << " " << Ez_pres[i][j][k_plano] << "\n";
		}
    z_cte << "\n";
}
z_cte.close();
	    
    }

    cout << "¡Simulación completada exitosamente!" << endl;
    cout << "Campo Ez en el centro: " << Ez_pres[pasos/2][pasos/2][pasos/2] << endl;
	
    return 0;
}
