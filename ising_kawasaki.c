#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define N 64 
#define PASOS 100000

double T_inicial = 0.25;
double T_final = 4.0;
double T_incremento = 0.25;

void inicializarRed_mitad(int red[N][N]);
void Guardar_Red(FILE *file, int red[N][N]);
double Energia_local(int red[N][N], int i1, int j1, int i2, int j2);
double calcularEnergia_Total_Inicial(int red[N][N]);
double magnetizacion_mitad_superior(int red[N][N]);
double magnetizacion_mitad_inferior(int red[N][N]);
void calcularDensidadesFila(int red[N][N], int fila, int *densidadpositivo, int *densidadnegativo);
void calcular_CalorEspecifico(double sumaEnergia, double sumaEnergiaCuadrada, int conteo, double T, double *cv);
void calcular_Susceptibilidad(double sumaMag, double sumaMag2, int conteo, double T, double *chi);
int algoritmoMetropolis(int red[N][N], double *sumaEnergia, double *sumaEnergiacuadrada, double *sumaMagnetizacionSuperior, double *sumaMagnetizacionInferior, int *configuracionesCambiadas, double *sumadensidadpositivo, double *sumadensidadnegativo, double *sumaMagnetizacion2, double T);

// --- MAIN AL PRINCIPIO ---
int main() {

    int red[N][N];
    srand(time(NULL)); // Semilla para números aleatorios
    static const char* tiempo_file = "tiempos.txt";
    static int primera_vez = 1;
    static clock_t tiempo_inicio, tiempo_fin;
    static double tiempo_total;

     tiempo_inicio = clock();

    // Archivos con nombre dinámico según N (no según temperatura)
    char fname_magsup[64], fname_maginfer[64], fname_denspos[64], fname_densneg[64], fname_cv[64], fname_chi[64];
    sprintf(fname_magsup, "promedio_magnetizacionsuperior%d.txt", N);
    sprintf(fname_maginfer, "promedio_magnetizacioninferior%d.txt", N);
    sprintf(fname_denspos, "promedio_densidadpositivo%d.txt", N);
    sprintf(fname_densneg, "promedio_densidadnegativo%d.txt", N);
    sprintf(fname_cv, "filecv%d.txt", N);
    sprintf(fname_chi, "susceptibilidad%d.txt", N);

    FILE *magnetizacionSuperiorFile = fopen(fname_magsup, "w");
    FILE *magnetizacionInferiorFile = fopen(fname_maginfer, "w");
    FILE *filedensidadpositivo = fopen(fname_denspos, "w");
    FILE *filedensidadnegativo = fopen(fname_densneg, "w");
    FILE *filecv = fopen(fname_cv, "w");
    FILE *filechi = fopen(fname_chi, "w");

    for (double T = T_inicial; T <= T_final + 1e-8; T += T_incremento) {
        // Inicializar la red: mitad positivos, mitad negativos
        inicializarRed_mitad(red);

        // Variables para acumular magnetización
        double sumaMagnetizacionSuperior = 0.0;
        double sumaMagnetizacionInferior = 0.0;
        int configuracionesCambiadas = 0;
        double sumaEnergia = 0.0;
        double sumaEnergiacuadrada = 0.0;
        double sumadensidadpositivo = 0.0;
        double sumadensidadnegativo = 0.0;
        double cv, chi, sumaMagnetizacion2 = 0.0;

        // Ejecutar el algoritmo de Monte Carlo
        algoritmoMetropolis(red, &sumaEnergia, &sumaEnergiacuadrada, 
            &sumaMagnetizacionSuperior, &sumaMagnetizacionInferior, 
            &configuracionesCambiadas, &sumadensidadpositivo, 
            &sumadensidadnegativo, &sumaMagnetizacion2, T);

        // Calcular promedios
        double promedioMagnetizacionSuperior = sumaMagnetizacionSuperior / configuracionesCambiadas;
        double promedioMagnetizacionInferior = sumaMagnetizacionInferior / configuracionesCambiadas;
        double promedioDensidadPositiva = sumadensidadpositivo / (configuracionesCambiadas*N);
        double promedioDensidadNegativa = sumadensidadnegativo / (configuracionesCambiadas*N);

        // Calor específico y susceptibilidad
        calcular_CalorEspecifico(sumaEnergia, sumaEnergiacuadrada, configuracionesCambiadas, T, &cv);
        calcular_Susceptibilidad(sumaMagnetizacionSuperior, sumaMagnetizacion2, configuracionesCambiadas, T, &chi);

        // Guardar resultados en los archivos (una línea por temperatura)
        fprintf(magnetizacionSuperiorFile, "%.6f %.2f\n", promedioMagnetizacionSuperior, T);
        fprintf(magnetizacionInferiorFile, "%.6f %.2f\n", promedioMagnetizacionInferior, T);
        fprintf(filedensidadnegativo, "%.6f %.2f\n", promedioDensidadNegativa, T);
        fprintf(filedensidadpositivo, "%.6f %.2f\n", promedioDensidadPositiva, T);
        fprintf(filecv, "%.6f %.2f\n", cv, T);
        fprintf(filechi, "%.10f %.2f\n", chi, T);

        printf("T=%.2f completado\n", T);
    }

    fclose(magnetizacionSuperiorFile);
    fclose(magnetizacionInferiorFile);
    fclose(filedensidadpositivo);
    fclose(filedensidadnegativo);
    fclose(filecv);
    fclose(filechi);

    printf("Simulación completada.\n");

    tiempo_fin = clock();
    tiempo_total = (double)(tiempo_fin - tiempo_inicio) / CLOCKS_PER_SEC;

    FILE *ftiempos = fopen(tiempo_file, "a");
    if (ftiempos != NULL) {
        fprintf(ftiempos, "%d\t%.6f\n", N, tiempo_total);
        fclose(ftiempos);
    } else {
        printf("No se pudo abrir el archivo de tiempos.\n");
    }

    return 0;
    return 0;
}

// --- FUNCIONES AUXILIARES CON NOMBRES COMO EN ising_kawasaki.c ---

void inicializarRed_mitad(int red[N][N]) {
    int espines_totales = (N-2) * N;
    int mitad = espines_totales / 2;
    int positivos = 0, negativos = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0) {
                red[i][j] = -1; // Borde superior con espines negativos
            } else if (i == N - 1) {
                red[i][j] = 1; // Borde inferior con espines positivos
            } else {
                if (positivos < mitad) {
                    red[i][j] = -1;
                    positivos++;
                } else {
                    red[i][j] = 1;
                    positivos++;
                }
            }
        }
    }
}

void Guardar_Red(FILE *file, int red[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%d", red[i][j]);
            if (j < N - 1) {
                fprintf(file, ", ");
            }
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

double Energia_local(int red[N][N], int i1, int j1, int i2, int j2 ) {
    int E = 0;

    // Vecinos del primer espín (i1, j1)
    int arriba1, abajo1, izquierda1, derecha1;

    if (i1 == 1) arriba1 = -1;
    else arriba1 = red[i1-1][j1];

    if (i1 == N - 2) abajo1 = 1;
    else abajo1 = red[i1+1][j1];

    if (j1 == 0) izquierda1 = red[i1][N-1];
    else izquierda1 = red[i1][j1-1];

    if (j1 == N - 1) derecha1 = red[i1][0];
    else derecha1 = red[i1][j1+1];

    // Vecinos del segundo espín (i2, j2)
    int arriba2, abajo2, izquierda2, derecha2;

    if (i2 == 1) arriba2 = -1;
    else arriba2 = red[i2-1][j2];

    if (i2 == N - 2) abajo2 = 1;
    else abajo2 = red[i2+1][j2];

    if (j2 == 0) izquierda2 = red[i2][N-1];
    else izquierda2 = red[i2][j2-1];

    if (j2 == N - 1) derecha2 = red[i2][0];
    else derecha2 = red[i2][j2+1];

    if ((i1 == i2 && (j1 + 1) % N == j2)) {
        derecha1 = 0;
        izquierda2 = 0;
    }
    if ((i1 == i2 && (j1 - 1 + N) % N == j2)) {
        izquierda1 = 0;
        derecha2 = 0;
    }
    if ((i1 + 1 == i2 && j1 == j2)) {
        arriba2 = 0;
        abajo1 = 0;
    }
    if ((i1 - 1 == i2 && j1 == j2)) {
        abajo2 = 0;
        arriba1 = 0;
    }

    E = -red[i1][j1] * (arriba1 + abajo1 + izquierda1 + derecha1);
    E += -red[i2][j2] * (arriba2 + abajo2 + izquierda2 + derecha2);

    return E;
}

double calcularEnergia_Total_Inicial(int red[N][N]) {
    double energia = 0.0;
    for (int i = 1; i < N-1; i++) {
        for (int j = 0; j < N; j++) {
            int sumaVecinos = red[i-1][j] + red[i+1][j] +
                              red[i][(j + 1) % N] + red[i][(j - 1 + N) % N];
            energia += red[i][j] * sumaVecinos;
        }
    }
    return -0.5 * energia;
}

double magnetizacion_mitad_superior(int red[N][N]) {
    int suma = 0;
    for (int i = 0; i < N / 2; i++) {
        for (int j = 0; j < N; j++) {
            suma += red[i][j];
        }
    }
    return (double) fabs(suma) / ((N / 2) * N);
}

double magnetizacion_mitad_inferior(int red[N][N]) {
    int suma = 0;
    for (int i = N / 2; i < N; i++) {
        for (int j = 0; j < N; j++) {
            suma += red[i][j];
        }
    }
    return (double) fabs(suma) / ((N / 2) * N);
}

void calcularDensidadesFila(int red[N][N], int fila, int *densidadpositivo, int *densidadnegativo) {
    int positivos = 0;
    int negativos = 0;
    for (int j = 0; j < N; j++) {
        if (red[fila][j] == 1) positivos++;
        else if (red[fila][j] == -1) negativos++;
    }
    *densidadpositivo = positivos;
    *densidadnegativo = negativos;
}

void calcular_CalorEspecifico(double sumaEnergia, double sumaEnergiaCuadrada, int conteo, double T, double *cv) {
    double promedio_E = sumaEnergia / conteo;
    double promedio_E2 = sumaEnergiaCuadrada / conteo;
    double varianza_E = promedio_E2 - (promedio_E * promedio_E);
    *cv = varianza_E / ((N-2)*N*T*T);
}

void calcular_Susceptibilidad(double sumaMag, double sumaMag2, int conteo, double T, double *chi) {
    double promedio_M = sumaMag / conteo;
    double promedio_M2 = sumaMag2 / conteo;
    double varianza_M = promedio_M2 - (promedio_M * promedio_M);
    *chi = varianza_M / (N * N * T);
}

int algoritmoMetropolis(int red[N][N], 
    double *sumaEnergia, double *sumaEnergiacuadrada, 
    double *sumaMagnetizacionSuperior, double *sumaMagnetizacionInferior, 
    int *configuracionesCambiadas, double *sumadensidadpositivo, 
    double *sumadensidadnegativo, double *sumaMagnetizacion2, double T) {

    double energias[PASOS];
    double energia_actual;

    int densidadpositivo = 0.0;
    int densidadnegativo = 0.0;
    int fila = N / 4;

    double energiaConfiguracion = 0.0;

    FILE *file = fopen("configuraciones.txt", "w");
    if (file == NULL) {
        printf("Error al abrir el fichero.\n");
        return 1;
    }

    FILE *energiaFile = fopen("energia_pmontecarlo.txt", "w");
    if (energiaFile == NULL) {
        printf("Error al abrir el fichero de energía.\n");
        return 1;
    }

    energia_actual = calcularEnergia_Total_Inicial(red);

    for (int j = 0; j < PASOS; j++) {
        for (int i = 0; i < N*N; i++) {
            int j1 = rand() % N;
            int i1 = rand() % (N-2) + 1;

            int deltaI = 0, deltaJ = 0;

            if (i1 == 1) {
                int direccion = rand() % 3;
                if (direccion == 0) deltaJ = 1;
                else if (direccion == 1) deltaJ = -1;
                else deltaI = 1;
            } else if (i1 == N - 2) {
                int direccion = rand() % 3;
                if (direccion == 0) deltaJ = -1;
                else if (direccion == 1) deltaJ = 1;
                else deltaI = -1;
            } else {
                int direccion = rand() % 4;
                if (direccion == 0) deltaJ = 1;
                else if (direccion == 1) deltaJ = -1;
                else if (direccion == 2) deltaI = 1;
                else if (direccion == 3) deltaI = -1;
            }

            int i2 = i1 + deltaI;
            int j2 = j1 + deltaJ;

            if (j2 < 0) j2 = N - 1;
            else if (j2 >= N) j2 = 0;

            if (red[i1][j1] == red[i2][j2]) continue;

            int energiaAntes = Energia_local(red, i1, j1, i2, j2);

            int temp = red[i1][j1];
            red[i1][j1] = red[i2][j2];
            red[i2][j2] = temp;

            int energiaDespues = Energia_local(red, i1, j1, i2, j2);

            int difE = energiaDespues - energiaAntes;

            double probabilidad = (difE > 0) ? exp(-difE / T) : 1.0;
            double r = (double)rand() / RAND_MAX;

            if (r > probabilidad) {
                temp = red[i1][j1];
                red[i1][j1] = red[i2][j2];
                red[i2][j2] = temp;
            } else {
                energia_actual += difE;
            }
        }

        energias[j] = energia_actual;
        fprintf(energiaFile, "%.6f %d\n", energia_actual, j + 1);
        Guardar_Red(file, red);

        if ((j + 1) % 100 == 0) {
            double magnetizacionSuperior = magnetizacion_mitad_superior(red);
            double magnetizacionInferior = magnetizacion_mitad_inferior(red);

            *sumaMagnetizacionSuperior += magnetizacionSuperior;
            *sumaMagnetizacionInferior += magnetizacionInferior;

            calcularDensidadesFila(red, fila, &densidadpositivo, &densidadnegativo);

            *sumadensidadpositivo += densidadpositivo;
            *sumadensidadnegativo += densidadnegativo;

            energiaConfiguracion = energias[j];
            *sumaEnergia += energiaConfiguracion;
            *sumaEnergiacuadrada += energiaConfiguracion * energiaConfiguracion;
            *sumaMagnetizacion2 += magnetizacionSuperior * magnetizacionSuperior;
            (*configuracionesCambiadas)++;
        }
    }

    fclose(file);
    fclose(energiaFile);
    return 0;
}