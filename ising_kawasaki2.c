#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

double T_ini = 0.25;
double T_fin = 5.0;
double T_step = 0.25;
int N;
int N_ini = 16;
int N_fin = 128;
int N_step = 16;


#define pasos 50000

void configurar_mitad(int matriz[N][N]);
void guardar_matriz(FILE *f, int matriz[N][N]);
double energia_local(int matriz[N][N], int x1, int y1, int x2, int y2);
double energia_total(int matriz[N][N]);
double mag_sup(int matriz[N][N]);
double mag_inf(int matriz[N][N]);
void densidades_fila(int matriz[N][N], int fila, int *pos, int *neg);
void calor_especifico(double sumaE, double sumaE2, int n, double T, double *cv);
void susceptibilidad(double sumaM, double sumaM2, int n, double T, double *chi);
int kawasaki(int matriz[N][N], double *sumaE, double *sumaE2, double *sumaMagSup, double *sumaMagInf, int *nCambios, double *sumaPos, double *sumaNeg, double *sumaMag2, double T);

int main() {
    srand(time(NULL)); // Semilla para números aleatorios
    // Puedes crear la carpeta desde C si quieres, pero aquí asumimos que ya existe
    for (int N_actual = N_ini; N_actual <= N_fin; N_actual += N_step) {
        N = N_actual;
        int red[N][N];

        char fname_sup[128], fname_inf[128], fname_pos[128], fname_neg[128], fname_cv[128], fname_chi[128], fname_energia[128];
        sprintf(fname_sup, "archivos/promedio_magnetizacionsuperior%d.txt", N);
        sprintf(fname_inf, "archivos/promedio_magnetizacioninferior%d.txt", N);
        sprintf(fname_pos, "archivos/promedio_densidadpositivo%d.txt", N);
        sprintf(fname_neg, "archivos/promedio_densidadnegativo%d.txt", N);
        sprintf(fname_cv, "archivos/filecv%d.txt", N);
        sprintf(fname_chi, "archivos/susceptibilidad%d.txt", N);
        sprintf(fname_energia, "archivos/promedio_energia%d.txt", N);

        FILE *fsup = fopen(fname_sup, "w");
        FILE *finf = fopen(fname_inf, "w");
        FILE *fpos = fopen(fname_pos, "w");
        FILE *fneg = fopen(fname_neg, "w");
        FILE *fcv = fopen(fname_cv, "w");
        FILE *fchi = fopen(fname_chi, "w");
        FILE *fenergia = fopen(fname_energia, "w");

        clock_t t_ini_clock = clock();

        for (double T = T_ini; T <= T_fin + 1e-8; T += T_step) {
            configurar_mitad(red);

            double sumaSup = 0.0, sumaInf = 0.0, sumaE = 0.0, sumaE2 = 0.0, sumaPos = 0.0, sumaNeg = 0.0, sumaMag2 = 0.0;
            int cambios = 0;
            double cv, chi;

            kawasaki(red, &sumaE, &sumaE2, &sumaSup, &sumaInf, &cambios, &sumaPos, &sumaNeg, &sumaMag2, T);

            double promSup = sumaSup / cambios;
            double promInf = sumaInf / cambios;
            double promPos = sumaPos / (cambios * N);
            double promNeg = sumaNeg / (cambios * N);
            double energia_media = sumaE / cambios;

            calor_especifico(sumaE, sumaE2, cambios, T, &cv);
            susceptibilidad(sumaSup, sumaMag2, cambios, T, &chi);

            fprintf(fsup, "%g %g\n", promSup, T);
            fprintf(finf, "%g %g\n", promInf, T);
            fprintf(fneg, "%g %g\n", promNeg, T);
            fprintf(fpos, "%g %g\n", promPos, T);
            fprintf(fcv, "%g %g\n", cv, T);
            fprintf(fchi, "%g %g\n", chi, T);
            fprintf(fenergia, "%g %g\n", energia_media, T);

            printf("N=%d T=%.2f Terminado\n", N, T);
        }

        fclose(fsup); fclose(finf); fclose(fpos); fclose(fneg); fclose(fcv); fclose(fchi); fclose(fenergia);

        clock_t t_fin_clock = clock();
        double t_total = (double)(t_fin_clock - t_ini_clock) / CLOCKS_PER_SEC;

        FILE *ft = fopen("archivos/tiemposN.txt", "a");
        if (ft != NULL) {
            fprintf(ft, "%d\t%.6f\n", N, t_total);
            fclose(ft);
        }

        printf("Bucle de temperatura para N=%d completado en %.g segundos.\n", N, t_total);
    }
    return 0;
}

void configurar_mitad(int matriz[N][N]) {
    int mitad = ((N-2) * N) / 2, pos = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == 0) matriz[i][j] = -1;
            else if (i == N-1) matriz[i][j] = 1;
            else if (pos < mitad) { matriz[i][j] = -1; ++pos; }
            else { matriz[i][j] = 1; ++pos; }
        }
    }
}

void guardar_matriz(FILE *f, int matriz[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(f, "%d", matriz[i][j]);
            if (j < N - 1) fprintf(f, ", ");
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

double energia_local(int matriz[N][N], int x1, int y1, int x2, int y2) {
    int E = 0;
    int arr1 = (x1 == 1) ? -1 : matriz[x1-1][y1];
    int ab1 = (x1 == N-2) ? 1 : matriz[x1+1][y1];
    int izq1 = (y1 == 0) ? matriz[x1][N-1] : matriz[x1][y1-1];
    int der1 = (y1 == N-1) ? matriz[x1][0] : matriz[x1][y1+1];

    int arr2 = (x2 == 1) ? -1 : matriz[x2-1][y2];
    int ab2 = (x2 == N-2) ? 1 : matriz[x2+1][y2];
    int izq2 = (y2 == 0) ? matriz[x2][N-1] : matriz[x2][y2-1];
    int der2 = (y2 == N-1) ? matriz[x2][0] : matriz[x2][y2+1];

    if ((x1 == x2 && (y1 + 1) % N == y2)) { der1 = 0; izq2 = 0; }
    if ((x1 == x2 && (y1 - 1 + N) % N == y2)) { izq1 = 0; der2 = 0; }
    if ((x1 + 1 == x2 && y1 == y2)) { arr2 = 0; ab1 = 0; }
    if ((x1 - 1 == x2 && y1 == y2)) { ab2 = 0; arr1 = 0; }

    E = -matriz[x1][y1] * (arr1 + ab1 + izq1 + der1);
    E += -matriz[x2][y2] * (arr2 + ab2 + izq2 + der2);
    return E;
}

double energia_total(int matriz[N][N]) {
    double energia = 0.0;
    for (int i = 1; i < N-1; ++i) {
        for (int j = 0; j < N; ++j) {
            int suma = matriz[i-1][j] + matriz[i+1][j] + matriz[i][(j+1)%N] + matriz[i][(j-1+N)%N];
            energia += matriz[i][j] * suma;
        }
    }
    return -0.5 * energia;
}

double mag_sup(int matriz[N][N]) {
    int suma = 0;
    for (int i = 0; i < N/2; ++i)
        for (int j = 0; j < N; ++j)
            suma += matriz[i][j];
    return fabs((double)suma) / ((N/2)*N);
}

double mag_inf(int matriz[N][N]) {
    int suma = 0;
    for (int i = N/2; i < N; ++i)
        for (int j = 0; j < N; ++j)
            suma += matriz[i][j];
    return fabs((double)suma) / ((N/2)*N);
}

void densidades_fila(int matriz[N][N], int fila, int *pos, int *neg) {
    int p = 0, n = 0;
    for (int j = 0; j < N; ++j) {
        if (matriz[fila][j] == 1) ++p;
        else if (matriz[fila][j] == -1) ++n;
    }
    *pos = p; *neg = n;
}

void calor_especifico(double sumaE, double sumaE2, int n, double T, double *cv) {
    double promE = sumaE / n;
    double promE2 = sumaE2 / n;
    double varE = promE2 - (promE * promE);
    *cv = varE / ((N-2)*N*T*T);
}

void susceptibilidad(double sumaM, double sumaM2, int n, double T, double *chi) {
    double promM = sumaM / n;
    double promM2 = sumaM2 / n;
    double varM = promM2 - (promM * promM);
    *chi = varM / (N * N * T);
}

int kawasaki(int matriz[N][N], double *sumaE, double *sumaE2, double *sumaMagSup, double *sumaMagInf, int *nCambios, double *sumaPos, double *sumaNeg, double *sumaMag2, double T) {
    double energia_actual = energia_total(matriz);
    int fila = N / 4;
    FILE *fmat = fopen("matrices.txt", "w");
    FILE *fener = fopen("energia_mc.txt", "w");
    if (!fmat || !fener) return 1;

    for (int paso = 0; paso < pasos; ++paso) {
        for (int k = 0; k < N*N; ++k) {
            int y1 = rand() % N;
            int x1 = rand() % (N-2) + 1;
            int dx = 0, dy = 0;
            if (x1 == 1) {
                int dir = rand() % 3;
                if (dir == 0) dy = 1;
                else if (dir == 1) dy = -1;
                else dx = 1;
            } else if (x1 == N-2) {
                int dir = rand() % 3;
                if (dir == 0) dy = -1;
                else if (dir == 1) dy = 1;
                else dx = -1;
            } else {
                int dir = rand() % 4;
                if (dir == 0) dy = 1;
                else if (dir == 1) dy = -1;
                else if (dir == 2) dx = 1;
                else dx = -1;
            }
            int x2 = x1 + dx, y2 = y1 + dy;
            if (y2 < 0) y2 = N-1;
            else if (y2 >= N) y2 = 0;
            if (matriz[x1][y1] == matriz[x2][y2]) continue;
            int eAntes = energia_local(matriz, x1, y1, x2, y2);
            int tmp = matriz[x1][y1];
            matriz[x1][y1] = matriz[x2][y2];
            matriz[x2][y2] = tmp;
            int eDespues = energia_local(matriz, x1, y1, x2, y2);
            int dE = eDespues - eAntes;
            double p = (dE > 0) ? exp(-dE / T) : 1.0;
            double r = (double)rand() / RAND_MAX;
            if (r > p) {
                tmp = matriz[x1][y1];
                matriz[x1][y1] = matriz[x2][y2];
                matriz[x2][y2] = tmp;
            } else {
                energia_actual += dE;
            }
        }
        fprintf(fener, "%.6f %d\n", energia_actual, paso + 1);
        guardar_matriz(fmat, matriz);
        if ((paso + 1) % 100 == 0) {
            double msup = mag_sup(matriz);
            double minf = mag_inf(matriz);
            *sumaMagSup += msup;
            *sumaMagInf += minf;
            int pos, neg;
            densidades_fila(matriz, fila, &pos, &neg);
            *sumaPos += pos;
            *sumaNeg += neg;
            *sumaE += energia_actual;
            *sumaE2 += energia_actual * energia_actual;
            *sumaMag2 += msup * msup;
            ++(*nCambios);
        }
    }
    fclose(fmat);
    fclose(fener);
    return 0;
}