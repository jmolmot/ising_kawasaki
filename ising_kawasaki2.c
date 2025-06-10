#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define TAM 64
#define PASOS_MC 100000

double TEMP_INI = 0.25;
double TEMP_FIN = 4.0;
double TEMP_STEP = 0.25;

// --- Prototipos ---
void configurar_mitad(int matriz[TAM][TAM]);
void guardar_matriz(FILE *f, int matriz[TAM][TAM]);
double energia_local(int matriz[TAM][TAM], int x1, int y1, int x2, int y2);
double energia_total(int matriz[TAM][TAM]);
double mag_sup(int matriz[TAM][TAM]);
double mag_inf(int matriz[TAM][TAM]);
void densidades_fila(int matriz[TAM][TAM], int fila, int *pos, int *neg);
void calor_especifico(double sumaE, double sumaE2, int n, double T, double *cv);
void susceptibilidad(double sumaM, double sumaM2, int n, double T, double *chi);
int kawasaki(int matriz[TAM][TAM], double *sumaE, double *sumaE2, double *sumaMagSup, double *sumaMagInf, int *nCambios, double *sumaPos, double *sumaNeg, double *sumaMag2, double T);

int main() {
    int red[TAM][TAM];
    srand(time(NULL));
    const char* archivo_tiempo = "tiempos.txt";
    clock_t t_ini = clock(), t_fin;
    double t_total;

    char fname_sup[64], fname_inf[64], fname_pos[64], fname_neg[64], fname_cv[64], fname_chi[64];
    sprintf(fname_sup, "mag_sup%d.txt", TAM);
    sprintf(fname_inf, "mag_inf%d.txt", TAM);
    sprintf(fname_pos, "dens_pos%d.txt", TAM);
    sprintf(fname_neg, "dens_neg%d.txt", TAM);
    sprintf(fname_cv, "cv%d.txt", TAM);
    sprintf(fname_chi, "chi%d.txt", TAM);

    FILE *fsup = fopen(fname_sup, "w");
    FILE *finf = fopen(fname_inf, "w");
    FILE *fpos = fopen(fname_pos, "w");
    FILE *fneg = fopen(fname_neg, "w");
    FILE *fcv = fopen(fname_cv, "w");
    FILE *fchi = fopen(fname_chi, "w");

    for (double T = TEMP_INI; T <= TEMP_FIN + 1e-8; T += TEMP_STEP) {
        configurar_mitad(red);

        double sumaSup = 0.0, sumaInf = 0.0, sumaE = 0.0, sumaE2 = 0.0, sumaPos = 0.0, sumaNeg = 0.0, sumaMag2 = 0.0;
        int cambios = 0;
        double cv, chi;

        kawasaki(red, &sumaE, &sumaE2, &sumaSup, &sumaInf, &cambios, &sumaPos, &sumaNeg, &sumaMag2, T);

        double promSup = sumaSup / cambios;
        double promInf = sumaInf / cambios;
        double promPos = sumaPos / (cambios * TAM);
        double promNeg = sumaNeg / (cambios * TAM);

        calor_especifico(sumaE, sumaE2, cambios, T, &cv);
        susceptibilidad(sumaSup, sumaMag2, cambios, T, &chi);

        fprintf(fsup, "%.6f %.2f\n", promSup, T);
        fprintf(finf, "%.6f %.2f\n", promInf, T);
        fprintf(fneg, "%.6f %.2f\n", promNeg, T);
        fprintf(fpos, "%.6f %.2f\n", promPos, T);
        fprintf(fcv, "%.6f %.2f\n", cv, T);
        fprintf(fchi, "%.10f %.2f\n", chi, T);

        printf("T=%.2f OK\n", T);
    }

    fclose(fsup); fclose(finf); fclose(fpos); fclose(fneg); fclose(fcv); fclose(fchi);

    t_fin = clock();
    t_total = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    FILE *ft = fopen(archivo_tiempo, "a");
    if (ft != NULL) {
        fprintf(ft, "%d\t%.6f\n", TAM, t_total);
        fclose(ft);
    }
    return 0;
}

// --- FUNCIONES PERSONALIZADAS ---

void configurar_mitad(int matriz[TAM][TAM]) {
    int mitad = ((TAM-2) * TAM) / 2, pos = 0;
    for (int i = 0; i < TAM; ++i) {
        for (int j = 0; j < TAM; ++j) {
            if (i == 0) matriz[i][j] = -1;
            else if (i == TAM-1) matriz[i][j] = 1;
            else if (pos < mitad) { matriz[i][j] = -1; ++pos; }
            else { matriz[i][j] = 1; ++pos; }
        }
    }
}

void guardar_matriz(FILE *f, int matriz[TAM][TAM]) {
    for (int i = 0; i < TAM; ++i) {
        for (int j = 0; j < TAM; ++j) {
            fprintf(f, "%d", matriz[i][j]);
            if (j < TAM - 1) fprintf(f, ", ");
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n");
}

double energia_local(int matriz[TAM][TAM], int x1, int y1, int x2, int y2) {
    int E = 0;
    int arr1 = (x1 == 1) ? -1 : matriz[x1-1][y1];
    int ab1 = (x1 == TAM-2) ? 1 : matriz[x1+1][y1];
    int izq1 = (y1 == 0) ? matriz[x1][TAM-1] : matriz[x1][y1-1];
    int der1 = (y1 == TAM-1) ? matriz[x1][0] : matriz[x1][y1+1];

    int arr2 = (x2 == 1) ? -1 : matriz[x2-1][y2];
    int ab2 = (x2 == TAM-2) ? 1 : matriz[x2+1][y2];
    int izq2 = (y2 == 0) ? matriz[x2][TAM-1] : matriz[x2][y2-1];
    int der2 = (y2 == TAM-1) ? matriz[x2][0] : matriz[x2][y2+1];

    if ((x1 == x2 && (y1 + 1) % TAM == y2)) { der1 = 0; izq2 = 0; }
    if ((x1 == x2 && (y1 - 1 + TAM) % TAM == y2)) { izq1 = 0; der2 = 0; }
    if ((x1 + 1 == x2 && y1 == y2)) { arr2 = 0; ab1 = 0; }
    if ((x1 - 1 == x2 && y1 == y2)) { ab2 = 0; arr1 = 0; }

    E = -matriz[x1][y1] * (arr1 + ab1 + izq1 + der1);
    E += -matriz[x2][y2] * (arr2 + ab2 + izq2 + der2);
    return E;
}

double energia_total(int matriz[TAM][TAM]) {
    double energia = 0.0;
    for (int i = 1; i < TAM-1; ++i) {
        for (int j = 0; j < TAM; ++j) {
            int suma = matriz[i-1][j] + matriz[i+1][j] + matriz[i][(j+1)%TAM] + matriz[i][(j-1+TAM)%TAM];
            energia += matriz[i][j] * suma;
        }
    }
    return -0.5 * energia;
}

double mag_sup(int matriz[TAM][TAM]) {
    int suma = 0;
    for (int i = 0; i < TAM/2; ++i)
        for (int j = 0; j < TAM; ++j)
            suma += matriz[i][j];
    return fabs((double)suma) / ((TAM/2)*TAM);
}

double mag_inf(int matriz[TAM][TAM]) {
    int suma = 0;
    for (int i = TAM/2; i < TAM; ++i)
        for (int j = 0; j < TAM; ++j)
            suma += matriz[i][j];
    return fabs((double)suma) / ((TAM/2)*TAM);
}

void densidades_fila(int matriz[TAM][TAM], int fila, int *pos, int *neg) {
    int p = 0, n = 0;
    for (int j = 0; j < TAM; ++j) {
        if (matriz[fila][j] == 1) ++p;
        else if (matriz[fila][j] == -1) ++n;
    }
    *pos = p; *neg = n;
}

void calor_especifico(double sumaE, double sumaE2, int n, double T, double *cv) {
    double promE = sumaE / n;
    double promE2 = sumaE2 / n;
    double varE = promE2 - (promE * promE);
    *cv = varE / ((TAM-2)*TAM*T*T);
}

void susceptibilidad(double sumaM, double sumaM2, int n, double T, double *chi) {
    double promM = sumaM / n;
    double promM2 = sumaM2 / n;
    double varM = promM2 - (promM * promM);
    *chi = varM / (TAM * TAM * T);
}

int kawasaki(int matriz[TAM][TAM], double *sumaE, double *sumaE2, double *sumaMagSup, double *sumaMagInf, int *nCambios, double *sumaPos, double *sumaNeg, double *sumaMag2, double T) {
    double energia_actual = energia_total(matriz);
    int fila = TAM / 4;
    FILE *fmat = fopen("matrices.txt", "w");
    FILE *fener = fopen("energia_mc.txt", "w");
    if (!fmat || !fener) return 1;

    for (int paso = 0; paso < PASOS_MC; ++paso) {
        for (int k = 0; k < TAM*TAM; ++k) {
            int y1 = rand() % TAM;
            int x1 = rand() % (TAM-2) + 1;
            int dx = 0, dy = 0;
            if (x1 == 1) {
                int dir = rand() % 3;
                if (dir == 0) dy = 1;
                else if (dir == 1) dy = -1;
                else dx = 1;
            } else if (x1 == TAM-2) {
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
            if (y2 < 0) y2 = TAM-1;
            else if (y2 >= TAM) y2 = 0;
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