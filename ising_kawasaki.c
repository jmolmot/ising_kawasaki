#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

int N_inicial = 32;   // Tamaño de la matriz
int N_final = 32;
double T = 2.0; // Temperatura inicial
double T_final = 2.0; // Temperatura final
double T_incremento = 0.1; // Incremento de temperatura 
int pasos_equilibrio = 1000; // Número de pasos de equilibrado
int pasos_medidas = 100000; // Número de pasos de medida

// Prototipos de funciones
int** crear_matriz(int n);
void copiar_matriz(int** matriz, int** matriz_copia, int n);
void inicializar_espin(int** matriz, int n);
void imprimir_matriz(int** matriz, int n, FILE*fichero);
void liberar_matriz(int** matriz_liberada, int n);
int indice_periodico(int i, int n);
int delta_E(int** matriz, int n, int i1, int j1, int i2, int j2);
void paso_kawasaki(int** matriz, int n, double t);
void paso_montecarlo(int** matriz, int p, int n, double t);
double energia_total(int** matriz, int n);
double densidad_media_y(int** matriz, int n);
double calor_especifico(int**matriz, int n, double t, int pasos_eq, int pasos_medida);
double magnetizacion_superior(int**matriz, int n);
double magnetizacion_inferior(int**matriz, int n);
double susceptibilidad_magnetica(int**matriz, int n, double t, int pasos_eq, int pasos_medida);
void histograma_spins_fila(int** matriz, int n, int* histo);

int main()
{
    srand(time(NULL));
    clock_t inicio = clock();

    FILE* fmag = fopen("magnetizaciones.txt", "w");
    FILE* fener = fopen("energias.txt", "w");
    FILE* fdens = fopen("densidades.txt", "w");
    FILE* fcv = fopen("cv.txt", "w");
    FILE* fchi = fopen("chi.txt", "w");
    FILE* fmat = fopen("matrices.txt", "w");

    for(int N=N_inicial; N<=N_final; N*=2)
    {
        int** spin = crear_matriz(N);
        int** spin_inicial = crear_matriz(N);
        inicializar_espin(spin, N);
        copiar_matriz(spin, spin_inicial, N);

        for(double T_actual=T; T_actual<=T_final; T_actual+=T_incremento)
        {
            copiar_matriz(spin_inicial, spin, N);
            paso_montecarlo(spin, pasos_equilibrio, N, T_actual);

            double suma_E = 0.0, suma_E2 = 0.0;
            double suma_M = 0.0, suma_M2 = 0.0;
            double suma_msup = 0.0, suma_minf = 0.0;
            double suma_dens = 0.0;
            double suma_msup_chi= 0.0, suma_minf_chi = 0.0;

            for(int paso=0; paso<pasos_medidas; ++paso)
            {
                paso_montecarlo(spin, 1, N, T_actual);

                // Magnetización superior e inferior
                double m_sup = magnetizacion_superior(spin, N);
                double m_inf = magnetizacion_inferior(spin, N);
                suma_msup += m_sup;
                suma_minf += m_inf;

                suma_msup_chi += m_sup;
                suma_minf_chi += m_inf;

                // Energía total
                double E = energia_total(spin, N);
                suma_E += E;
                suma_E2 += E*E;


                // Densidad media
                suma_dens += densidad_media_y(spin, N);

                // Guardado en archivos
                fprintf(fmag, "%g\t%g\t%g\t%d\t%g\n", T_actual, m_sup, m_inf, paso, (m_sup+m_inf)/2.0);
                fprintf(fener, "%g\t%g\t%d\n", T_actual, E, paso);
                fprintf(fdens, "%g\t%d", T_actual, paso);
                for(int i=0; i<N; ++i)
                {
                    double dens_col = 0.0;
                    for(int j=0; j<N; ++j)
                        if(spin[i][j] == 1) dens_col += 1.0;
                    dens_col /= N;
                    fprintf(fdens, "\t%g", dens_col);
                }
                fprintf(fdens, "\n");

                // Calor específico y susceptibilidad magnética instantáneos
                double prom_E = suma_E / (paso+1);
                double prom_E2 = suma_E2 / (paso+1);
                double cv = (prom_E2 - prom_E*prom_E) / (T_actual*T_actual*N*N);
                fprintf(fcv, "%g\t%d\t%g\n", T_actual, paso, cv);

                double prom_msup_chi = suma_msup_chi / (paso+1);
                double prom_minf_chi = suma_minf_chi / (paso+1);
                double chi = (prom_minf_chi-prom_msup_chi*prom_msup_chi) / (T_actual*N*N);
                fprintf(fchi, "%g\t%d\t%g\n", T_actual, paso, chi);

                // Guardar la matriz completa en matrices.txt
                 for(int i=0; i<N; ++i)
                {
                    for(int j=0; j<N; ++j) 
                    {
                        fprintf(fmat, "%d", spin[i][j]);
                        if (j < N-1) fprintf(fmat, ", ");
                    }
                    fprintf(fmat, "\n");
                }
                fprintf(fmat, "\n");
            }
        }
        int* histo = (int*)calloc(N+1, sizeof(int));
        histograma_spins_fila(spin, N, histo);

        FILE* fhisto = fopen("histograma.txt", "w");
        fprintf(fhisto, "#spins+1\tcuentas   \n");
        for(int k=0; k<=N; ++k)
        {
            fprintf(fhisto, "%d\t%d\n", k, histo[k]);
        }
        fclose(fhisto);
        free(histo);

        liberar_matriz(spin, N);
        liberar_matriz(spin_inicial, N);
    }

    fclose(fmag);
    fclose(fener);
    fclose(fdens);
    fclose(fcv);
    fclose(fchi);
    fclose(fmat);

    clock_t fin = clock();
    double tiempo_total = (double)(fin - inicio) / CLOCKS_PER_SEC;
    printf("Tiempo total de ejecución: %.2f segundos\n", tiempo_total);

    return 0;
}

int** crear_matriz(int n)
{
    int** matriz = (int**)malloc(n * sizeof(int*));
    if (matriz == NULL)
    {
        perror("Error al asignar memoria para la matriz");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n; ++i)
    {
        matriz[i] = (int*)malloc(n * sizeof(int));
        if (matriz[i] == NULL)
        {
            perror("Error al asignar memoria para una fila");
            for (int j = 0; j < i; j++)
                free(matriz[j]);
            free(matriz);
            exit(EXIT_FAILURE);
        }
    }

    return matriz;
}

void copiar_matriz(int** matriz, int** matriz_copia, int n)
{
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            matriz_copia[i][j] = matriz[i][j];
}

void inicializar_espin(int** matriz, int n)
{
    for(int j=0; j<n; ++j)
    {
        matriz[0][j]=1;
        matriz[n-1][j]=-1;
    }
    for (int i = 1; i < n-1; ++i)
        for (int j = 0; j < n; ++j)
            matriz[i][j] = (rand() % 2) ? 1 : -1;
}

void imprimir_matriz(int** matriz, int n, FILE*fichero)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
            fprintf(fichero, "%c", matriz[i][j] == 1 ? '+' : '-');
        fprintf(fichero, "\n");
    }
}

void liberar_matriz(int** matriz_liberada, int n)
{
    for (int i = 0; i < n; ++i)
        free(matriz_liberada[i]);
    free(matriz_liberada);
}

int indice_periodico(int i, int n)
{
	int num;
    if( i>=n)
    {
         num = 0;
    }
    else if( i<0)
    {
        num = n-1;
    } 
    else
    {
        num = i;
    }
    return num;
}

int delta_E(int** matriz, int n, int i1, int j1, int i2, int j2)
{
    int s1 = matriz[i1][j1];
    int s2 = matriz[i2][j2];
    int delta = 0;

    int vecinos[4][2] = {{-1,0}, {1,0}, {0,-1}, {0,1}};

    for(int k = 0; k < 4; ++k)
    {
        int ni1 = indice_periodico_filas(i1 + vecinos[k][0], n);
        int nj1 = indice_periodico(j1 + vecinos[k][1], n);

        int ni2 = indice_periodico_filas(i2 + vecinos[k][0], n);
        int nj2 = indice_periodico(j2 + vecinos[k][1], n);

        if (ni1 != i2 || nj1 != j2)
            delta += 2 * s1 * matriz[ni1][nj1];

        if (ni2 != i1 || nj2 != j1)
            delta += 2 * s2 * matriz[ni2][nj2];
    }

    return delta;
}

void paso_kawasaki(int** matriz, int n, double t)
{
    int i = 1 + rand() % (n - 2); // Solo filas internas
    int j = rand() % n;
    int s = matriz[i][j];

    int vecinos[4][2] = {{-1,0}, {1,0}, {0,-1}, {0,1}};
    int candidatos[4][2];
    int num_candidatos = 0;

    for (int k = 0; k < 4; ++k)
    {
        int ni = i + vecinos[k][0];
        int nj = j + vecinos[k][1];

        // Asegurarse de que ni está en [1, n-2] y nj en [0, n-1]
        if (ni >= 1 && ni < n - 1)
        {
            // Usar índice periódico solo en la dirección j
            nj = indice_periodico(nj, n);

            if (matriz[ni][nj] != s)
            {
                candidatos[num_candidatos][0] = ni;
                candidatos[num_candidatos][1] = nj;
                num_candidatos++;
            }
        }
    }

    if (num_candidatos == 0) return;

    int elegido = rand() % num_candidatos;
    int i2 = candidatos[elegido][0];
    int j2 = candidatos[elegido][1];

    int dE = delta_E(matriz, n, i, j, i2, j2);
    double probabilidad = exp(-dE / t);
    double r = (double)rand() / RAND_MAX;

    if (dE < 0 || r < probabilidad)
    {
        int aux = matriz[i][j];
        matriz[i][j] = matriz[i2][j2];
        matriz[i2][j2] = aux;
    }
}

void paso_montecarlo(int** matriz, int pasos, int n, double t)
{
    for (int p=0; p < pasos; ++p)
    {
        for(int k=0; k < n*n; ++k)
        {
            paso_kawasaki(matriz, n, t);
        }
    }
}

int indice_periodico_filas(int i, int n) {
    if (i == 0) return 2;           // La fila 0 solo se conecta con la 2
    if (i == n-1) return n-2;       // La fila n-1 solo se conecta con la n-2
    if (i < 0) return n-2;          // Para filas internas
    if (i >= n) return 1;           // Para filas internas
    return i;
}

double energia_total(int** matriz, int n)
{
    double E = 0.0;

    for(int i=0; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            int s = matriz[i][j];

            // Vecino a la derecha (condiciones periódicas en columnas)
            int vecino_derecha = matriz[i][indice_periodico(j+1, n)];

            // Vecino abajo (condiciones periódicas solo entre filas internas)
            int vecino_abajo = 0;
            int i_abajo = indice_periodico_filas(i+1, n);
            // Solo filas internas se conectan periódicamente entre sí
            if (i != n-1) {
                vecino_abajo = matriz[i_abajo][j];
            }
            // Si i == n-1 (última fila), vecino_abajo = 0

            E -= s * (vecino_derecha + vecino_abajo);
        }
    }
    return E;
}

double densidad_media_y(int** matriz, int n)
{
    double suma = 0.0;
    for(int i=0; i<n; ++i)
    {
        double cuenta_fila=0.0;
        for (int j=0; j<n; ++j)
        {
            if(matriz[i][j] == 1)
            {
                cuenta_fila += 1.0;
            }
        }
        double densidad_fila = cuenta_fila / n;
        suma += densidad_fila;
    }
    return suma / n;
}

double calor_especifico(int**matriz, int n, double t, int pasos_eq, int pasos_medida)
{
    paso_montecarlo(matriz, pasos_eq, n, t);

    double suma_E=0.0;
    double suma_E2= 0.0;

    for(int i=0; i<pasos_medida; ++i)
    {
        paso_montecarlo(matriz, 1, n, t);
        
        double E= energia_total(matriz, n);
        suma_E += E;
        suma_E2 += E*E;
    }

    double promedio_E = suma_E / pasos_medida;
    double promedio_E2 = suma_E2 / pasos_medida;

    double cN= (promedio_E2 - promedio_E*promedio_E) / (t*t * n*n);

    return cN;
}

double magnetizacion_superior(int**matriz, int n)
{
    double suma = 0.0;

    for(int i=0; i<n/2; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            suma += matriz[i][j];
        }
    }

    return suma/(n*n/2.0);
}

double magnetizacion_inferior(int**matriz, int n)
{
    double suma = 0.0;

    for(int i=n/2; i<n; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            suma += matriz[i][j];
        }
    }

    return suma/(n*n/2.0);
}

double susceptibilidad_magnetica(int**matriz, int n, double t, int pasos_eq, int pasos_medida)

{
    paso_montecarlo(matriz, pasos_eq, n, t);
    double suma_m=0.0;
    double suma_m2=0.0;

    for(int i=0; i<pasos_medida; ++i)
    {
        paso_montecarlo(matriz, 1, n, t);
        
        double m= magnetizacion_superior(matriz, n);
        suma_m += m;
        suma_m2 += m*m;
    }

    double promedio_m = suma_m / pasos_medida;
    double promedio_m2 = suma_m2 / pasos_medida;

    double susceptibilidad= (promedio_m2 - promedio_m*promedio_m) / (t * n * n);

    return susceptibilidad;
}

void histograma_spins_fila(int** matriz, int n, int* histo) 
{
    // Inicializa el histograma
    for (int k = 0; k <= n; ++k)
        histo[k] = 0;

    // Cuenta cuántos +1 hay en cada fila y actualiza el histograma
    for (int i = 0; i < n; ++i) {
        int cuenta = 0;
        for (int j = 0; j < n; ++j)
            if (matriz[i][j] == 1)
                cuenta++;
        histo[cuenta]++;
    }
}