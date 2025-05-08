#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define N_inicial 2   // Tamaño de la matriz
#define N_final 32
#define T 1.50 // Temperatura inicial
#define T_final 3.0 // Temperatura final
#define T_incremento 0.1 // Incremento de temperatura 
#define P 10000 //Número de pasos de Monte Carlo
#define pasos_equilibrio 1000 // Número de pasos de equilibrado
#define pasos_medidas 1000 // Número de pasos de medida

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

int main()
{
    srand(time(NULL)); // Inicializa semilla de aleatoriedad

    FILE*fichero = fopen("resultados.txt", "w");

    for(int N=N_inicial; N<=N_final; N*=2)
    {

        int** spin = crear_matriz(N);
        int** spin_inicial = crear_matriz(N); // Matriz para guardar el estado inicial
        if (spin == NULL)
    {
        printf("Error al reservar memoria para la matriz.\n");
        return 1;
    }
        fprintf(fichero, "\n# ===== RESULTADOS PARA N = %d =====\n", N);
        
        inicializar_espin(spin, N);
        copiar_matriz(spin, spin_inicial, N); // Copia el estado inicial

        fprintf(fichero, "#T\tE_por_espin\tDensidad_y\tCalor_especifico\tSusceptibilidad\n");

        for(double T_actual=T; T_actual<=T_final+T_incremento; T_actual+=T_incremento)
        {
            copiar_matriz(spin_inicial, spin, N); // Copia el estado inicial en cada iteración
            paso_montecarlo(spin, P, N, T_actual);

            double E = energia_total(spin, N)/(N*N);
            double densidad_y = densidad_media_y(spin, N);
            double cN = calor_especifico(spin, N, T_actual, pasos_equilibrio, pasos_medidas);
            double chi = susceptibilidad_magnetica(spin, N, T_actual, pasos_equilibrio, pasos_medidas);

            fprintf(fichero, "%g\t%g\t%g\t%g\t%g\n", T_actual, E, densidad_y, cN, chi);

            // Escribir matriz de espines
            fprintf(fichero, "Matriz de espines (T=%g):\n", T_actual);
            imprimir_matriz(spin, N, fichero);
            fprintf(fichero, "\n");
        
            printf("T=%g hecho.\n", T_actual);
        }
        liberar_matriz(spin, N);
    }

    fclose(fichero);
    printf("Resultados guardados en resultados.txt\n");

    /*

    //Matriz antes del cambio
    printf("Matriz de espines inicial:\n");
    imprimir_matriz(spin, N);

    printf("\n");

    printf("Matriz despues del paso Monte Carlo:\n");
    paso_montecarlo(spin, P, N, T);
    imprimir_matriz(spin, N);

    printf("\n");

    double energia = energia_total(spin, N);
    printf("Energia total: %f\n", energia);
    printf("Energia por espin: %f\n", energia / (N * N));

    printf("\n");

    double rho_y = densidad_media_y(spin, N);
    printf("Densidad media de partículas (+1) en dirección y: %g\n", rho_y);

    printf("\n");

    double cN= calor_especifico(spin, N, T, pasos_eq, pasos_medida);
    printf("Calor específico: %g\n", cN);

    printf("\n");

    double m_sup = magnetizacion_superior(spin, N);
    printf("Magnetización superior: %g\n", m_sup);
    double m_inf = magnetizacion_inferior(spin, N);
    printf("Magnetización inferior: %g\n", m_inf);
    double promedio_magnetizacion = (m_sup + m_inf) / 2.0;
    printf("Promedio de magnetización: %g\n", promedio_magnetizacion);

    printf("\n");
    double m = (m_sup + m_inf) / 2.0;

    double chi = susceptibilidad_magnetica(spin, N, T, pasos_eq, pasos_medida);
    printf("Susceptibilidad magnética: %g\n", chi);

    */

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

    for (int i = 0; i < n; i++)
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

    printf("Red de %d x %d creada\n", n, n);
    return matriz;
}

void copiar_matriz(int** matriz, int** matriz_copia, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matriz_copia[i][j] = matriz[i][j];
}

void inicializar_espin(int** matriz, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matriz[i][j] = (rand() % 2) ? 1 : -1;
}

void imprimir_matriz(int** matriz, int n, FILE*fichero)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            fprintf(fichero, "%c", matriz[i][j] == 1 ? '+' : '-');
        fprintf(fichero, "\n");
    }
}

void liberar_matriz(int** matriz_liberada, int n)
{
    if (matriz_liberada == NULL) return;
    for (int i = 0; i < n; i++)
        free(matriz_liberada[i]);
    free(matriz_liberada);
    printf("Memoria liberada correctamente.\n");
}

int indice_periodico(int i, int n)
{
    if( i>=n)
    {
        return 0;
    }
    if( i<0)
    {
        return n-1;
    }
    return i;
}

int delta_E(int** matriz, int n, int i1, int j1, int i2, int j2)
{

    int s1 = matriz[i1][j1];
    int s2 = matriz[i2][j2];
    int delta = 0;

    int vecinos[4][2] = {{-1,0}, {1,0}, {0,-1}, {0,1}};

    for(int k = 0; k < 4; k++)
    {
        int ni1 = indice_periodico(i1 + vecinos[k][0], n);
        int nj1 = indice_periodico(j1 + vecinos[k][1], n);

        int ni2 = indice_periodico(i2 + vecinos[k][0], n);
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
    int i1, j1, i2, j2;

    // Selección aleatoria de dos espines
    do
    {
        i1 = rand() % n;
        j1 = rand() % n;

        int direccion = rand() % 4; // Dirección aleatoria

        if (direccion == 0)
        {
            i2 = indice_periodico(i1 - 1, n);
            j2 = j1;
        }  
        else if (direccion == 1)
        {
            i2 = indice_periodico(i1 + 1, n);
            j2 = j1;
        }
        else if (direccion == 2)
        {
            i2 = i1;
            j2 = indice_periodico(j1 - 1, n);
        }
        else
        {
            i2 = i1;
            j2 = indice_periodico(j1 + 1, n);
        }
    
    }while (matriz[i1][j1] == matriz[i2][j2]); // Asegurarse de que los espines sean diferentes

    int dE = delta_E(matriz, n, i1, j1, i2, j2);

    if(dE<0 || (rand()/(double)RAND_MAX) < exp(-dE/t))
    {
        // Intercambiar espines
        int aux = matriz[i1][j1];
        matriz[i1][j1] = matriz[i2][j2];
        matriz[i2][j2] = aux;

    }
}

void paso_montecarlo(int** matriz, int pasos, int n, double t)
{
    for (int p=0; p < pasos; p++)
    {
        for(int k=0; k < n*n; k++)
        {
            paso_kawasaki(matriz, n, t);
        }
    }
}

double energia_total(int** matriz, int n)
{
    double E = 0.0;
    
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            int s = matriz[i][j];
            
            int vecino_derecha= matriz[i][indice_periodico(j+1, n)];
            int vecino_abajo= matriz[indice_periodico(i+1, n)][j];

            E -= s * (vecino_derecha + vecino_abajo); 
        }
    }
    return E;
}

double densidad_media_y(int** matriz, int n)
{
    double suma = 0.0;
    for(int i=0; i<n; i++)
    {
        double cuenta_fila=0.0;
        for (int j=0; j<n; j++)
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

    for(int i=0; i<pasos_medida; i++)
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

    for(int i=0; i<n/2; i++)
    {
        for(int j=0; j<n; j++)
        {
            suma += matriz[i][j];
        }
    }

    return suma/(n*n/2.0);
}

double magnetizacion_inferior(int**matriz, int n)
{
    double suma = 0.0;

    for(int i=n/2; i<n; i++)
    {
        for(int j=0; j<n; j++)
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

    for(int i=0; i<pasos_medida; i++)
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