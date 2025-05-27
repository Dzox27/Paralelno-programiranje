#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define MAX_ITER 1000
#define TOL 1e-6
#define NUM_TRIALS 10

// Sekvencijalna implementacija množenja matrice i vektora
void matrix_vector_multiplication_sequential(double **A, double *x, double *y, int N)
{
    // Iteracija kroz redove matrice
    for (int i = 0; i < N; i++)
    {
        y[i] = 0.0; // Inicijalizacija i-tog elementa rezultujućeg vektora

        // Množenje i-tog reda matrice sa vektorom x
        for (int j = 0; j < N; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Paralelna implementacija pomoću OpenMP-a
void matrix_vector_multiplication_parallel(double **A, double *x, double *y, int N)
{
// Svaka nit paralelno obrađuje različit red matrice
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        y[i] = 0.0; // Inicijalizacija i-tog elementa rezultujućeg vektora

        // Množenje i-tog reda matrice sa vektorom x
        for (int j = 0; j < N; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Sekvencijalna implementacija Euklidske norme (L2 norme) vektora
double vector_norm_sequential(double *x, int N)
{
    double suma = 0.0;

    // Sabiranje kvadrata svih elemenata vektora
    for (int i = 0; i < N; i++)
        suma += x[i] * x[i];

    // Vraća kvadratni korijen sume (norma vektora)
    return sqrt(suma);
}

// Paralelna implementacija Euklidske norme vektora koristeći OpenMP
double vector_norm_parallel(double *x, int N)
{
    double suma = 0.0;

    // Paralelna petlja: više niti sabira djelove sume
    // reduction(+: suma) osigurava da se sve djelimične sume korektno objedine
#pragma omp parallel for reduction(+ : suma)
    for (int i = 0; i < N; i++)
        suma += x[i] * x[i];

    // Vraća kvadratni korijen ukupne sume
    return sqrt(suma);
}

// Sekvencijalna normalizacija vektora
void normalize_sequential(double *x, int N)
{
    // Izračunavanje Euklidske norme vektora
    double norm = vector_norm_sequential(x, N);

    // Svaki element vektora se dijeli sa normom
    // Time se vektor skalira tako da njegova dužina postane 1
    for (int i = 0; i < N; i++)
        x[i] /= norm;
}

// Paralelna normalizacija vektora koristeći OpenMP
void normalize_parallel(double *x, int N)
{
    // Prvo izračunaj normu vektora (paralelno)
    double norm = vector_norm_parallel(x, N);

// Svaki element vektora se paralelno dijeli sa normom
#pragma omp parallel for
    for (int i = 0; i < N; i++)
        x[i] /= norm;
}

// Sekvencijalna implementacija skalarnog proizvoda dva vektora
double dot_product_sequential(double *a, double *b, int N)
{
    double rezultat = 0.0;

    // Sabira proizvod odgovarajućih elemenata vektora a i b
    for (int i = 0; i < N; i++)
        rezultat += a[i] * b[i];

    // Vraća ukupni zbir — skalarni proizvod
    return rezultat;
}

// Paralelna implementacija skalarnog proizvoda pomoću OpenMP-a
double dot_product_parallel(double *a, double *b, int N)
{
    double rezultat = 0.0;

// Svaka nit računa dio proizvoda i koristi reduction da se sve sume objedine
#pragma omp parallel for reduction(+ : rezultat)
    for (int i = 0; i < N; i++)
        rezultat += a[i] * b[i];

    return rezultat;
}

// Sekvencijalna implementacija Power Iteration metode
void power_iteration_sequential(double **A, double *eigenvector, double *eigenvalue, int N)
{
    // Alocira memoriju za trenutni i sledeći vektor (iteracije)
    double *b_k = malloc(sizeof(double) * N);
    double *b_k1 = malloc(sizeof(double) * N);

    // Inicijalna pretpostavka: vektor sa svim jedinicama
    for (int i = 0; i < N; i++)
        b_k[i] = 1.0;

    // Glavna iteracija algoritma
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        // Množenje matrice A i vektora b_k: b_k1 = A * b_k
        matrix_vector_multiplication_sequential(A, b_k, b_k1, N);

        // Normalizacija rezultujućeg vektora
        normalize_sequential(b_k1, N);

        // Računanje razlike između stare i nove iteracije (za provjeru konvergencije)
        double diff = 0.0;
        for (int i = 0; i < N; i++)
        {
            diff += fabs(b_k1[i] - b_k[i]); // suma apsolutnih razlika
            b_k[i] = b_k1[i];               // ažurira b_k za sledeću iteraciju
        }

        // Ako je razlika manja od tolerancije, prekidamo iteraciju
        if (diff < TOL)
            break;
    }

    // Kopira konvergentni vektor kao svojstveni vektor
    for (int i = 0; i < N; i++)
        eigenvector[i] = b_k[i];

    // Računa svojstvenu vrijednost koristeći skalarni proizvod
    double *Ax = malloc(sizeof(double) * N);
    matrix_vector_multiplication_sequential(A, eigenvector, Ax, N);
    *eigenvalue = dot_product_sequential(Ax, eigenvector, N);

    // Oslobađa dinamički alociranu memoriju
    free(b_k);
    free(b_k1);
    free(Ax);
}

// Paralelna implementacija Power Iteration metode koristeći OpenMP
void power_iteration_parallel(double **A, double *eigenvector, double *eigenvalue, int N)
{
    // Alociranje memorije za trenutni i sledeći iterativni vektor
    double *b_k = malloc(sizeof(double) * N);
    double *b_k1 = malloc(sizeof(double) * N);

    // Inicijalizacija vektora sa jedinicama
    for (int i = 0; i < N; i++)
        b_k[i] = 1.0;

    // Glavna iteracija metode
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        // Paralelno množenje matrice i vektora: b_k1 = A * b_k
        matrix_vector_multiplication_parallel(A, b_k, b_k1, N);

        // Paralelna normalizacija rezultujućeg vektora
        normalize_parallel(b_k1, N);

        // Računanje ukupne razlike između prethodnog i trenutnog vektora
        double diff = 0.0;

// Svaka nit računa dio razlike, koristeći 'reduction' da sabere rezultate
#pragma omp parallel for reduction(+ : diff)
        for (int i = 0; i < N; i++)
        {
            diff += fabs(b_k1[i] - b_k[i]); // suma apsolutnih razlika
            b_k[i] = b_k1[i];               // ažurira b_k za sledeću iteraciju
        }

        // Ako je razlika ispod zadate tolerancije, prekidamo iteraciju
        if (diff < TOL)
            break;
    }

    // Kopiranje rezultujućeg vektora kao svojstvenog vektora
    for (int i = 0; i < N; i++)
        eigenvector[i] = b_k[i];

    // Izračunavanje svojstvene vrijednosti 
    double *Ax = malloc(sizeof(double) * N);
    matrix_vector_multiplication_parallel(A, eigenvector, Ax, N);
    *eigenvalue = dot_product_parallel(Ax, eigenvector, N);

    // Oslobađanje memorije
    free(b_k);
    free(b_k1);
    free(Ax);
}

// Funkcija za kreiranje simetrične matrice dimenzija NxN
double **create_matrix(int N)
{
    double **A = malloc(sizeof(double *) * N); // Alokacija redova matrice

    for (int i = 0; i < N; i++)
    {
        A[i] = malloc(sizeof(double) * N); // Alokacija kolona za red i

        // Popunjavanje donje trougaone matrice i kopiranje vrijednosti u gornju (simetrična matrica)
        for (int j = 0; j <= i; j++)
        {
            A[i][j] = A[j][i] = rand() / (double)RAND_MAX;
        }
    }

    return A;
}

// Oslobađanje memorije matrice
void free_matrix(double **A, int N)
{
    for (int i = 0; i < N; i++)
        free(A[i]); // Oslobađanje svakog reda

    free(A); // Oslobađanje pokazivača na redove
}

// Pokretanje benchmark testa za određenu dimenziju matrice
void run_benchmark(int N)
{
    // Kreira matricu i alocira memoriju za svojstveni vektor i vrijednost
    double **A = create_matrix(N);
    double *eigvec = malloc(sizeof(double) * N);
    double eigval;

    double time_seq = 0.0, time_omp = 0.0;

    // Višestruko pokretanje testa za pouzdaniji prosjek
    for (int t = 0; t < NUM_TRIALS; t++)
    {
        // Mjerenje vremena sekvencijalne verzije
        double start = omp_get_wtime();
        power_iteration_sequential(A, eigvec, &eigval, N);
        time_seq += omp_get_wtime() - start;

        // Mjerenje vremena paralelne verzije
        start = omp_get_wtime();
        power_iteration_parallel(A, eigvec, &eigval, N);
        time_omp += omp_get_wtime() - start;
    }

    // Prikaz prosječnih vremena i akceleracije
    printf("Matrix size: %4d | Seq avg: %.4fs | OMP avg: %.4fs | Speedup: %.2fx\n",
           N, time_seq / NUM_TRIALS, time_omp / NUM_TRIALS, time_seq / time_omp);

    // Oslobađanje memorije
    free(eigvec);
    free_matrix(A, N);
}

// Glavna funkcija koja pokreće benchmarke za više dimenzija matrica
int main()
{
    srand(42); // Inicijalizacija generatora slučajnih brojeva za konzistentne rezultate

    printf("Benchmarking Power Iteration (avg of %d trials each)\n", NUM_TRIALS);
    printf("--------------------------------------------------------\n");

    // Niz dimenzija matrica koje će se testirati
    int sizes[] = {100, 200, 400, 800, 1200, 1600, 2000, 4000, 10000};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    // Pokretanje benchmarka za svaku dimenziju
    for (int i = 0; i < num_sizes; i++)
    {
        run_benchmark(sizes[i]);
    }

    return 0;
}

