#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define SINGULAR_matriz_ERROR 1
#define MAX_IT_ERROR 2

/**********************************************************************************************************************************
                                                    INICIALIZAÇÃO DA MATRIZ
**********************************************************************************************************************************/
typedef struct{
    double **dados;
    const unsigned int m; // linhas
    const unsigned int n; // colunas
} matriz;

/**********************************************************************************************************************************
                                                    FUNÇÃO PARA INICIAR OS DADOS
**********************************************************************************************************************************/
double **ini_dados(unsigned int m, unsigned int n) {
    double **dados = (double **)malloc(sizeof(double *) * m);
    for (int i = 0; i < m; i++) {
        dados[i] = (double *)malloc(sizeof(double) * n);
    }
    return dados;
}

/**********************************************************************************************************************************
                                                    FUNÇÃO CRIA MATRIZ VAZIA
**********************************************************************************************************************************/
matriz *matriz_vazia(unsigned int m, unsigned int n) {
    matriz *M = (matriz *)malloc(sizeof(matriz));
    matriz s = (matriz){NULL, m, n};
    memcpy(M, &s, sizeof(matriz));
    M->dados = ini_dados(m, n);
    return M;
}

/**********************************************************************************************************************************
                                                    FUNÇÃO LIBERA MATRIZ DA MEMÓRIA
**********************************************************************************************************************************/
void matriz_libera(matriz *M) {
    for (int i = 0; i < M->m; i++) {
        free(M->dados[i]);
    }
    free(M->dados); free(M);
}

/**********************************************************************************************************************************
                                                    FUNÇÃO PREENCHE UM VALOR
**********************************************************************************************************************************/
matriz *matriz_preencher(unsigned int m, unsigned int n, int valor) {
    matriz *M = NULL;
    if (M = matriz_vazia(m, n))
        for (int i = 0; i < M->m; i++) {
            for (int j = 0; j < M->n; j++) {
                M->dados[i][j] = valor;
            }
        }   
    return M;
}

/**********************************************************************************************************************************
                                                    FUNÇÃO GAUSS
**********************************************************************************************************************************/
int gauss_seidel(const matriz *A, const matriz *x, const matriz *b, double e, int itmax){
    int n = A->m;
    double novo_valor;

    for (int it = 0; it < itmax; it++) {
        double max_delta = 0;
        for (int i = 0; i < n; i++) {
            novo_valor = b->dados[i][0];
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                novo_valor -= A->dados[i][j] * x->dados[j][0];
            }
            if (fabs(A->dados[i][i]) < DBL_EPSILON) {
                return SINGULAR_matriz_ERROR;
            }
            novo_valor /= A->dados[i][i];
            max_delta = MAX(max_delta, fabs(x->dados[i][0] - novo_valor));
            x->dados[i][0] = novo_valor;
        }
        if (max_delta < e) {
            return EXIT_SUCCESS;
        }
    }
    return MAX_IT_ERROR;
}

/**********************************************************************************************************************************
                                                        FUNÇÃO SUBTRAIR MATRIZES
**********************************************************************************************************************************/
matriz *matriz_sub(const matriz *A, const matriz *B) {
    // Verifica se as matrizes são da mesma ordem
    if ((A->m != B->m) || (A->n != B->n))
        return NULL;
    // Cria a matriz C e efetua a subtração
    matriz *C = matriz_vazia(A->m, A->n);
    for (int i = 0; i < C->m; i++) {
        for (int j = 0; j < C->n; j++) {
            C->dados[i][j] = A->dados[i][j];
        }
    }
    for (int i = 0; i < C->m; i++) {
        for (int j = 0; j < C->n; j++) {
            C->dados[i][j] -= B->dados[i][j];
        }
    }
    return C;
}

/**********************************************************************************************************************************
                                                        FUNÇÃO MULTIPLICA MATRIZES
**********************************************************************************************************************************/
matriz *matriz_mult(const matriz *A, const matriz *B) {
    // Verifica se é possivel efetuar a multiplicação (pela ordem das matrizes)
    if (A->n != B->m)
        return NULL;
    // Cria a matriz C e efetua a multiplicação
    matriz *C = matriz_preencher(A->m, B->n, 0);
    for (int i = 0; i < C->m; i++) {
        for (int j = 0; j < C->n; j++) {
            for (int k = 0; k < A->n; k++) {
                C->dados[i][j] += A->dados[i][k] * B->dados[k][j];
            }
        }
    }
    return C;
}


/**********************************************************************************************************************************
                                                        FUNÇÃO CALCULA NORMA
**********************************************************************************************************************************/
double norm_calc(const matriz *A) {
    double max_sum = 0;
    for (int i = 0; i < A->m; i++) {
        double sum = 0;
        for (int j = 0; j < A->n; j++) {
            sum += fabs(A->dados[i][j]);
        }
        if (sum > max_sum) {
            max_sum = sum;
        }
    }
    return max_sum;
}


/**********************************************************************************************************************************
                                                        FUNÇÃO DIFERENÇA DE NORMA
**********************************************************************************************************************************/
double norm_dif(const matriz *A, const matriz *B){
    matriz *C = matriz_sub(A, B);
    double norm = norm_calc(C);
    matriz_libera(C);
    return norm;
}

/**********************************************************************************************************************************
                                                        FUNÇÃO CRIA MATRIZ COM VALORES DO ENUNCIADO
**********************************************************************************************************************************/
matriz *cria_matriz(int n) {
    matriz *A = matriz_preencher(n, n, 0);
    for (int i = 0; i < n; i++) {
        A->dados[i][i] = 3;
        if (i + 1 < n) {
            A->dados[i][i + 1] = -0.75;
            A->dados[i + 1][i] = -0.75;
        }
        if (i + 3 < n) {
            A->dados[i][i + 3] = -0.75;
            A->dados[i + 3][i] = -0.75;
        }
    }
    return A;
}

/**********************************************************************************************************************************
                                                        FUNÇÃO PRIMEIRO CASO
***********************************************************************************************************************************/
void caso1(int n, double e, int itmax) {
    matriz *A = NULL; matriz *x = NULL; matriz *b = NULL;
    A = cria_matriz(n);
    x = matriz_preencher(n, 1, 0);
    b = matriz_preencher(A->m, 1, 0);
    for (int i = 0; i < b->m; i++) {
        for (int j = 0; j < A->n; j++) {
            b->dados[i][0] += A->dados[i][j];
        }
    }

    int gs_error;
    if (gs_error = gauss_seidel(A, x, b, e, itmax)) {
        if (gs_error == 1)
            fprintf(stderr, "Função Gauss seidel: Divisão por 0.\n");
        if(gs_error == 2)
            fprintf(stderr, "Função Gauss seidel: Máximo de interações.\n");
    }
    matriz *Ax = matriz_mult(A, x);

    matriz *ones = matriz_preencher(n, 1, 1);
    printf("||x - 1||∞ = %e\n", norm_dif(x, ones));
    printf("||b - Ax||∞ = %e\n", norm_dif(b, Ax));

    matriz_libera(A); matriz_libera(x); matriz_libera(b); matriz_libera(Ax); matriz_libera(ones);
}
/**********************************************************************************************************************************
                                                        FUNÇÃO SEGUNDO CASO
***********************************************************************************************************************************/
void caso2(int n, double e, int itmax) {
    matriz *A = NULL; matriz *x = NULL; matriz *b = NULL;
    A = cria_matriz(n);
    x = matriz_preencher(n, 1, 0);
    b = matriz_preencher(n, 1, 1);
    for (int i = 0; i < n; i++) {
        b->dados[i][0] = 2 / (i + 1);
    }

    int gs_error;
    if (gs_error = gauss_seidel(A, x, b, e, itmax)) {
        if (gs_error == 1)
            fprintf(stderr, "Função Gauss seidel: Divisão por 0.\n");
        if (gs_error == 2)
            fprintf(stderr, "Função Gauss seidel: Máximo de interações.\n");
    }

    matriz *ones = matriz_preencher(n, 1, 1);

    printf("||x - 1||∞ = %e\n", norm_dif(x, ones));

    matriz_libera(A); matriz_libera(x); matriz_libera(b);
}

/**********************************************************************************************************************************
                                                        FUNÇÃO PRINCIPAL
***********************************************************************************************************************************/
int main(){
    double e = 1e-10;
    int itmax = 100000;

    // 1º caso para n = 100:
    printf("1º caso de teste para n = 100:\n");
    caso1(100, e, itmax);

    // 1º caso para n = 200:
    printf("1º caso de teste para n = 200:\n");
    caso1(200, e, itmax);

    // 2º caso para n = 1000:
    printf("2º caso de teste para n = 1000:\n");
    caso2(1000, e, itmax);

    return 0;
}