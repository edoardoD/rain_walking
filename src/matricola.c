#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>

#define GRAPH_INDEX(x, y, col) ((x) * (col) + (y))
#define IS_DRY -1
/*definizione delle strutture necessarie*/
typedef enum
{
    false,
    true
} bool;

/* struttura arco */
typedef struct Edge
{
    int src;       /* nodo sorgente */
    int dst;       /* nodo destinazione */
    double weight; /* peso dell'arco*/
    struct Edge *next;
} Edge;

/* struttura grafo */
typedef struct
{
    int n;        /* numero di nodi               */
    int m;        /* numero di archi              */
    Edge **edges; /* array di liste di adiacenza  */
    int *in_deg;  /* grado entrante dei nodi      */
    int *out_deg; /* grado uscente dei nodi       */
} Graph;
/*fine definizione delle strutture*/

/*funzioni per il grafo*/
int graph_n_nodes(const Graph *g)
{
    assert(g != NULL);

    return g->n;
}

Graph *graph_create(int n)
{
    int i;
    Graph *g = (Graph *)malloc(sizeof(*g));
    assert(g != NULL);
    assert(n > 0);

    g->n = n;
    g->m = 0;
    g->edges = (Edge **)malloc(n * sizeof(Edge *));
    assert(g->edges != NULL);
    g->in_deg = (int *)malloc(n * sizeof(*(g->in_deg)));
    assert(g->in_deg != NULL);
    g->out_deg = (int *)malloc(n * sizeof(*(g->out_deg)));
    assert(g->out_deg != NULL);
    for (i = 0; i < n; i++)
    {
        g->edges[i] = NULL;
        g->in_deg[i] = g->out_deg[i] = 0;
    }
    return g;
}
Edge *graph_adj(const Graph *g, int v)
{
    assert(g != NULL);
    assert((v >= 0) && (v < graph_n_nodes(g)));

    return g->edges[v];
}
/***
 * Edoardo Desiderio
 * controlla se esite già un nodo che abbia sorgente e destinatario passato,
 * @returns 1 se esiste, 0 altrimenti
 */
static int exsit_edge(Graph *g, int src, int dst)
{
    Edge *e;
    for (e = graph_adj(g, src); e != NULL; e = e->next)
    {
        if (e->src == src && e->dst == dst)
        {
            return 1;
        }
    }
    return 0;
}

void graph_add_edge(Graph *g, int src, int dst, double weight)
{
    Edge *new_edge = (Edge *)malloc(sizeof(Edge));

    assert((src >= 0) && (src < graph_n_nodes(g)));
    assert((dst >= 0) && (dst < graph_n_nodes(g)));
    assert(new_edge != NULL);

    if (!exsit_edge(g, src, dst))
    {
        new_edge->src = src;
        new_edge->dst = dst;
        new_edge->weight = weight;
        new_edge->next = graph_adj(g, src);
        g->edges[src] = new_edge; /*aggiorno la lista di adiacenza*/
        g->out_deg[src]++;
        g->in_deg[dst]++;

        if (!exsit_edge(g, dst, src))
        {
            graph_add_edge(g, dst, src, weight);
        }
    }
    else
    {
        fprintf(stderr, "ERRORE: arco duplicato (%d,%d)\n", src, dst);
    }
}

int graph_out_degree(const Graph *g, int v)
{
    assert(g != NULL);
    assert((v >= 0) && (v < graph_n_nodes(g)));
    return g->out_deg[v];
}

void graph_print(const Graph *g)
{
    int i;

    assert(g != NULL);

    for (i = 0; i < g->n; i++)
    {
        const Edge *e;
        int out_deg = 0; /* ne approfittiamo per controllare la
                            correttezza dei gradi uscenti */
        printf("[%2d] -> ", i);
        for (e = graph_adj(g, i); e != NULL; e = e->next)
        {
            printf("(%d, %d, %f) -> ", e->src, e->dst, e->weight);
            out_deg++;
        }
        assert(out_deg == graph_out_degree(g, i));
        printf("NULL\n");
    }
}
/***
 * Edoardo Desiderio
 * legge le righe e le colonne della matrice
 * @return un puntatore di 2 elementi 
 * ***/
int *read_dimension(FILE *filein)
{
    int *dim = (int *)malloc(2 * sizeof(int));
    fscanf(filein, "%d %d", &dim[0], &dim[1]);
    return dim;
}

/***
 * Edoardo Desiderio
 * funzione che legge la matrice dal file e la carica in memoria.
 * Durante il ciclo while mentre carica in memoria la matrice
 * distingue se la cella risulta bagnata oppure no, nel caso
 * assegna il valore IS_DRY
 * @param filein
 * @param n righe
 * @param m colonne
 * @returns la matrice, e nella posizione n+1 il numero di palazzi
 ***/
int **load_matrix(FILE *filein, int n, int m)
{
    int i = 0, j = 0, skyscapers=0;
    int  tmp;
    int delta = INT_MIN;
    int **matrix;
    char c;
    /*inizializzo la matrice con calloc in maniera da pulire le memeorie*/
    matrix = (int **)calloc(n+1, sizeof(int *));
    c = fgetc(filein);
    while (c != EOF)
    {
        if (isdigit(c))
        {
            tmp = atoi(&c);
            if (tmp > 0)
            {
                /*mi serve per poter allocare solo i nodi necessari
                al grafo*/
                skyscapers++;
                /*bug sul delta, un palazzo da 9 piani copre anche le celle 
                che non coprirebbero palazzi piccoli successivi*/
                if (tmp+j > delta)
                {
                    delta = tmp + j <= m ? tmp + j : m;
                }
            }
            else
            {
                if (j <= delta)
                {
                    tmp = IS_DRY;
                }
            }
            matrix[i][j] = tmp;
            j++;
            if (j == m)
            {
                i++;
                j = 0;
            }
        }
        else
        {
            matrix[i] = (int *)calloc(m, sizeof(int));
            delta = INT_MIN;
        }
        c = fgetc(filein);
    }
   
    fclose(filein);
    matrix[n+1] = (int*)malloc(sizeof(int));
    matrix[n+1][0] = skyscapers;
    /*printf("skyscapers: %d\n", matrix[n+1][0]);*/
    return matrix;
}
/***
 * Edoardo Desiderio.
 * funzione che riempie il grafo con i nodi e gli archi,
 * un arco entrante in una cella IS_DRY avrà peso 1.0, altrimenti
 * verrà calcolato il peso con un numero double randomico sufficientemente
 * in modo da far pesare molto il fatto che si passi in una cella bagnata
 * @param g grafo
 * @param matrix matrice
 * @param n righe
 * @param m colonne
 *  
 * ***/
void fill_graph(Graph* g,int** matrix, int n,int m){
    
    int i,j;
    double wet_weigth;
    srand(time(NULL));
    wet_weigth = (double)rand() / (double)RAND_MAX;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (matrix[i][j] <= 0)
            {
                if (i + 1 < n && matrix[i + 1][j] <= 0)
                {
                    graph_add_edge(g, GRAPH_INDEX(i, j, m), GRAPH_INDEX(i + 1, j, m), matrix[i + 1][j] == IS_DRY ? 1.0 : wet_weigth);
                }
                if (j + 1 < m && matrix[i][j + 1] <= 0)
                {
                    graph_add_edge(g, GRAPH_INDEX(i, j, m), GRAPH_INDEX(i, j + 1, m), matrix[i][j + 1] == IS_DRY ? 1.0 : wet_weigth);
                }
            }
        }
    }
}

/***
 * @author Edoardo Desiderio
 * entry point del programma, avvia le varie funzioni
 * necessarie alla risoluzione del problema. 
 * Il programma accetta un input corrispondente alla 
 * matrice descritta dal problema
 * @param argc
 * @param argv
 ***/
int main(int argc, char const *argv[])
{
    int *m_n;
    int m, n;
    int skyscapers;
    int **matrix;
    FILE *filein = stdin;
    Graph *g;
    if (argc != 2)
    {
        fprintf(stderr, "Invocare il programma con: %s matrice\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "-") != 0)
    {
        filein = fopen(argv[1], "r");
        if (filein == NULL)
        {
            fprintf(stderr, "Can not open %s\n", argv[1]);
            return EXIT_FAILURE;
        }
    }
    printf("matrice caricata da %s\n", argv[1]);

    /*leggo la prima riga del file per le dimensioni*/
    m_n = read_dimension(filein);
    m = m_n[0];
    n = m_n[1];
    free(m_n);
    matrix = load_matrix(filein, n, m);
    skyscapers = matrix[n+1][0];
    free(matrix[n+1]);

    
    g = graph_create(n * m - skyscapers);
    assert(g != NULL);
    fill_graph(g, matrix, n, m);
    graph_print(g);
    return EXIT_SUCCESS;
}
