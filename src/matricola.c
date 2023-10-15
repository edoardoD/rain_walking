#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

#define GRAPH_INDEX(x, y, col) ((x) * (col) + (y))
/*definizione delle strutture necessarie*/
typedef enum
{
    false,
    true
} bool;

typedef struct Node
{
    int id;
    bool is_wet;
} Node;

/* struttura arco */
typedef struct Edge
{
    Node src;      /* nodo sorgente */
    Node dst;      /* nodo destinazione */
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

/***
 * @author Edoardo Desiderio
 * dato il file input la funzione legge la matrice descritta
 * e fornisce un grafo non or***/
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

/***
 * @author Edoardo Desiderio
 * caricamento della matrice, dato lo specifico caso
 * sarebbe comodo ricondursi a costruire un grafo che permetta
 * di risolvere il problema nel miglior compromesso possibile
 * fra utilizzo di memoria e costo computazionale
 * @param argc
 * @param argv
 ***/
int main(int argc, char const *argv[])
{
    int n, m;
    unsigned int i = 0, j = 0, skyscapers = 0;
    int *v,**matrix;
    char c;
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
    fscanf(filein, "%d %d", &n, &m);

    /*inizializzo la matrice con calloc in maniera da pulire le memeorie*/
    matrix = (int **)calloc(n, sizeof(int *));
    c = fgetc(filein);
    while (c != EOF)
    {
        if (isdigit(c))
        {

            matrix[i][j] = atoi(&c);
            if (matrix[i][j] > 0)
            {
                /*mi serve per poter allocare solo i nodi necessari*/
                skyscapers++;
            }
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
        }
        c = fgetc(filein);
    }
    fclose(filein);
    
    printf("skyscapers: %d\n", skyscapers);
    g = graph_create(n * m - skyscapers);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%d ", GRAPH_INDEX(i, j, m));
            v = matrix[i];
        }
        printf("\n");
    }

    return EXIT_SUCCESS;
}
