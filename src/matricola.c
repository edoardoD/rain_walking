#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#include <limits.h>

#define GRAPH_INDEX(x, y, col) ((x) * (col) + (y))
#define IS_DRY -1


/*min heap */
const int NODE_UNDEF = -1;

typedef struct {
    int key;
    double prio;
} HeapElem;

typedef struct {
    HeapElem *heap;
    int *pos; /* TODO: non inizializzato! Usare solo dopo una implementazione funzionante di minheap_change_prio() in cui la chiave da modificare si trova con una normale ricerca lineare sull'array heap[] */
    int n; /* quante coppie (chiave, prio) sono effettivamente presenti nello heap */
    int size; /* massimo numero di coppie (chiave, prio) che possono essere contenuti nello heap */
} MinHeap;

MinHeap *Q;

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


void minheap_print(const MinHeap *h)
{
    int i, j, width = 1;

    assert(h != NULL);

    printf("\n** Contenuto dello heap:\n\n");
    printf("n=%d size=%d\n", h->n, h->size);
    printf("Contenuto dell'array heap[] (stampato a livelli:\n");
    i = 0;
    while (i < h->n) {
        j = 0;
        while (j<width && i < h->n) {
            printf("h[%2d]=(%2d, %6.2f) ", i, h->heap[i].key, h->heap[i].prio);
            i++;
            j++;
        }
        printf("\n");
        width *= 2;
    }
    printf("\nContenuto dell'array pos[]:\n");
    for (i=0; i<h->size; i++) {
        printf("pos[%d]=%d ", i, h->pos[i]);
    }
    printf("\n\n** Fine contenuto dello heap\n\n");
}

void minheap_clear( MinHeap *h )
{
    int i;
    assert(h != NULL);
    for (i=0; i<h->size; i++) {
        h->pos[i] = -1;
    }
    h->n = 0;
}

/* Costruisce un min-heap vuoto che può contenere al massimo
   `size` elementi */
MinHeap *minheap_create(int size)
{
    MinHeap *h = (MinHeap*)malloc(sizeof(*h)) ;
    assert(h != NULL);
    assert(size > 0);
    h->size = size;
    h->heap = (HeapElem*)malloc(size * sizeof(*(h->heap)));
    assert(h->heap != NULL);
    h->pos = (int*)malloc(size * sizeof(*(h->pos)));
    assert(h->pos != NULL);
    minheap_clear(h);
    return h;
}

void minheap_destroy( MinHeap *h )
{
    assert(h != NULL);

    h->n = h->size = 0;
    free(h->heap);
    free(h->pos);
    free(h);
}

/* Restituisce true (nonzero) se l'indice `i` è valido per lo heap
   `h`, false (zero) altrimenti. */
static int valid(const MinHeap *h, int i)
{
    assert(h != NULL);

    return ((i >= 0) && (i < h->n));
}

/* Scambia heap[i] con heap[j] */
static void swap(MinHeap *h, int i, int j)
{
    HeapElem tmp;

    assert(h != NULL);
    assert(valid(h, i));
    assert(valid(h, j));
    assert(h->pos[h->heap[i].key] == i);
    assert(h->pos[h->heap[j].key] == j);

    tmp = h->heap[i];
    h->heap[i] = h->heap[j];
    h->heap[j] = tmp;

    h->pos[h->heap[i].key] = i;
    h->pos[h->heap[j].key] = j;
}

/* Restituisce l'indice del padre del nodo i */
static int parent(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return (i+1)/2 - 1;
}

/* Restituisce l'indice del figlio sinistro del nodo `i`. Ritorna un
   indice non valido se `i` non ha figlio sinistro. */
static int lchild(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return 2*i + 1;
}

/* Restituisce l'indice del figlio destro del nodo `i`. Ritorna un
   indice non valido se `i` non ha figlio destro. */
static int rchild(const MinHeap *h, int i)
{
    assert(valid(h, i));

    return 2*i + 2;
}

/* Restituisce l'indice del figlio di `i` con priorità minima. Se `i`
   non ha figli, restituisce -1 */
static int min_child(const MinHeap *h, int i)
{
    int l, r, result = -1;

    assert(valid(h, i));

    l = lchild(h, i);
    r = rchild(h, i);
    if (valid(h, l)) {
        result = l;
        if (valid(h, r) && (h->heap[r].prio < h->heap[l].prio)) {
            result = r;
        }
    }
    return result;
}

/* Scambia l'elemento in posizione `i` con il padre fino a quando
   raggiunge la posizione corretta nello heap */
static void move_up(MinHeap *h, int i)
{
    int p;

    assert(valid(h, i));

    p = parent(h, i);
    while ( valid(h, p) && (h->heap[i].prio < h->heap[p].prio) ) {
        swap(h, i, p);
        i = p;
        p = parent(h, i);
    }
}

/* Scambia l'elemento in posizione `i` con il figlio avente priorità
   minima, fino a quando l'elemento raggiunge la posizione
   corretta. Questa funzione corrisponde a Min-Heapify() */
static void move_down(MinHeap *h, int i)
{
    int done = 0;

    assert(valid(h, i));

    /* L'operazione viene implementata iterativamente, sebbene sia
       possibile una implementazione ricorsiva probabilmente più
       leggibile. */
    do {
        const int dst = min_child(h, i);
        if (valid(h, dst) && (h->heap[dst].prio < h->heap[i].prio)) {
            swap(h, i, dst);
            i = dst;
        } else {
            done = 1;
        }
    } while (!done);
}

/* Restituisce true (nonzero) se lo heap è vuoto */
int minheap_is_empty(const MinHeap *h)
{
    assert(h != NULL);

    return (h->n == 0);
}

/* Restituisce true (nonzero) se lo heap è pieno, cioè è stata
   esaurita la capienza a disposizione */
int minheap_is_full(const MinHeap *h)
{
    assert(h != NULL);

    return (h->n == h->size);
}

/* Restituisce il numero di elementi presenti nello heap */
int minheap_get_n(const MinHeap *h)
{
    assert(h != NULL);

    return h->n;
}

/* Restituisce la chiave associata alla priorità minima */
int minheap_min(const MinHeap *h)
{
    assert( !minheap_is_empty(h) );

    return h->heap[0].key;
}

/* Come minheap_min(), ma restituisce la coppia (chiave, prio) */
HeapElem minheap_min2( const MinHeap *h)
{
    assert( !minheap_is_empty(h) );

    return h->heap[0];
}

/* Inserisce una nuova coppia (key, prio) nello heap. */
void minheap_insert(MinHeap *h, int key, double prio)
{
    int i;

    assert( !minheap_is_full(h) );
    assert((key >= 0) && (key < h->size));
    assert(h->pos[key] == -1);

    i = h->n++;
    h->pos[key] = i;
    h->heap[i].key = key;
    h->heap[i].prio = prio;
    move_up(h, i);
}

/* Rimuove la coppia (chiave, priorità) con priorità minima;
   restituisce la chiave associata alla priorità minima. */
int minheap_delete_min(MinHeap *h)
{
    int result;

    assert( !minheap_is_empty(h) );

    result = minheap_min(h);
    swap(h, 0, h->n-1);
    assert( h->heap[h->n - 1].key == result );
    h->pos[result] = -1;
    h->n--;
    if (!minheap_is_empty(h)) {
        move_down(h, 0);
    }
    return result;
}

/* Come minheap_delete_min(), ma restituisce la coppia (chiave, prio) */
HeapElem minheap_delete_min2(MinHeap *h)
{
    HeapElem result;

    assert( !minheap_is_empty(h) );

    result = minheap_min2(h);
    swap(h, 0, h->n-1);
    assert( h->heap[h->n - 1].key == result.key );
    h->pos[result.key] = -1;
    h->n--;
    if (!minheap_is_empty(h)) {
        move_down(h, 0);
    }
    return result;
}

/* Modifica la priorità associata alla chiave key. La nuova priorità
   può essere maggiore, minore o uguale alla precedente. */
void minheap_change_prio(MinHeap *h, int key, double newprio)
{
    int j;
    double oldprio;

    assert(h != NULL);
    assert(key >= 0 && key < h->size);
    
    j = h->pos[key];
    assert( valid(h, j) );
    oldprio = h->heap[j].prio;
    h->heap[j].prio = newprio;
    if (newprio > oldprio) {
        move_down(h, j);
    } else {
        move_up(h, j);
    }
}

/*--------------FINE MINHEAP---------------*/

/* INIZIO CODICE SUL GRAFO */

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
    int i = 0, j = 0;
    int  tmp;
    int delta = INT_MIN;
    int **matrix;
    int skyscapers=0;
    char c;
    /*inizializzo la matrice con calloc in maniera da pulire le memeorie*/
    matrix = (int **)calloc(n+1, sizeof(*matrix));
   
    c = fgetc(filein);
    while (c != EOF)
    {
        if (isdigit(c))
        {
            tmp = atoi(&c);
            if (tmp > 0)
            {
                skyscapers++;
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
            matrix[i] = (int *)calloc(m, sizeof(**matrix));
            delta = INT_MIN;
        }
        c = fgetc(filein);
    }
   
    matrix[n+1] = (int*)malloc(sizeof(**matrix));
    *matrix[n+1] = skyscapers;
    

    return matrix;
}
/*Genera double random sicuramente maggiore del valore affidato per archi entranti in 
celle asciutte */
static double generate_random_number()
{
    double random_number; 
    srand(time(NULL));
    random_number = ((double)rand() / (double)RAND_MAX) * 100.0;
    while (random_number <= 3.0) {
        random_number = ((double)rand() / (double)RAND_MAX) * 100.0;
    }
    return random_number;
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
    double wet_weigth = generate_random_number();
    
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

/*FINE CODICE SUL GRAFO*/

/*algoritmo di dykstra*/
void dijkstra( const Graph *g, int s, double *d, int *p, const Edge **sp)
{
    int i;
    int v,u;
    const Edge *e;

    /*inizializzazione contenuti*/
    for ( i = 0; i < g->n; i++)
    {
        sp[i]=NULL;
        d[i] = HUGE_VAL;
        p[i] = NODE_UNDEF;
    }
    d[s] = 0 ;
    p[s] = NODE_UNDEF;
    /*inizializzazione coda di priorità*/
    Q = minheap_create(g->n);

    
    for ( i = 0; i < g->n; i++)
    {
        minheap_insert(Q, i, d[i]);
    }

    /*ciclo principale*/
    while (!minheap_is_empty(Q))
    {
        u = minheap_delete_min(Q);
        e = graph_adj(g, u);
        while (e != NULL)
        {
            v = e->dst;
            if (d[v] > d[u] + e->weight)
            {
                d[v] = d[u] + e->weight;
                p[v] = u;
                sp[v] = e;
                minheap_change_prio(Q,v,d[v]);
            }
            e = e->next;
        }
    }
    /*beauty print array Q->pos*/
}
void print_path(const int *p, int src, int dst)
{
    if (src == dst)
    {
        printf("%d",src);
    }
    else if (p[dst] == NODE_UNDEF)
    {
        printf("Non raggiungibile");
    }
    else
    {
        print_path(p,src,p[dst]);
        printf("->%d",dst);
    }
    
}
void print_dist( const Graph *g, int src, int dst, const int *p, const double *d )
{
    const int n = graph_n_nodes(g);
    int v, from_idx = 0, to_idx = n-1;

    assert(p != NULL);
    assert(d != NULL);

    if (dst >= 0) {
        from_idx = to_idx = dst;
    }

    printf(" src  dst            d path\n");
    printf("---- ---- ------------ -------------------------\n");
    for (v=from_idx; v<=to_idx; v++) {
        printf("%4d %4d %12.4f ", src, v, d[v]);
        print_path(p, src, v);
        printf("\n");
    }

}

void matrix_destroy(int **matrix, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void print_matrix(int **matrix, int m, int n)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%d\t", GRAPH_INDEX(i,j,n));
        }
        printf("\n");
    }
}
/**/
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
    int m, n,skyscapers;
    int **matrix;
    FILE *filein = stdin;
    const Edge **sp; /* sp[v] è il puntatore all'arco che collega v
                        con il suo predecessore nell'albero dei
                        cammini minimi */
    double *d;       /* d[v] è la distanza minima dalla sorgente al
                        nodo v */
    int *p;         /* p[v] è il predecessore di v nel cammino
                        minimo dalla sorgente */
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
    skyscapers = *matrix[n+1];
    free(matrix[n+1]);
    matrix[n+1] = NULL;
    print_matrix(matrix, n, m);
    printf("skyscapers: %d\n", skyscapers);
    g = graph_create(n * m);
    assert(g != NULL);
    fill_graph(g, matrix, n, m);
    
    matrix_destroy(matrix,n);
    
    printf("Grafo creato\n");

    d = (double*)malloc(n * sizeof(*d)); assert(d != NULL);
    p = (int*)malloc(n * sizeof(*p)); assert(p != NULL);
    sp = (const Edge**)malloc(n * sizeof(*sp)); assert(sp != NULL);

    dijkstra(g, 0, d, p, sp);
    print_path(p,0,GRAPH_INDEX(n-1,m-1,m));
    
    return EXIT_SUCCESS;
}
