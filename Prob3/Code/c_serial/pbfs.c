#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>
#include <mpi.h>

#define IF if((wrank)==0)
#define Fprintf fprintf

typedef struct graphstruct { // A graph in compressed-adjacency-list (CSR) form
	int nv;            // number of vertices
	int64_t ne;            // number of edges
	int *nbr;          // array of neighbors of all vertices
	int *firstnbr;     // index in nbr[] of first neighbor of each vtx
} graph;

int read_graph_from_file(graph **G, int *startvtx, const char* filename);
void print_CSR_graph (const graph *G);
void bfs    (const int s, const graph *G, int *nlevelsp, int **levelsizep);

static inline double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1.e-6;
}

int read_graph_from_file(graph **G, int *startvtx, const char* filename) {

	FILE *fp = fopen (filename, "r");
	if ( fp == NULL ) {
		Fprintf(stderr, "error to open file %s: %s\n", filename, strerror(errno));
		return -1;
	}
	*G = (graph *) malloc(sizeof(graph));
	fread(startvtx, 1, sizeof(int), fp);
	fread(&((*G)->ne), 1, sizeof(int64_t), fp);
	fread(&((*G)->nv), 1, sizeof(int), fp);
	(*G)->nbr = (int *) malloc((*G)->ne*sizeof(int));
	(*G)->firstnbr = (int *) malloc(((*G)->nv+1)*sizeof(int));
	fread((*G)->nbr, (*G)->ne, sizeof(int), fp);
	fread((*G)->firstnbr, (*G)->nv+1, sizeof(int), fp);

	fclose(fp);
	return 0;
}

void print_CSR_graph (const graph *G) {
	int vlimit = 20;
	int elimit = 50;
	int e,v;
	Fprintf(stdout, "\nGraph has %d vertices and %"PRId64" edges.\n",G->nv,G->ne);
	Fprintf(stdout, "firstnbr =");
	if (G->nv < vlimit) vlimit = G->nv;
	for (v = 0; v <= vlimit; v++) printf(" %d",G->firstnbr[v]);
	if (G->nv > vlimit) Fprintf(stdout, " ...");
	Fprintf(stdout, "\n");
	Fprintf(stdout, "nbr =");
	if (G->ne < elimit) elimit = G->ne;
	for (e = 0; e < elimit; e++) printf(" %d",G->nbr[e]);
	if (G->ne > elimit) Fprintf(stdout, " ...");
	Fprintf(stdout, "\n\n");
}

void bfs   (const int s,
            const graph *G,
            int *nlevelsp,
            int **levelsizep, int wrank, int wsize) {
	int *levelsize;
	int thislevel; 
	int *queue, back, front; 
	int i, v, w, e; 
	int *level = (int *) malloc(G->nv*sizeof(int)); 
	int *llevel = (int *) malloc(G->nv*sizeof(int)); 
	levelsize = *levelsizep = (int *) malloc(G->nv*sizeof(int)); 
	queue = (int *) malloc(G->nv*sizeof(int));
    IF gqueue = (int *) malloc(G->nv*sizeof(int)*wsize);
	IF bmap = (int *) calloc(G->nv, sizeof(int));
	lqueue = (int *) malloc(G->nv*sizeof(int));
    int *num, *sdex, lback, lfront; 
    num = (int *)malloc(wsize);
    sdex = (int *)malloc(wsize);


	// initially, queue is empty, all levels are -1
	back = 0;   // position next element will be added to queue
	front = 0;  // position next element will be removed from queue
    lback = 0;
    lfront = 0;
	for (v = 0; v < G->nv; v++) level[v] = -1;

	// assign the starting vertex level 0 and put it on the queue to explore
    //
	thislevel = 0;
	level[s] = 0;
	levelsize[0] = 1;
	queue[back++] = s;
    int mynum; 
	// loop over levels, then over vertices at this level, then over neighbors
	while (levelsize[thislevel] > 0) {
		levelsize[thislevel+1] = 0;

        // ** Assign the Job to the cores
        v = 0; 
        for(i=0; i<wsize; i++)
        {
            sdex[i] = v;
            num[i] = levelsize[thislevel]/wsize;
            if (i < levelsize[thislevel]%wsize)
                num[i] += 1; 
            v += num[i];
        }

        MPI_Scatterv(queue+front, num, sdex, MPI_INT,\
                lqueue, num[wrank], MPI_INT, 0, MPI_COMM_WORLD);
        lfront = 0;
        lback = num[wrank];

        // ** local queue - 
        // ** level - Update information of other
        // ** levelsize - sum queue
        for (i = 0; i < num[wrank]; i++) {
            v = lqueue[lfront++];
			for (e = G->firstnbr[v]; e < G->firstnbr[v+1]; e++) {
				w = G->nbr[e];          // w is the current neighbor of v
				if (llevel[w] == -1) {   // w has not already been reached
					llevel[w] = thislevel+1;
					levelsize[thislevel+1]++;
					lqueue[back++] = w;    // put w on queue to explore
				}
			}
		}

        mynum = num[wrank];
        //queue, level size
        MPI_Gather(levelsize+thislevel+1, 1, MPI_INT,\
                num, 1, MPI_INT, 0, MPI_COMM_WORLD);
        v = 0;
        for(i=0; i<wsize; i++)
        {
            sdex[i] = v;
            v += num[i];
        }

        MPI_Gatherv(lqueue+mynum, levelsize[thislevel+1], MPI_INT, \
                gqueue, num, sdex, MPI_INT, MPI_COMM_WORLD);
        MPI_Reduce(levelsize+thislevel+1, levelsize+thislevel+1, 1, MPI_INT,\
            MPI_SUM, 0, MPI_COMM_WORLD);

        
        //Remove Overlap
        //
        IF
        {
            for(i=0; i<v ;i++)
            {
                if(bmap[gqueue[i]]==0)
                    queue[back++]=gqueue[i];
                else
                    levelsize[thislevel+1]--;
            }
        }
        //level
        MPI_Allreduce(level, llevel, G->nv, MPI_INT, MP_MAX, MPI_COMM_WORLD);
        memcpy(level, llevel, G->nv * sizeof(int));
        // ** Fuse the data in the levels
		thislevel = thislevel+1;
	}
	*nlevelsp = thislevel;
	free(queue);
	free(level);

    // Added
    IF free(gqueue);
    IF free(bmap);
    free(num);
    free(sdex);
    free(lqueue);
}

int main (int argc, char* argv[]) {
	graph *G;
	int *levelsize;
	int nlevels;
	int startvtx;
	int i, reached;
	double starttime, elapsedtime;
	const char* filename;
    int wrank, wsize;
	if (argc == 2) {
		filename = argv[1];
	} else {
		Fprintf(stderr, "usage: bfs <filename>\n");
		exit(1);
	}
//	starttime = get_time();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

	starttime = MPI_Wtime();
	if (read_graph_from_file(&G, &startvtx, filename)) return -1;
	Fprintf(stdout, "Elapsed Time to read and construct graph: %f\n", get_time() - starttime);
	print_CSR_graph (G);

	Fprintf(stdout, "Starting vertex for BFS is %d.\n\n",startvtx);
	starttime = get_time();
    // Begin the code
	bfs (startvtx, G, &nlevels, &levelsize, wrank, wsize);
    elapsedtime = MPI_Wtime() - starttime;
	//elapsedtime = get_time() -starttime;

	reached = 0;
	for (i = 0; i < nlevels; i++) reached += levelsize[i];
	Fprintf(stdout, "Breadth-first search from vertex %d reached %d levels and %d vertices.\n",
		startvtx, nlevels, reached);
	for (i = 0; i < nlevels; i++) printf("level %d vertices: %d\n", i, levelsize[i]);
	Fprintf(stdout, "\n");
	Fprintf(stdout, "GTEPs: %f\n", (reached/elapsedtime)/1000000);
	Fprintf(stdout, "Elapsed Time: %f\n", elapsedtime);

	free(levelsize);
	free(G);
    MPI_Finalize();
}
