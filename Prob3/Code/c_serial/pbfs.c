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
#define Fprintf if((wrank)==0)fprintf
#define Dprintf if(0)printf

typedef struct graphstruct { // A graph in compressed-adjacency-list (CSR) form
	int nv;            // number of vertices
	int64_t ne;            // number of edges
	int *nbr;          // array of neighbors of all vertices
	int *firstnbr;     // index in nbr[] of first neighbor of each vtx
} graph;

int read_graph_from_file(graph **G, int *startvtx, const char* filename, int wrank);
void print_CSR_graph (const graph *G, int wrank);
void bfs    (const int s, const graph *G, int *nlevelsp, int **levelsizep, int wrank, int wsize);

static inline double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1.e-6;
}

/*
int read_graph_from_file(graph **G, int *startvtx, const char* filename, int wrank) {

	//FILE *fp = fopen (filename, "r");
	//if ( fp == NULL ) {
//		fprintf(stderr, "error to open file %s: %s\n", filename, strerror(errno));
//		return -1;
//	}
	*G = (graph *) malloc(sizeof(graph));
	//fread(startvtx, 1, sizeof(int), fp);
	//fread(&((*G)->ne), 1, sizeof(int64_t), fp);
	//fread(&((*G)->nv), 1, sizeof(int), fp);
	*startvtx=0;
	double poworder=10;
	(*G)->ne=100000000;
	(*G)->nbr = (int *) malloc((*G)->ne*sizeof(int));
	int *tempnbr=(int *) malloc(((*G)->ne)*sizeof(int));
	//fread((*G)->nbr, (*G)->ne, sizeof(int), fp);
	//fread((*G)->firstnbr, (*G)->nv+1, sizeof(int), fp);
	double or;
	int count=0,edge=0,vtxd=1;
	while(edge<(*G)->ne){
	//	vtxd=(*G)->ne/(*G)->nv;	
		while(1){
			vtxd=1;
			or=(double)((double)rand()/(double)RAND_MAX)*poworder;
			for(int j=0;j<or;j++)vtxd*=2;	
			if(1./(double)vtxd>(double)rand()/(double)RAND_MAX) break;
		}
		tempnbr[count]=edge;	
		count++;
		edge+=vtxd;
	}
	(*G)->nv=count;
	(*G)->firstnbr = (int *) malloc(((*G)->nv+1)*sizeof(int));
	
	for(int i=0;i<count;i++){
		(*G)->firstnbr[i]=tempnbr[i];
	}	
	count=0;
	while(count<(*G)->ne){
		(*G)->nbr[count]=rand()%(*G)->nv;	
		count++;
	}
    printf("%ld\t%d\n", (*G)->ne, (*G)->nv);
    FILE *fp = fopen (filename, "w");
	if ( fp == NULL ) {
		fprintf(stderr, "error to open file %s: %s\n", filename, strerror(errno));
		return -1;
	}
	fwrite(startvtx, 1, sizeof(int), fp);
	fwrite(&((*G)->ne), 1, sizeof(int64_t), fp);
	fwrite(&((*G)->nv), 1, sizeof(int), fp);
	fwrite((*G)->nbr, (*G)->ne, sizeof(int), fp);
	fwrite((*G)->firstnbr, (*G)->nv+1, sizeof(int), fp);
	free(tempnbr);
	fclose(fp);
    printf("%ld\t%d\n", (*G)->ne, (*G)->nv);
	return 0;
}
*/
/*
*/
int read_graph_from_file(graph **G, int *startvtx, const char* filename, int wrank) {

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

void print_CSR_graph (const graph *G, int wrank) {
	int vlimit = 20;
	int elimit = 50;
	int e,v;
	Fprintf(stdout, "\nGraph has %d vertices and %"PRId64" edges.\n",G->nv,G->ne);
	Fprintf(stdout, "firstnbr =");
	if (G->nv < vlimit) vlimit = G->nv;
	IF for (v = 0; v <= vlimit; v++) printf(" %d",G->firstnbr[v]);
	if (G->nv > vlimit) Fprintf(stdout, " ...");
	Fprintf(stdout, "\n");
	Fprintf(stdout, "nbr =");
	if (G->ne < elimit) elimit = G->ne;
	IF for (e = 0; e < elimit; e++) printf(" %d",G->nbr[e]);
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
	queue = (int *) malloc(G->nv*sizeof(int));
	int *level = (int *) malloc(G->nv*sizeof(int)); 
	levelsize = *levelsizep = (int *) malloc(G->nv*sizeof(int)); 

    int *gqueue, *bmap, *lqueue;
    IF gqueue = (int *) malloc(G->nv*sizeof(int)*wsize);
	IF bmap = (int *) calloc(G->nv, sizeof(int));
	lqueue = (int *) malloc(G->nv*sizeof(int));
    int *num, *sdex, lback, lfront; 
    num = (int *)malloc(wsize);
    sdex = (int *)malloc(wsize);
    int lvsize;
    double time[10];


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
    for (i=0; i<10;i++)
        time[i]=0;
    time[0] = MPI_Wtime();

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
        front += v;

        lfront = 0;
        lback = num[wrank];

        //MPI_Barrier(MPI_COMM_WORLD);
        // ** local queue - 
        // ** level - Update information of other
        // ** levelsize - sum queue
        time[1] += MPI_Wtime() - time[0];
        time[0] = MPI_Wtime();
        for (i = 0; i < num[wrank]; i++) {
            v = lqueue[lfront++];
			for (e = G->firstnbr[v]; e < G->firstnbr[v+1]; e++) {
				w = G->nbr[e];          // w is the current neighbor of v
				if (level[w] == -1) {   // w has not already been reached
					levelsize[thislevel+1]++;
					lqueue[lback++] = w;    // put w on queue to explore
				}
			}
        }

        time[2] += MPI_Wtime() - time[0];
        time[0] = MPI_Wtime();
        mynum = num[wrank];

        //queue, level size
        MPI_Gather(levelsize+thislevel+1, 1, MPI_INT,\
                num, 1, MPI_INT, 0, MPI_COMM_WORLD);
        IF
        {
            v = 0;
            for(i=0; i<wsize; i++)
            {
                sdex[i] = v;
                v += num[i];
            }
        }
        
        MPI_Gatherv(lqueue+mynum, levelsize[thislevel+1], MPI_INT, gqueue, num, sdex, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Reduce(levelsize+thislevel+1, &lvsize, 1, MPI_INT,\
            MPI_SUM, 0, MPI_COMM_WORLD);
        time[3] += MPI_Wtime() - time[0];
        time[0] = MPI_Wtime();
         
        //Remove Overlap
        //
        IF
        {
            lvsize = v;
            for(i=0; i<v ;i++)
            {
                if(bmap[gqueue[i]]==0)
                {
                    queue[back++]=gqueue[i];
                    bmap[gqueue[i]]++;
                }
                else
                    lvsize--;
            }
            levelsize[thislevel+1] = lvsize;
        }
        MPI_Bcast(levelsize+thislevel+1, 1, MPI_INT, 0, MPI_COMM_WORLD);
        time[4] += MPI_Wtime() - time[0]; //0.6s
        time[0] = MPI_Wtime();

        //level
        MPI_Bcast(&back, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&front, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(queue+front, back-front, MPI_INT, 0, MPI_COMM_WORLD);
        time[5] += MPI_Wtime() - time[0]; // 0.35s
        time[0] = MPI_Wtime();
        for(i=front; i<back; i++)
            level[queue[i]] = thislevel + 1;
        time[6] += MPI_Wtime() - time[0]; // 0.35s
        time[0] = MPI_Wtime();
        // ** Fuse the data in the levels
		thislevel = thislevel+1;
	}
    IF
    {
        for(i=1;i<7;i++)
            printf("Elapsed %f\n", time[i]);
    }

	*nlevelsp = thislevel;
	free(queue);
    // Added
    IF free(gqueue);
    IF free(bmap);
    free(lqueue);
    free(num);
    free(sdex);
	free(level);
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
		fprintf(stderr, "usage: bfs <filename>\n");
		exit(1);
	}
//	starttime = get_time();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    if (read_graph_from_file(&G, &startvtx, filename, wrank)) return -1;

	print_CSR_graph (G, wrank);

	Fprintf(stdout, "Starting vertex for BFS is %d.\n\n",startvtx);
	starttime = MPI_Wtime();
    // Begin the code
	bfs (startvtx, G, &nlevels, &levelsize, wrank, wsize);
    elapsedtime = MPI_Wtime() - starttime;
	//elapsedtime = get_time() -starttime;
    if(wrank==0)
    {
        reached = 0;
        for (i = 0; i < nlevels; i++) reached += levelsize[i];
        fprintf(stdout, "Breadth-first search from vertex %d reached %d levels and %d vertices.\n",
            startvtx, nlevels, reached);
        for (i = 0; i < nlevels; i++) fprintf(stdout, "level %d vertices: %d\n", i, levelsize[i]);
        fprintf(stdout, "\n");
        fprintf(stdout, "GTEPs: %f\n", (reached/elapsedtime)/1000000);
        fprintf(stdout, "Elapsed Time: %f\n", elapsedtime);
    }
	free(levelsize);
	free(G);
    MPI_Finalize();
}
