#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAX_LINE_LEN 1024
#define MAX_DIM 64

double euclideanDist(double* x, double*y, int dim) { /* Calculate Euclidean distance between 2 Vectors */
    double sum = 0.0;
    int i;
    for (i = 0; i < dim; i++)
    {
        sum += pow(x[i] - y[i], 2);
    }
    return sqrt(sum);
}

int calcMin(double* x, double** centroids, int k, int dim) /* Find the closest centroid to Vector X */
{
    int minIndex = 0;
    double minDist = euclideanDist(x, centroids[0], dim);
    int i;
    for (i = 1; i < k; i++)
    {
        double dist = euclideanDist(x, centroids[i], dim);
        if (dist < minDist) {
            minDist = dist;
            minIndex = i;
        }
    }
    return minIndex;
}

/* Split the String like Strtok (include string.h is not allowed) */
char* my_string(char* str, char param) {
    static char* next = 0;
    char* t;
    
    if (str) next = str;
    if (!next || !*next) return 0;
    t = next;
    
    while (*next && *next != param) next++;
    
    if (*next) 
    {
        *next = '\0';
        next++;
    }
    
    return t;
}

/* Function to copy memory (replacement for memcpy since string.h isn't allowed) */
void my_memory(void* dest, const void* src, int n) {
    int i;
    char* d = (char*)dest;
    const char* s = (const char*)src;
    
    for (i = 0; i < n; i++) 
    {
        d[i] = s[i];
    }
}
/* Read info from Terminal */
double** readData(int argc, char* argv[], int* k, int* iter, int* numVectors, int* dim) {
    char line[MAX_LINE_LEN];
    int capacity = 100;
    int vectorCount = 0;
    double** vectors;
    double* vector;
    char* token;
    int i;

    *k = atoi(argv[1]);
    *iter = (argc == 2) ? 400 : atoi(argv[2]); /* Check if iter is given or define as 400 */

    if (*iter < 2 || *iter >= 1000) {
        printf("Incorrect maximum iteration!\n");
        exit(1); /*Error 2 - Stop Program*/
    }

    vectors = malloc(capacity * sizeof(double*)); /* Allocate memory to store all vectors*/
    
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '\n' || line[0] == '\0') continue; /* skip empty lines */
        
        vector = malloc(MAX_DIM * sizeof(double));
        if (!vector) {
            printf("Memory allocation failed.\n");
            exit(1);
        }

        i = 0;
        token = my_string(line, ','); /* Split line by commas */
        while (token && i < MAX_DIM) { /* Convert each line to a Vector */
            vector[i++] = atof(token);
            token = my_string(NULL, ',');
        }

        if (vectorCount == capacity) { /* Check if we need to reallocate memory */
            capacity *= 2;
            vectors = realloc(vectors, capacity * sizeof(double*));
            if (!vectors) {
                printf("Memory reallocation failed.\n");
                exit(1);
            }
        }

        vectors[vectorCount++] = vector;
    }

    if (vectorCount < *k) {
        printf("Incorrect number of clusters!\n"); /*Error 1 - Stop Program*/
        exit(1);
    }

    *numVectors = vectorCount;
    return vectors;
}

double** initCentr(double** vectors, int k, int dim) { /* Initialize centroids */
    int i;
    double** centroids;
    centroids = malloc(k * sizeof(double*)); /* Allocate memory for centroids */
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(dim * sizeof(double));
        my_memory(centroids[i], vectors[i], dim * sizeof(double)); /* Copy vector */
    }
    return centroids;
}

void assignClusters(double** vectors, int numVectors, double** centroids, int k, int dim, double**** clustersOut, int** clusterSizesOut) 
/* Derive new clusters from centroids*/
{
    double*** clusters;
    int* clusterSizes;
    int* clusterCaps;
    int i;
    int v;
    clusters = malloc(k * sizeof(double**)); /* Denote a list with k bins which all contain all relevant vectors to the cluster*/
    clusterSizes = calloc(k, sizeof(int));
    clusterCaps = malloc(k * sizeof(int));
    
    for (i = 0; i < k; i++) 
    {
        clusterCaps[i] = 10;
        clusters[i] = malloc(clusterCaps[i] * sizeof(double*));
    }

    for (v = 0; v < numVectors; v++) /* Find closest Vectors*/
    {
        int closest;
        closest = calcMin(vectors[v], centroids, k, dim);

        if (clusterSizes[closest] >= clusterCaps[closest]) 
        {
            clusterCaps[closest] *= 2;
            clusters[closest] = realloc(clusters[closest], clusterCaps[closest] * sizeof(double*));
        }

        clusters[closest][clusterSizes[closest]++] = vectors[v];
    }

    *clustersOut = clusters;
    *clusterSizesOut = clusterSizes;

    free(clusterCaps);
}





double* newCentr(double** cluster, int size, int dim) { /* Calculate the new Average (Centroid value) after assignments */
    int i, j;
    double* centroid;
    double* new_centroid;
    if (size < 2) { /* If there is only one vector return it*/
        centroid = malloc(dim * sizeof(double));
        for (i = 0; i < dim; i++) {
            centroid[i] = cluster[0][i];
        }
        return centroid;
    }

    /* Allocate new centroid */
    new_centroid = calloc(dim, sizeof(double)); /* denote with zeros */

    for (i = 0; i < dim; i++) { /* For each 'intro' of the vector calc the new average */
        for (j = 0; j < size; j++) {
            new_centroid[i] += cluster[j][i];
        }
        new_centroid[i] /= size;
    }

    return new_centroid; /* Return the new centroid according to new vectors assignment*/
}

int has_converged(double** old_centroids, double** new_centroids, int k, int dim) {
    int i;
    for (i = 0; i < k; i++) {
        if (euclideanDist(old_centroids[i], new_centroids[i], dim) >= 0.001) {
            return 0;
        }
    }
    return 1;
}

int updateCentr(double*** clusters, int* clusterSizes, double*** centroids, int k, int dim) {
    /* Update centroids according to new clusters */
    double** old_centroids;
    double** new_centroids;
    int i;
    int flag;
    old_centroids = *centroids;
    new_centroids = malloc(k * sizeof(double*)); /* Allocate memory for new centroids */

    for (i = 0; i < k; i++) { /* For each cluster calculate the new centroid */
        new_centroids[i] = newCentr(clusters[i], clusterSizes[i], dim);
    }

    flag = !has_converged(old_centroids, new_centroids, k, dim); /* Check if the centroids have converged */

    for (i = 0; i < k; i++) { /* Free old centroids */
        free(old_centroids[i]);
    }
    free(old_centroids);
    *centroids = new_centroids;

    return flag;
}

void print_cents(double** centroids, int k, int dim) { /* Helper function for printing centroids*/
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < dim; j++) {
            printf("%.4f", centroids[i][j]);
            if (j < dim - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int k, iter, numVectors, dim;
    int i, c;
    double** vectors;
    double** centroids;
    vectors = readData(argc, argv, &k, &iter, &numVectors, &dim); /* Read */
    centroids = initCentr(vectors, k, dim); /* Initialize centroids as first k Vectors */

    for (i = 0; i < iter; i++) { /* Update and check centroids num_of_iteration times*/
        double*** clusters;
        int* clusterSizes;
        int flag;

        assignClusters(vectors, numVectors, centroids, k, dim, &clusters, &clusterSizes); /* Assign Clusters according to changes */
        flag = updateCentr(clusters, clusterSizes, &centroids, k, dim); /* Update centroids values (averages) */

        for (c = 0; c < k; c++) { /* Free clusters */
            free(clusters[c]);
        }
        free(clusters);
        free(clusterSizes);

        if (!flag) break;
    }

    print_cents(centroids, k, dim); /* Print */

    return 0;
}
