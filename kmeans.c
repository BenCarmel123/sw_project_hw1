#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LEN 1024
#define MAX_DIM 64

double euclideanDist(double* x, double*y, int dim) {
    double sum = 0.0;
    int i;
    for (i = 0; i < dim; i++)
    {
        sum += pow(x[i] - y[i], 2);
    }

    return sqrt(sum);
}

int clacMin(double* x, double** centroids, int k, int dim)
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

double** readData(int argc, char* argv[], int* k, int* iter, int* numVectors, int* dim) {
    char line[MAX_LINE_LEN];
    int capacity = 100;
    int vectorCount = 0;

    /* Handle command line args */
    *k = atoi(argv[1]);
    *iter = (argc == 2) ? 400 : atoi(argv[2]);

    if (*iter < 2 || *iter >= 1000) {
        printf("Incorrect maximum iteration!\n");
        exit(1);
    }

    /* Allocate initial space for vectors */
    double** vectors = malloc(capacity * sizeof(double*));
    if (!vectors) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '\n' || line[0] == '\0') continue; /* skip empty lines */

        /* Allocate space for 1 vector */
        double* vector = malloc(MAX_DIM * sizeof(double));
        if (!vector) {
            printf("Memory allocation failed.\n");
            exit(1);
        }

        int i = 0;
        char* token = strtok(line, ",");
        while (token && i < MAX_DIM) {
            vector[i++] = atof(token);
            token = strtok(NULL, ",");
        }

        if (vectorCount == 0) {
            *dim = i; /* save vector dimensionality */
        }

        if (i != *dim) {
            printf("Inconsistent vector dimension.\n");
            exit(1);
        }

        if (vectorCount == capacity) {
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
        printf("Incorrect number of clusters!\n");
        exit(1);
    }

    *numVectors = vectorCount;
    return vectors;
}

double** initCentr(double** vectors, int k, int dim) {
    double** centroids = malloc(k * sizeof(double*));
    int i;
    for (i = 0; i < k; i++) {
        centroids[i] = malloc(dim * sizeof(double));
        if (!centroids[i]) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
        memcpy(centroids[i], vectors[i], dim * sizeof(double)); /* copy vector */
    }
    return centroids;
}

void assignClusters(
    double** vectors, int numVectors,
    double** centroids, int k, int dim,
    double**** clustersOut, int** clusterSizesOut
) {
    /* Allocate array of cluster arrays */
    double*** clusters = malloc(k * sizeof(double**));
    int* clusterSizes = calloc(k, sizeof(int));
    int* clusterCaps = malloc(k * sizeof(int));

    int i;
    for (i = 0; i < k; i++) {
        clusterCaps[i] = 10;
        clusters[i] = malloc(clusterCaps[i] * sizeof(double*));
    }

    int v;
    for (v = 0; v < numVectors; v++) {
        int closest = clacMin(vectors[v], centroids, k, dim);

        if (clusterSizes[closest] >= clusterCaps[closest]) {
            clusterCaps[closest] *= 2;
            clusters[closest] = realloc(clusters[closest], clusterCaps[closest] * sizeof(double*));
        }

        clusters[closest][clusterSizes[closest]++] = vectors[v];
    }

    *clustersOut = clusters;
    *clusterSizesOut = clusterSizes;

    free(clusterCaps);
}

double* newCentr(double** cluster, int size, int dim) {
    int i, j;
    if (size < 2) {
        /* Return a copy of the only vector */
        double* centroid = malloc(dim * sizeof(double));
        for (i = 0; i < dim; i++) {
            centroid[i] = cluster[0][i];
        }
        return centroid;
    }

    /* Allocate new centroid */
    double* new_centroid = calloc(dim, sizeof(double)); /* zero-initialized */

    for (i = 0; i < dim; i++) {
        for (j = 0; j < size; j++) {
            new_centroid[i] += cluster[j][i];
        }
        new_centroid[i] /= size;
    }

    return new_centroid;
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

int updateCentr(
    double*** clusters,    /* array of k clusters */
    int* clusterSizes,     /* number of vectors in each cluster */
    double*** centroids,   /* pointer to current centroids (updated in-place) */
    int k,                 /* number of clusters */
    int dim                /* vector dimension */
) {
    double** old_centroids = *centroids;
    double** new_centroids = malloc(k * sizeof(double*));

    if (!new_centroids) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    int i;
    for (i = 0; i < k; i++) {
        new_centroids[i] = newCentr(clusters[i], clusterSizes[i], dim);
    }

    int flag = !has_converged(old_centroids, new_centroids, k, dim);

    for (i = 0; i < k; i++) {
        free(old_centroids[i]);
    }
    free(old_centroids);
    *centroids = new_centroids;

    return flag;
}

void print_cents(double** centroids, int k, int dim) {
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
    double** vectors = readData(argc, argv, &k, &iter, &numVectors, &dim);
    double** centroids = initCentr(vectors, k, dim);

    for (i = 0; i < iter; i++) {
        double*** clusters;
        int* clusterSizes;
        int flag;

        assignClusters(vectors, numVectors, centroids, k, dim, &clusters, &clusterSizes);
        flag = updateCentr(clusters, clusterSizes, &centroids, k, dim);

        for (c = 0; c < k; c++) {
            free(clusters[c]);
        }
        free(clusters);
        free(clusterSizes);

        if (!flag) break;
    }

    print_cents(centroids, k, dim);

    for (i = 0; i < numVectors; i++) free(vectors[i]);
    free(vectors);

    for (i = 0; i < k; i++) free(centroids[i]);
    free(centroids);

    return 0;
}
