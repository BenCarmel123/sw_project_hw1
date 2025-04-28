import sys
import math

def euclideanDist(x,y):
    return math.sqrt(sum((xi - yi) ** 2 for xi, yi in zip(x, y)))

def calcMin(x,centroids):
    distances = [euclideanDist(x, centroid) for centroid in centroids]
    return distances.index(min(distances))

def readData(args): # read info from terminal
    k = int(args[1])
    iter = 400 if len(args) == 2 else int(args[2])
    if iter not in range (2,1000): 
        print('Incorrect maximum iteration!')
        return
    vectors = []
    for line in sys.stdin:
        if line.strip():
            vector = [float(num) for num in line.strip().split(',')]
            vectors.append(vector)
    if len(vectors) < k: 
        print('Incorrect number of clusters!')
        return
    return (vectors, k, iter)

def initCentr(vectors,k): # initalize centroids
    return [vectors[i] for i in range(k)]
    
def assignClusters(vectors, centroids, k): # derive new clusters from centroids
    clusters = [[] for _ in range(k)]
    for vector in vectors:
        i = calcMin(vector, centroids)
        clusters[i].append(vector)
    return clusters

def newCentr(cluster, dim): 
    if len(cluster) < 2:
        return cluster[0]
    new_centroid = [0] * dim
    size = len(cluster)
    for i in range(dim):
        new_centroid[i] = (sum(cluster[j][i] for j in range(size)))/size 
    return new_centroid 

def updateCentr(vectors, centroids, clusters): # update centroids according to new clusters
    old_centroids = centroids
    new_centroids = []
    flag = False
    for index, centroid in enumerate(centroids):
        new_centroids.append(newCentr(clusters[index], len(vectors[0])))
    flag = not has_converged(old_centroids, new_centroids)
    return (new_centroids, flag)

def has_converged(old_centroids, new_centroids):
    for old, new in zip(old_centroids, new_centroids):
        if euclideanDist(old, new) >= 0.001:
            return False
    return True

def print_cents(centroids):
    for centroid in centroids:
        print(",".join(f"{float(x):.4f}" for x in centroid))
    
def main(): # main function 
    vectors, k ,iter = readData(sys.argv)
    centroids = initCentr(vectors, k)
    for i in range(iter):
        clusters = assignClusters(vectors, centroids, k)
        centroids, flag = updateCentr(vectors, centroids, clusters)
        if not flag:
            break
    print_cents(centroids)
    

if __name__  == "__main__" :
    main()

            



