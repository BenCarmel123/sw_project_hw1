import sys
import math

def euclideanDist(x,y):
    baseDist = [math.pow((x[i] - y[i]),2) for i in range(len(x))]
    return math.sqrt(sum(baseDist))

def calcMin(x,centeroids):
    minDist = min([euclideanDist(x,centroid) for centroid in centeroids])
    return minDist.index()

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
    clusters = [[] * k]
    for vector in vectors:
        i = calcMin(vector, centroids)
        clusters[i].append(vector)
    return clusters

def newCentr(cluster, dim): 
    if len(cluster) < 2:
        return cluster[0]
    new_centroid = [0 * dim]
    size = len(cluster)
    for i in range(dim):
        new_centroid[i] = math.avg(cluster[j][i] for j in range(size))
    return new_centroid 

def updateCentr(vectors, centroids, clusters): # update centroids according to new clusters
    old_centroids, new_centroids = centroids
    flag = False
    for centroid, index in enumerate(centroids):
        new_centroids[index] = newCentr(clusters[index], len(vectors[0]))
    if euclideanDist(old_centroids, new_centroids) > 0.001:
        flag = True
    return (new_centroids, flag)
    
def main(): # main function 
    vectors, k ,iter = readData(sys.argv)
    centroids = initCentr(vectors, k)
    for i in range(iter):
        clusters = assignClusters(vectors, centroids, k)
        centroids, flag = updateCentr(vectors, centroids, clusters)
        if not flag:
            return centroids
    return centroids 

if __name__  == "__main__" :
    main()
            



