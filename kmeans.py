import sys
import math

def euclideanDist(x,y): #Calculate euclidean distance between 2 vectors
    return math.sqrt(sum((xi - yi) ** 2 for xi, yi in zip(x, y)))

def calcMin(x,centroids): # Find the closest centoid to vector X
    distances = [euclideanDist(x, centroid) for centroid in centroids]
    return distances.index(min(distances))

def readData(args): # read info from terminal
    k = int(args[1])
    iter = 400 if len(args) == 2 else int(args[2]) #Check if iter is given or define as 400
    if iter not in range (2,1000): 
        print('Incorrect maximum iteration!')
        sys.exit(1) #Error 2 - stop program
    vectors = [] #Denote an empty list which will contain all vectors
    for line in sys.stdin:
        if line.strip():
            vector = [float(num) for num in line.strip().split(',')] #Convert each line to a vector (list)
            vectors.append(vector) # add to vectors list
    if len(vectors) < k: 
        print('Incorrect number of clusters!')
        sys.exit(1) #Error 1 - stop program
    return (vectors, k, iter)

def initCentr(vectors,k): # initalize centroids
    return [vectors[i] for i in range(k)]
    
def assignClusters(vectors, centroids, k): # derive new clusters from centroids
    clusters = [[] for _ in range(k)] # Denote a list with k bins which will contain all relevant vectors to the cluster
    for vector in vectors:
        i = calcMin(vector, centroids)
        clusters[i].append(vector) # Check which centroid is the closest and assign the vector
    return clusters

def newCentr(cluster, dim): # Calculate the new average (centroid value) after assignments
    if len(cluster) < 2: # If there is only 1 vector return it
        return cluster[0]
    new_centroid = [0] * dim # Denote with zeros 
    size = len(cluster)
    for i in range(dim): # For each 'intro' of the vector calc the new average
        new_centroid[i] = (sum(cluster[j][i] for j in range(size)))/size 
    return new_centroid # Return the new centroid according to new vectors assignment

def updateCentr(vectors, centroids, clusters): # update centroids according to new clusters
    old_centroids = centroids
    new_centroids = []
    flag = False # Boolean which signs if has converged or not. if True we continue running, if False we finished.
    for index, centroid in enumerate(centroids):
        new_centroids.append(newCentr(clusters[index], len(vectors[0]))) # Fill new_centroids with new values according to calculation
    flag = not has_converged(old_centroids, new_centroids) # Update Flag bool
    return (new_centroids, flag)

def has_converged(old_centroids, new_centroids): # Checks if the euclidean distance of all vectors is < 0.001
    for old, new in zip(old_centroids, new_centroids):
        if euclideanDist(old, new) >= 0.001:
            return False
    return True

def print_cents(centroids): #Helper function for printing centroids
    for centroid in centroids:
        print(",".join(f"{float(x):.4f}" for x in centroid))
    
def main(): # main function 
    vectors, k ,iter = readData(sys.argv) # Read 
    centroids = initCentr(vectors, k) # Initialize centroids as first k vectors
    for i in range(iter): # Update and check centroids num_of_iterations times
        clusters = assignClusters(vectors, centroids, k) # Assign clusters according to changes
        centroids, flag = updateCentr(vectors, centroids, clusters) # Update centroids values (averages)
        if not flag: # Break if has converged before max iterations number
            break
    print_cents(centroids) # Print centroids
    
if __name__  == "__main__" :
    main()

