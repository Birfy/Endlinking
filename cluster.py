import numpy as np
import random
import sys

chainlength = int(sys.argv[1])
dfname = sys.argv[2]
outfl = 'result.data'
cluster_size = int(sys.argv[3])

def readsize(dfname):
    with open(dfname, 'r') as df:
        lines = df.readlines()
    for line in lines:
        content = line.split()
        if content and content[-1] == 'xhi':
            return 2*float(content[1])

def readdata(dfname, chainlen):
    X=[]
    Xi=[]
    with open(dfname, 'r') as df:
        lines = df.readlines()
    for line in lines:
        content = line.split()
        if len(content) == 9:
            # print(content)
            if (int(content[0]) % chainlen == 0 or int(content[0]) % chainlen == 1) and int(content[2]) != 3 and int(content[2]) != 4 :
                X.append([float(content[i]) for i in range(3,6)])
                Xi.append(int(content[0]))
    return np.array(X), np.array(Xi)

def initmeans(n):
    M=[]
    for i in range(n):
        M.append([size*(random.random()-0.5),size*(random.random()-0.5),size*(random.random()-0.5)])
    return np.array(M)

def SetDistMat(X, means):
    distmat_dtype = [('key',int), ('dist',float)]
    distmat = np.empty((n,k),dtype=distmat_dtype)
    for i in range(n):
        distmat[i,:] = [(c[0], GetDist(X[i], c[1])) for c in enumerate(means)]
        distmat[i,:] = np.sort(distmat[i,:], order='dist')
    return distmat

def GetDist(x, c):
    dist = np.linalg.norm(x-c-boxl*np.around((x-c)/boxl))
    return dist

def Get_plst(assigned, distmat, full):
    plst = []
    for i in range(n):
        if (i not in assigned):
            j = 0
            while j<k:
                if (not full[distmat[i,j][0]]):
                    bestkey = distmat[i,j][0]
                    mindist = distmat[i,j][1]
                    break
                else:
                    j += 1
            for j in range(k-1,-1,-1):
                if (not full[distmat[i,j][0]]):
                    maxdist = distmat[i,j][1]
                    break
            plst.append((i, bestkey, maxdist-mindist))
    plst.sort(key=lambda t:t[2])
    return plst

def InitialAssignment(distmat):
    clusters = {}
    full = np.zeros(k,dtype=bool) # a boolean array that records which clusters are full
    assigned = [] # a list of objects who has been assigned to a cluster
    plst = Get_plst(assigned, distmat, full)
    while (len(plst)):
        temp = plst.pop()
        try:
            if (len(clusters[temp[1]])<cluster_size):
                clusters[temp[1]].append(temp[0])
                assigned.append(temp[0])
            else:
                full[temp[1]] = True
                plst = Get_plst(assigned, distmat, full)
        except KeyError:
            clusters[temp[1]] = [temp[0]]
            assigned.append(temp[0])
    return clusters

def CalcMeans(X, oldmeans, clusters):
    means = np.zeros((k,3))
    keys = sorted(clusters.keys())
    for key in keys:
        for i in clusters[key]:
            means[key] += X[i]-boxl*np.around((X[i]-oldmeans[key])/boxl)
        means[key] /= len(clusters[key])
        means[key] -= boxl*np.around(means[key]/boxl)
    return means

def SortObj(X, clusters, means, distmat):
    objlst = [] # list of objects ordered in asceding delta of the current 
                # assignment and the best possible alternate assignment
    keys = sorted(clusters.keys())
    for key in keys:
        for i in clusters[key]:
            currdist = GetDist(X[i],means[key])
            mindist = distmat[i,0][1]
            objlst.append((i, key, currdist-mindist))
    objlst.sort(key=lambda t:t[2], reverse=True)
    return objlst

def Transfer(obj, clufrom, cluto, clusters):
    clusters[clufrom].remove(obj)
    clusters[cluto].append(obj)
    return clusters

def WriteResult(file, X, means, clusters):
    with open(file, 'w') as fl:
        # keys = sorted(clusters.keys())
        # i = 1
        # for key in keys:
        #     for obj in clusters[key]:
        #         fl.write("%d\t%d\t%f\t%f\t%f\t%d\n"\
        #                 %(obj,Xi[obj], X[obj][0], X[obj][1], X[obj][2], key)) 
                # i = i + 1
        for c in enumerate(means):
            fl.write("%d\t%f\t%f\t%f"%(c[0], c[1][0], c[1][1], c[1][2]))
            for obj in clusters[c[0]]:
                fl.write("\t%d"%(Xi[obj]))
            fl.write('\n') 
            # i = i + 1
    return

# This function will perform statistical analysis to the clustering results
def ClusterStat(X, means, clusters):
    # Average distance between means
    means_avg = 0.
    for i in range(k-1):
        for j in range(i+1,k):
            means_avg += GetDist(means[i], means[j])
    means_avg /= (k*(k-1)/2.)
    # Average distance between obj and mean in a cluster
    obj2mean_avg = np.zeros(k)
    # Variance of the distances between obj and mean in a cluster
    obj2mean_var = np.zeros(k)
    keys = sorted(clusters.keys())
    for key in keys:
        for i in clusters[key]:
            obj2mean = GetDist(X[i], means[key])
            obj2mean_avg[key] += obj2mean 
            obj2mean_var[key] += obj2mean*obj2mean
        obj2mean_avg[key] /= len(clusters[key])
        obj2mean_var[key] /= len(clusters[key])
        obj2mean_var[key] = np.sqrt(obj2mean_var[key])
    # Average within cluster distances between objects
    winclu_avg = np.zeros(k)
    # Average of within cluster distances of all clusters
    winclu_grandavg = 0.
    for key in keys:
        for i in clusters[key]:
            x = X[i]
            for j in clusters[key]:
                if j>i:
                    winclu_avg[key] += GetDist(x, X[j])
        s = len(clusters[key])
        winclu_avg[key] /= (s*(s-1)/2) 
        winclu_grandavg += winclu_avg[key]
    winclu_grandavg /= k
    # write the summary 
    print("average distance among means: %f"%means_avg)
    #print("average distance from objects to the mean of a cluster:")
    #for i in range(k):
    #    print("cluster %i: %f"%(i, obj2mean_avg[i]))
    #print("variance of distances from objects to the mean of a cluster:")
    #for i in range(k):
    #    print("cluster %i: %f"%(i, obj2mean_var[i]))
    #print("within-cluster average distances:")
    #for i in range(k):
    #    print("cluster %i: %f"%(i, winclu_avg[i]))
    print("grand average of within-cluster average distances: %f"%winclu_grandavg)
    return 



X, Xi = readdata(dfname, chainlength)
size = readsize(dfname)
boxl = np.array([size, size, size])


n = len(X)
k = int(len(X)/cluster_size)

# Set up the database of objects
# X = readdata(dfname, chainlength)
# Choose initial means with K-means
means = initmeans(k)
# Set up initial clusters
distmat = SetDistMat(X, means) 
clusters = InitialAssignment(distmat) 
## debug code
#keys = sorted(clusters.keys())
#for key in keys:
#    print("cluster %i:"%key)
#    print(clusters[key])
## end of debug
# Iteration step
for iter in range(100):
    active = 0 # indicate the number of transfers in the current iteration
    tranlst = (-1)*np.ones(k, dtype='int') # set up transfer list for each cluster
    # Compute the cluster means
    oldmeans = means.copy()
    means = CalcMeans(X, oldmeans, clusters)
    # Get statistics about the clustering
    #ClusterStat(X, means, clusters)
    ## debug code
    #print("old means:")
    #print(oldmeans)
    #print("new means:")
    #print(means)
    ## end of debug
    # For each object, compute the distances to the cluster means
    distmat = SetDistMat(X, means)
    # Sort objects based on the delta of the current assignment and the best 
    # possible alternate assignment
    objlst = SortObj(X, clusters, means, distmat)
    ##debug code
    #print(objlst)
    ##return
    #end of debug
    # For each element by prioty:
    while (len(objlst)):
        (i, key, temp) = objlst.pop()
        obj2key = GetDist(X[i], means[key])
        transferred = False #record if any transfering has occured to i 
        if (key == distmat[i,0][0]):
            ##debug
            #print("%i is already the opt cluster for obj %i. no transfer"%(clu, i))
            ##end of debug
            continue
        # For each other clusters by element gain:
        else:
            for j in range(k):
                clu = distmat[i,j][0] # the key of another cluster
                objgain = obj2key - distmat[i,j][1] # gain by transfering i from cluster key to clu
                if (clu==key): # already in the cluster
                    continue
                if (len(clusters[clu]) < cluster_size):
                    active += 1
                    transferred = True
                    clusters = Transfer(i, key, clu, clusters)
                    ##debug
                    #print("cluster %i not full. transfer obj %i from cluster %i to it."%(clu, i, key))
                    ##end of debug
                    break
                elif (tranlst[clu] != -1): # if the tranlst of another cluster is not empty
                    # distance between the obj in the tranlst and the current cluster
                    tran2key = GetDist(X[tranlst[clu]], means[key])
                    tran2clu = GetDist(X[tranlst[clu]], means[clu])
                    # gain by transfering the obj in tranlst from cluster clu to key
                    trangain = tran2clu - tran2key
                    if (objgain + trangain > 0): # transfer if the sum of gains are positive, ie net gain
                        active += 2
                        transferred = True
                        clusters = Transfer(i, key, clu, clusters)
                        clusters = Transfer(tranlst[clu], clu, key, clusters)
                        ##debug
                        #print("obj %i is transfered from cluster %i to %i"%(i, key, clu))
                        #print("obj %i is transfered from cluster %i to %i"%(tranlst[clu], clu, key))
                        #print("objgain: %f, trangain: %f"%(objgain, trangain))
                        ##end of debug
                        tranlst[clu] = -1 # reset the tranlst to empty
                        break
            if (not transferred):
                tranlst[key] = i
                ##debug
                #print("add obj %i in cluster %i to the transfer list"%(i, key))
                ##end of debug
    # nothing is transferred during this iteration, return the clustering result
    if (not active):
            break
    #debug code
    print("number of transfers in iter %i: %i\n"%(iter+1, active))
    #end of debug
print("K-means clustering converged in %d iterations!\n"%(iter+1))
# Output the clustering results
WriteResult(outfl, X, means, clusters)
ClusterStat(X, means, clusters)
# print(X)
