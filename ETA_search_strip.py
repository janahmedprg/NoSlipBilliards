import numpy as np
import math
import matplotlib.pyplot as plt

pi = np.pi

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def mat_norm(PHI):
    sqsum = 0
    I = [[1,0,0],[0,1,0],[0,0,1]]
    for i in range(0,3):
        for j in range(0,3):
            sqsum += (float(PHI[i,j])-float(I[i][j]))**2
    return (sqsum**0.5)

def rotate(x):
    return np.matrix([[1,0,0],[0,cos(x),sin(x)],[0,-sin(x),cos(x)]],dtype=float)

def mat_mul(mat_list):
    PHI = [[1,0,0],[0,1,0],[0,0,1]]
    for i in range(0,len(mat_list)):
        PHI = np.matmul(mat_list[i],PHI)
    return PHI

wall_list = [[1,0,1,0,1,0,1,0,1,0]]
ETAStart = 0
ETAEnd = 1
ETAPoints = 10002
min_dis = 9000
min_eta = 0
# check = []
etaVSdis = [[],[]]

for eta in np.linspace(ETAStart,ETAEnd,ETAPoints):
    Tr = [[-cos(eta*pi),-sin(eta*pi),0],[-sin(eta*pi),cos(eta*pi),0],[0,0,-1]]
    mat_list = []
    for i in range(0,len(wall_list[0])):
        # if wall_list[0][i] == 0:
        #     # mat_list.append(rotate(-pi*1/6))
        #     # mat_list.append(rotate(pi))
        #     mat_list.append(Tr)
        #     # mat_list.append(rotate(pi))
        #     # mat_list.append(rotate(pi*1/6))
        # elif wall_list[0][i] == 1:
        #     mat_list.append(rotate(pi*2/3))
        #     mat_list.append(Tr)
        #     mat_list.append(rotate(-pi*2/3))
        # elif wall_list[0][i] == 2:
        #     mat_list.append(rotate(-pi*2/3))
        #     mat_list.append(Tr)
        #     mat_list.append(rotate(pi*2/3))
        if wall_list[0][i] == 1:
            mat_list.append(rotate(pi))
            mat_list.append(Tr)
            mat_list.append(rotate(pi))
            # check.append('R1')
            # check.append('T1')
            # check.append('Ri1')
        elif wall_list[0][i] == 0:
            mat_list.append(Tr)
            # check.append('T0')
    PHI = mat_mul(mat_list)
    dis = mat_norm(PHI)
    if dis<min_dis:
        min_dis = dis
        min_eta = eta
    etaVSdis[0].append(eta)
    etaVSdis[1].append(dis)

plt.clf()
plt.scatter(etaVSdis[0],etaVSdis[1],s=1)
plt.savefig('eta_X_distance_range'+str(ETAStart)+'-'+str(ETAEnd)+'_iter'+str(ETAPoints)+'.png')
##print(str(best_eta) + ' is the best Eta from the run')
print(min_eta)
print('done')
