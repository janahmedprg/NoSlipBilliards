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

# [1,0,1,0],[1,0,2,0],[1,0,3,0],[1,0,4,0],[1,2,3,0],[1,2,4,0],[1,2,1,0],[1,3,2,0],[1,3,4,0],[1,3,1,0],[2,1,4,0],[2,1,2,0],[2,1,3,0],[2,0,2,0],[2,0,1,0],[2,0,3,0],[2,0,4,0]
# [2,0,2,0,2,0],[1,0,1,0,1,0],[1,2,3,2,1,0]

################################################
################################################
wall_list = [[2,0,2,0,2,0],[1,0,1,0,1,0],[1,2,3,2,1,0]]
ETAStart = 0.01
ETAEnd = 0.99
ETAPoints = 2000
################################################
################################################
etaVSdis = [[],[],[]]
color = 1
lbls =[]
eta_check = []
etastr ='Eta candidates:\n'

for i in range(0,len(wall_list)):
    lbls.append(str(wall_list[i]))
for combin in range(0,len(wall_list)):
    min_eta = 0
    min_dis = 9000
    repeat = False
    for eta in np.linspace(ETAStart,ETAEnd,ETAPoints):
        Tr = [[-cos(eta*pi),-sin(eta*pi),0],[-sin(eta*pi),cos(eta*pi),0],[0,0,-1]]
        mat_list = []
        for i in range(0,len(wall_list[combin])):
            if wall_list[combin][i] == 0:
                mat_list.append(Tr)
            elif wall_list[combin][i] == 1:
                mat_list.append(rotate(pi*2/5))
                mat_list.append(Tr)
                mat_list.append(rotate(-pi*2/5))
            elif wall_list[combin][i] == 2:
                mat_list.append(rotate(-pi*6/5))
                mat_list.append(Tr)
                mat_list.append(rotate(pi*6/5))
            elif wall_list[combin][i] == 3:
                mat_list.append(rotate(pi*6/5))
                mat_list.append(Tr)
                mat_list.append(rotate(-pi*6/5))
            elif wall_list[combin][i] == 4:
                mat_list.append(rotate(-pi*2/5))
                mat_list.append(Tr)
                mat_list.append(rotate(pi*2/5))

        PHI = mat_mul(mat_list)
        dis = mat_norm(PHI)
        if dis<min_dis:
            min_dis = dis
            min_eta = eta
        etaVSdis[0].append(eta)
        etaVSdis[1].append(dis)
        etaVSdis[2].append(color/len(wall_list))

    for i in range(0,len(eta_check)):
        if min_eta == eta_check[i]:
            repeat = True

    if min_eta != ETAStart and min_eta != ETAEnd and repeat == False:
        print(min_eta)
        etastr += str(min_eta)+'\n'
        eta_check.append(min_eta)
    color += 1

sc=plt.scatter(etaVSdis[0],etaVSdis[1],c=etaVSdis[2],s=1,cmap='gist_rainbow')
plt.legend(handles=sc.legend_elements()[0], labels=lbls,bbox_to_anchor=(1.04,1), loc="upper left")
# plt.figtext(0,-0.1, etastr, fontsize=10)
plt.savefig('P_5_ES_range'+str(ETAStart)+'-'+str(ETAEnd)+'_iter'+str(ETAPoints)+'.png',bbox_inches = "tight")
##print(str(best_eta) + ' is the best Eta from the run')
print('done')
