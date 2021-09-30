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
#[1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0],[1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0,1,0,2,0],

################################################
################################################
wall_list = [[1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2],[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]]#,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2]]#,[2,0,1,2,1,0,2,1,0,1,2,0,1,2,1,0,2,1]]
for ii in range(0,len(wall_list)):
    print(len(wall_list[ii]))
ETAStart = 0.05
ETAEnd = 0.95
ETAPoints = 5000
wedge_Langle = 120
wedge_Rangle = 120
################################################
################################################
etaVSdis = [[],[],[]]
color = 1
lbls =[]
eta_check = []
etastr ='Eta candidates:\n'
phi_Rangle = math.radians(wedge_Rangle)
phi_Langle = math.radians(wedge_Langle)

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
                # check.append('T0')
            elif wall_list[combin][i] == 1:
                mat_list.append(rotate(-phi_Rangle+pi))
                mat_list.append(Tr)
                mat_list.append(rotate(phi_Rangle-pi))
                # check.append('R1')
                # check.append('T1')
                # check.append('Ri1')
            elif wall_list[combin][i] == 2:
                mat_list.append(rotate(phi_Langle+pi))
                mat_list.append(Tr)
                mat_list.append(rotate(-phi_Langle-pi))
                # check.append('R1')
                # check.append('T1')
                # check.append('Ri1')

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

    if min_eta != ETAStart and min_eta != ETAEnd and repeat == False and min_dis<0.3:
        print(min_eta)
        etastr += str(min_eta)+'\n'
        eta_check.append(min_eta)
    color += 1

print(eta_check)
sc=plt.scatter(etaVSdis[0],etaVSdis[1],c=etaVSdis[2],s=1,cmap='gist_rainbow')
plt.legend(handles=sc.legend_elements()[0], labels=lbls,bbox_to_anchor=(1.04,1), loc="upper left")
# plt.figtext(0,-0.1, etastr, fontsize=10)
plt.savefig('truncated_wedge_ES_range'+'LAng'+str(round(wedge_Langle,1))+'RAng'+str(round(wedge_Rangle,1))+'_'+str(ETAStart)+'-'+str(ETAEnd)+'_iter'+str(ETAPoints)+'.png',bbox_inches = "tight")
plt.show()
print('done')
