# we import the python package numpy to use for our matrix calculation
# "as np" allows us to use the shortcut np instead of numpy
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import special
import matplotlib.animation as animation

pi = np.pi
def RR(x):
    return float(x)

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def arctan(x):
    return np.arctan(x)

#a 3X3 matrix that that rotates the second and third coordinates by rotation_angle
# This will be used to change from the coordinates at the wall to standard coordinates
def rotation(rotation_angle):
    return np.matrix([[-1,0,0],[0,cos(rotation_angle),-sin(rotation_angle)],[0,sin(rotation_angle),cos(rotation_angle)]])


# This gives the collision transformation in coodinates with the normal upward in the direction of the third coordinate
# First coordinate is rotational velocity, second tangential, and third is normal
# This might look unnecessarily complicated Notice we have g=gamma=0, which makes it very simple
def reflection(g):
    return np.matrix([[-(1-g**2)/(1+g**2),-2*g/(1+g**2),0],[-2*g/(1+g**2),(1-g**2)/(1+g**2),0],[0,0,-1]])


# This finds the new position, given by the next collision point
# state=the current state, one row of the states matrix
# g is the mass disribution, always zero for now

def new_pos(state,e):
    #print(state)
    A=RR((1-e)*(state[4])**2+(state[5])**2)
    B=RR(2*(1-e)*(state[1]*state[4])+2*state[2]*state[5])
    C=RR((1-e)*(state[1])**2+(state[2])**2)+(e-1)
    #print(A,B,C)
    t=RR((-B+(B**2-4*A*C)**.5)/(2*A))
    state=[state[0]+t*state[3],state[1]+t*(state[4]),state[2]+t*(state[5]),state[3],state[4],state[5]]
    return state



def normalshift(x,y,e):
    if y>0:
        if x<0:
            return RR(-pi/2-arctan(y/(x*(1-e))))
        if x>0:
            return RR(-arctan(y/(x*(1-e)))+pi/2)
        else:
            return RR(0)
    if y<0:
        if x<0:
            return RR(-pi/2-arctan(y/(x*(1-e))))
        if x>0:
            return RR(-arctan(y/(x*(1-e)))+pi/2)
        else:
            return RR(pi)
    else:
        if x<0:
            return RR(-pi/2)
        else:
            return RR(pi/2)

# This finds the new velocity after the collision
# state=the current state, one row of the states matrix
# g is the mass disribution, always zero for now
# e=eccentricity squared
def new_vel(state,g,e):
    normal_angle=normalshift(state[1],state[2],e)
    #print(normal_angle)
    #OLD:normal_angle=RR(arctan(state[2]/(state[1]*(1-e)))-pi/2)
    #print(normal_angle)
    velocity=np.matrix([ [state[3] ],[state[4] ],[state[5] ] ])
    # To get the new velocity, we rotate to standard position, apply the collision matrix, then rotate back
    velocity_rel=reflection(g)*rotation(normal_angle)*velocity
    velocity=rotation(-normal_angle)*reflection(g)*rotation(normal_angle)*velocity
    state=[state[0],state[1],state[2],velocity[0,0],velocity[1,0],velocity[2,0],velocity_rel[1,0],velocity_rel[2,0]]
    return state

# Each iteration involves finding the new position and then finding the new velocity after the collision
# state=the current state, one row of the states matrix
# g is the mass disribution, always zero for now
# e=eccentricity squared
def F(state, g,e):
    state=new_pos(state, e)

    state=new_vel(state, g,e)
    return state

def arcpart(ec,xypoints,epsilon):
    xylist = [[],[]]
    for x in np.linspace(0.001,0.999,xypoints):
        xylist[0].append(x-epsilon)
        y = ((1-x**2)*(1-ec))**0.5
        xylist[1].append(y)
    return xylist

def arclength(ec,B,X,Y):
    arc = special.ellipeinc(np.pi/2,ec)
    if X == 0 and Y > 0:
        return arc
    elif X == 0 and Y < 0:
        return 3*arc
    if X > 0 and Y >= 0:
        return special.ellipeinc(np.arctan((Y/X)),ec)
    elif X < 0 and Y >= 0:
        return 2*arc + special.ellipeinc(np.arctan((Y/X)),ec)
    elif X < 0 and Y < 0:
        return 2*arc + special.ellipeinc(np.arctan((Y/X)),ec)
    elif X > 0 and Y < 0:
        return 4*arc + special.ellipeinc(np.arctan((Y/X)),ec)


def rotate(angle):
     ax.view_init(elev=5, azim=angle)

# The following two features allow varying the collision model and mass distribution
# collision_type=selector(['Specular','No-Slip'],buttons=True)
# gamma=slider(0,1,label="mass distribution",default=1/2**.5),
# spin=slider(-1,1,default=0, label="angular velocity"),
# (They are not implemented correctly yet, so they are here)

# This allows us to run interactively, changing the parameters

################################################################################
################################################################################
collisions = 200
eclist = [0.15]#[0.15,0.35,0.5,0.75,0.95] # square of eccentricity
xypoints = 10
ANGLEPoints = 10
ANGLEStart = 0.0001
ANGLEEnd = np.pi-0.01
SPINStart = -0.8
SPINEnd = 0.8
SPINPoints = 5
gammalist = [1/(2**0.5)]#[0.1,0.4,1/(2**0.5),0.9,1,2,3]
################################################################################
epsilon = 0.0001
maxpoints = ANGLEPoints*xypoints*SPINPoints
angleList = list(np.linspace(ANGLEStart,ANGLEEnd,ANGLEPoints))


for ec in eclist:
    for gamma in gammalist:
        xylist = arcpart(ec,xypoints,epsilon)
        b = (1-ec)**0.5
        print('focus','('+str(-ec**0.5)+',0)')
        print('focus','('+str(ec**0.5)+',0)')
        halfarc = 2*special.ellipeinc(np.pi/2,ec)
        fname = 'ellipse_ec'+str(round(ec,3))+'_gamma'+str(round(gamma,3))+'_spin'+str(round(SPINStart,2))+'-'+str(round(SPINEnd,2))+'_iter'+str(maxpoints*collisions)
        fname+='_greenred'

        noslip_xyz = [[],[],[],[]]
        color = 1
        count = 0

        for spin in np.linspace(SPINStart,SPINEnd,SPINPoints):
            for x0,y0 in zip(xylist[0],xylist[1]):
                angleList.append(pi+arctan(y0/(x0-ec**0.5)))
                angleList.append(pi+arctan(y0/(x0+ec**0.5)))
                for angle in angleList:
                    vx = cos(angle)
                    vy = sin(angle)
                    E = vx**2 + vy**2 + spin**2
                    spin = spin/(E**0.5)
                    vx = vx/(E**0.5)
                    vy = vy/(E**0.5)

                    if x0**2+y0**2/(1-ec)>1:
                        print("WARNING! WARNING! Outside of the ellipse!")


                    states=np.zeros(shape=(collisions,8))
                    xf=0 #x velocity in frame
                    yf=0 #y velocity in frame
                    states[0]=[0,x0,y0,spin,vx,vy,xf,yf]

                    # Finding where the first trajectory intersects the x-axis
                    s = -y0/vy
                    root = x0+vx*s # intersection point

                    if abs(round(root,12))==round(ec**0.5,12):
                        color = '#00ff00'
                        count += 1
                        print(count)
                        # continue
                    elif abs(root) < ec**0.5:
                        color = '#ff0000'
                        # continue
                    elif abs(root) > ec**0.5:
                        color = '#0000ff'
                        # continue
                    else:
                        print("WARNING")

                    for i in range(0,collisions-1):
                        states[i+1]=F(states[i], gamma, ec)
                        alength = arclength(ec,b,states[i+1][1],states[i+1][2])
                        # if states[i+1][6] < 0:
                        #     continue
                        # if states[i+1][3] < 0:
                        #     continue
                        if alength <= halfarc:
                            noslip_xyz[0].append(states[i+1][6])
                            noslip_xyz[1].append(states[i+1][3])
                            noslip_xyz[2].append(alength)
                            noslip_xyz[3].append(color) #/maxpoints)
                    # color += 1

                angleList.pop()
                angleList.pop()

        fig=plt.figure(figsize=(20,20))

        ################ 2D Spin vs Angle VPP ###########################################
        sc = plt.scatter(noslip_xyz[0],noslip_xyz[1],s=1,c=noslip_xyz[3],cmap='gist_rainbow')
        fname += '_2dspin'
        plt.axis('off')
        ################ 2D Arclength vs Angle VPP ##################################
        # sc = plt.scatter(noslip_xyz[0],noslip_xyz[2],s=1,c=noslip_xyz[3],cmap='gist_rainbow')
        # fname += '_2darcl'
        # plt.axis('off')
        ################ 3D PHASE SPACE ###################################
        # ax = plt.axes(projection = '3d')
        # ax.set_box_aspect(aspect = (1,1,1))
        # # ax.view_init(10,110)
        # fname += '_3d'
        # noslip = ax.scatter(noslip_xyz[0],noslip_xyz[1],noslip_xyz[2],s=1,c=noslip_xyz[3],cmap='gist_rainbow')
        ################ Save #############################################
        plt.savefig(fname+'.png',transparent=True)
        ################ GIF Animation ####################################
        # print('animation')
        # ani = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 360, 9), interval=1000)
        # ani.save(fname+'.gif', writer=animation.PillowWriter(fps=3))
        ################ Show #############################################
        # plt.show()
        plt.close('all')
