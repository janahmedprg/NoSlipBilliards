# A very general billiards program capable of handling specular, no-slip, and rolling billiards

# Imports
from flask import send_file
from io import BytesIO
import  functions_framework
import matplotlib.pyplot as plt
import numpy as np
import math


# Constants
large=10000000
res = 100
eps=0.0000000001
pi = np.pi

def cos(x):
    return np.cos(x)

def sin(x):
    return np.sin(x)

def arccos(x):
    return np.arccos(x)

def arcsin(x):
    return np.arcsin(x)

# an old rotation matrix
def rot(x):
    return np.matrix([[-1,0,0],[0,cos(x),-sin(x)],[0,sin(x),cos(x)]],dtype=float)

# basic functions concerning wall segments
def qb(wall_segment):
    return (wall_segment[0],wall_segment[1])
def qe(wall_segment):
    return (wall_segment[2],wall_segment[3])
def diff(u1,u2):
    return (u2[0]-u1[0],u2[1]-u1[1])

def innprod(u1,u2):
    return u1[0]*u2[0]+u1[1]*u2[1]

# this returns the (smallest) angle between two vectors
def vecang(u1,u2):
    if innprod(u1,u2)/(veclen(u1)*veclen(u2))>1:
        #print(innprod(u1,u2)/(veclen(u1)*veclen(u2)))
        return 0.000000000001
    if innprod(u1,u2)/(veclen(u1)*veclen(u2))<-1:
        #print(innprod(u1,u2)/(veclen(u1)*veclen(u2)))
        return np.pi-0.000000000001

    return arccos(innprod(u1,u2)/(veclen(u1)*veclen(u2)))

def ub(wall_segment):
    w=(wall_segment[4],wall_segment[5])
    l=veclen(w)
    return (wall_segment[4]/l,wall_segment[5]/l)
# Testing simple functions above:
# show(veclen(diff(qb(easy_walls[1]),qe(easy_walls[1]))))

# normal vector from qb and qe
def normal(u1,u2):
    w=diff(u1,u2)
    l=veclen(w)
    return (-w[1]/l,w[0]/l)

#show(normal(qb(easy_walls[1]),qe(easy_walls[1])))
def curv(u1,u2,u3):
    # u1=qb, u2=qe, u3=ub
    norm=normal(u1,u2)
    return 2*(u3[0]*norm[0]+u3[1]*norm[1])/veclen(diff(u1,u2))
#show(curv(qb(easy_walls[1]),qe(easy_walls[1]),ub(easy_walls[1])))

def ue(u1,u2,u3):
    # u1=qb, u2=qe, u3=ub
    norm=normal(u1,u2)
    l=(2*(u3[0]*norm[0]+u3[1]*norm[1]))
    w=(norm[0]*l,norm[1]*l)
    return diff(w,u3)
#show( ue( qb(easy_walls[1]),qe(easy_walls[1]),ub(easy_walls[1])) )

def theta(u1,u2,u3,u4):
    # qb,qe,ub,ue
    if innprod(diff(u1,u2),u3)>=0:
        return arccos(innprod(u3,u4))
    if innprod(diff(u1,u2),u3)<0:
        return 2*pi-arccos(innprod(u3,u4))
#show(theta(qb(easy_walls[1]), qe(easy_walls[1]), ub(easy_walls[1]), ue(qb(easy_walls[1]),qe(easy_walls[1]),ub(easy_walls[1]))))

def center(u1,u2,u3):
    #qb,qe,ub
    cu=curv(u1,u2,u3)+eps
    return diff((-u3[1]/cu,u3[0]/cu),u1)
#show(center(qb(easy_walls[1]), qe(easy_walls[1]), ub(easy_walls[1])))

def pdistance(u1,u2):
    return ((u1[0]-u2[0])**2+(u1[1]-u2[1])**2)**(1/2)

def circ_or(v,w):
    if v[0]*w[1]-v[1]*w[0]>0:
        #print('counter clockwise')
        return 0
    else:
        #print('clockwise')
        return 1

def pointonarc(x,y,wall):
    u1=qb(wall)
    u2=qe(wall)
    u3=ub(wall)
    u4=ue(u1,u2,u3)
    C=center(u1,u2,u3)
    thet=theta(u1,u2,u3,u4)
    v1=diff(C,(x,y))
    v2=diff(C,u1)
    ori=circ_or(diff(C,u1),diff(C,u2))
    #print('ori',ori)
    ang=vecang(v1,v2)
    #print('theta', thet, 'ang',ang)
    if circ_or(v1,v2)==0:
        angdis=2*pi-ang
    else:
        angdis=ang
    if ori==1 and (thet>pi+eps or thet<pi-eps):
        if thet>2*pi-angdis:
            #print('1')
            return True
        if thet<2*pi-angdis:
            return False
    if ori==0 and (thet>pi+eps or thet<pi-eps):
        if thet>angdis:
            #print('2')
            return True
        if thet<angdis:
            return False
    if thet==pi:
        if innprod(v1,u3)>0:
            return True
        else:
            return False
        
def get_tangent(wall,x,y):
    if wall[6]==1:
        return getangle(diff(qb(wall),qe(wall)))
    if wall[6]==0:
        #print('circ tan')
        u1=qb(wall)
        u2=qe(wall)
        u3=ub(wall)
        C=center(u1,u2,u3)
        perp=getangle(diff(C,(x,y)))
        return perp+pi/2
    else:
        return False

def getoutnorm(vx,vy,ang):
    if vx*np.cos(ang)+vy*np.sin(ang)>0:
        return (np.cos(ang), np.sin(ang))
    else:
        return (-np.cos(ang), -np.sin(ang))
    

def veclen(u):
    return (u[0]**2+u[1]**2)**.5

def getangle(A):
    # return the angle of the vector A relative to the pos x-axis
    if A[0]>=0:
        if A[1]>=0:
            if (A[1]-eps)/veclen(A)>1:
                return arcsin(1)
            elif (A[1]-eps)/veclen(A)<-1:
                return arcsin(-1)
            else:
                return arcsin((A[1]-eps)/veclen(A))
        if A[1]<0:
            if (A[1]+eps)/veclen(A)>1:
                return arcsin(1)
            elif (A[1]+eps)/veclen(A)<-1:
                return arcsin(-1)
            else:
                return arcsin((A[1]+eps)/veclen(A))

    else:
        if (A[1]-eps)/veclen(A)>1:
            return pi-arcsin(1)
        elif (A[1]-eps)/veclen(A)<-1:
            return pi-arcsin(-1)
        else:
            return pi-arcsin((A[1]-eps)/veclen(A))

def wall_collide(pos,vel,wall):
    #print('in collide')
    # first check if the trajectory intersects the segment for linear boundaries
    # Note: works for only some circles--need better solution
    # It will excluse some legitimate circles if the trajectory passes in and out the arc
    if vel[0]==0:
            m2=10000000000000
    else:
            m2=vel[1]/vel[0]

    if wall[6]==1:
        v1=diff(pos,(wall[0],wall[1]))
        v2=diff(pos,(wall[2],wall[3]))
        #print(v1,v2,(vel[1],vel[2]))
        #print(vecang((vel[1],vel[2]),v1),vecang((vel[1],vel[2]),v2),vecang(v1,v2))
        #print(vecang((vel[1],vel[2]),v1)+vecang((vel[1],vel[2]),v2))
        angle_difference=vecang((vel[0],vel[1]),v1)+vecang((vel[0],vel[1]),v2)-vecang(v1,v2) #note cannot use vel directly as it is 3D
        #print(angle_difference)
        if angle_difference>0.000001:
            #print('does not collide')
            return (0,0,large+1)
        
        # find the intersection point
        if wall[2]==wall[0]:
            m1=(wall[3]-wall[1])/(eps)
        else:
            m1=(wall[3]-wall[1])/(wall[2]-wall[0])
        if m2==m1:
            m2=m2-eps
        x=(m2*pos[0]-pos[1]-m1*wall[0]+wall[1])/(m2-m1)
        y=pos[1]+m2*(x-pos[0])
        if ((x-pos[0])**2+(y-pos[1])**2)<eps:
            #print('not supposed to be here')
            return (0,0,large+1)
        else:
            #print('found a collision')
            return (x,y,((x-pos[0])**2+(y-pos[1])**2))
    # m2 is the slope of vel, handling infinite case by making it large
    
    if wall[6]==0:
        #print('checking arc')
        u1=qb(wall)
        u2=qe(wall)
        u3=ub(wall)
        C=center(u1,u2,u3)
        radi=pdistance(C,u1)
        #print('radius',radi)
        A=1+m2**2
        B=-2*C[0]+2*m2*(pos[1]-m2*pos[0]-C[1])
        c=C[0]**2+m2**2*pos[0]**2+2*m2*pos[0]*(C[1]-pos[1])+(pos[1]-C[1])**2-radi**2
        x1=(-B+(B**2-4*A*c)**(1/2))/(2*A)
        x2=(-B-(B**2-4*A*c)**(1/2))/(2*A)
        # This handles the cases where the trajectory misses the circle entirely
        if math.isnan(x1):
            return (0,0,large+1)
        #if x2.real==False:
        #    return (0,0,large+1)
        else:
            y1=pos[1]+m2*(x1-pos[0])
            y2=pos[1]+m2*(x2-pos[0])
            #print('x1,y1,pointonarc(x1,y1,wall)',x1,y1,pointonarc(x1,y1,wall))
            #print('x2,y2,pointonarc(x2,y2,wall)',x2,y2,pointonarc(x2,y2,wall))
            d2=pdistance(pos,(x2,y2))
            d1=pdistance(pos,(x1,y1))
            #print(d1,d2)
            pointon1=pointonarc(x1,y1,wall) and d1>eps
            pointon2=pointonarc(x2,y2,wall) and d2>eps
            if pointon1:
                #print('made it here')
                # to this point there is nothing to prevent the point being the original point
                # also, we need to make such the intersection is in forward time
                #print('here', d1>eps, ((vel[0]>0 and x1-pos[0]>0) or (vel[0]<0 and x1-pos[0]<0)) , d1<d2 or pointon2==False or  (vel[0]<0 and x2-pos[0]>0) or (vel[0]>0 and x2-pos[0]<0) )
                if d1>eps and ((vel[0]>0 and x1-pos[0]>0) or (vel[0]<0 and x1-pos[0]<0)) and (d1<d2 or pointon2==False or (vel[0]<0 and x2-pos[0]>0) or (vel[0]>0 and x2-pos[0]<0)):
                    #print('returning x1')
                    return (x1,y1,d1)
            #print(x2,y2,pointon2)
            if pointon2:
                #print('made it to x2 part')
                #print('d2>eps,d2,eps',d2>eps,d2,eps)
                if d2>eps and ((vel[0]>0 and x2-pos[0]>0) or (vel[0]<0 and x2-pos[0]<0)):
                    #print('returning x2')
                    return (x2,y2,d2)
            #print(x1,y1,d1)
            #print(x2,y2,d2)
        #print('emergency exit')
        return(0,0,large+1)

def new_positions(walls,rotational_position,x,y,rotational_velocity,vx,vy):
    #print('in next wall')
    best_dist=large
    bestx=x
    besty=y
    best_wall=0
    for i in range(0,len(walls)):
        #print('wall', i)
        #note: dist=distance squared
        (xnew,ynew,dist)=wall_collide((x,y),(vx,vy),walls[i])
        #print('dist',dist)
        if eps<dist<best_dist:
            bestx=xnew
            besty=ynew
            best_dist=dist
            best_wall=i
    return ((best_dist/(rotational_velocity**2+vx**2+vy**2)**(1/2))*rotational_velocity+rotational_position,bestx,besty,best_wall)

def reflect(wall,x,y,vr,vx,vy,eta):
    V=np.zeros(shape=(3,1),dtype=float)
    V[0,0]=vr
    V[1,0]=vx
    V[2,0]=vy
    #show(V)
    # omega is the angle of the line or tangent to the circle, oriented by the direction
    omega=get_tangent(wall,x,y)
    #print(norman)
    
    a=np.cos(np.pi*eta)
    b=np.sin(np.pi*eta)
    U=np.matrix([[-a,b,0],[b,a,0],[0,0,-1]],dtype=float) #rough collision with upward normal

    V=np.matmul(np.matmul(rot(omega),U),np.matmul(rot(-omega),V))
    cosphi= V[1,0]*np.cos(omega)+V[2,0]*np.sin(omega)
    return (V[0,0],V[1,0],V[2,0],cosphi)

# Print table
def draw_table(wall_segments):
    #print(len(wall_segments))
#     newgraph=Graphics()
    for i in range(0,len(wall_segments)):
        #print('i',i)
        if wall_segments[i][6]==1:
            # draw line
            plt.plot([wall_segments[i][0],wall_segments[i][2]],[wall_segments[i][1],wall_segments[i][3]],color='black')
        else:
            u1=qb(wall_segments[i])
            u2=qe(wall_segments[i])
            u3=ub(wall_segments[i])
            u4=ue(u1,u2,u3)
            #print(u1,u2,u3,u4)
            C=center(u1,u2,u3)
            #print('center',C)
            arc_angle=theta(u1,u2,u3,u4)
            #print('total angle', arc_angle)
            radius=1/abs(curv(u1,u2,u3))
            #print('radius', radius)
            x0=u1[0]
            y0=u1[1]
            startvec=diff(C,u1)
            #print('start vector',startvec)
            alpha=getangle(startvec)
            #print('alpha',alpha)
            for j in range(0,res+1):
                omega=(-1)**(circ_or(u3,diff(u1,u2)))*j/res*arc_angle
                #print(omega)
                #print(alpha)
                #print(radius*cos(alpha+omega),radius*sin(alpha+omega))
                #print(j,alpha+omega)
                x1=radius*cos(alpha+omega)+C[0]
                y1=radius*sin(alpha+omega)+C[1]
                #print(x1,y1)
                plt.plot([x0,x1], [y0,y1], color='black')
                x0=x1
                y0=y1
                #show(newgraph)

def BilliardMapIteration(wall_segments, state,eta):
    stoperror=False
    OP1 = state[1]
    OP2 = state[2]
    (state[0],state[1],state[2],Pwall)=new_positions(wall_segments,state[0],state[1],state[2],state[3],state[4],state[5])
    (state[3],state[4],state[5],cosphi)=reflect(wall_segments[Pwall],state[1],state[2],state[3],state[4],state[5],eta)
    #print(V0,V1,V2)
#     pathgraph=Draw_paths(OP1,OP2,P1,P2)

    if OP1==state[1] and OP2==state[2]:
        stoperror=True
    return (stoperror,cosphi)


#########################################################
# TABLE DEFINITION
#########################################################
# Tables are lists of walls oriented counterclockwise
# Each wall is a list
# [{0}starting x,
# {1} starting y,
# {2} ending x,
# {3} ending y,
# {4} direction vector x,
# {5} direction vector y,
# {6} shape: 0=arc 1=line,



#############################################################
# Automatic Table Generators
#############################################################
# curvilinear polygon: a regular n-gon with the sides arcs curving (in or out) according to crvang
# all walls no-slip
def make_crvngon(n,crvang):
    tab=[[1,0,cos(2*pi/n),sin(2*pi/n),cos(((2*pi+n*pi)/(2*n))+crvang),sin(((2*pi+n*pi)/(2*n))+crvang),0]]
    for i in range(1,n):
        tab+=[[cos(i*2*pi/n),sin(i*2*pi/n),cos((i+1)*2*pi/n),sin((i+1)*2*pi/n),cos(((2*pi+n*pi)/(2*n))+crvang+i*2*pi/n),sin(((2*pi+n*pi)/(2*n))+crvang+i*2*pi/n),0]]
    return tab

def make_ngon(n):
    tab=[[1,0,cos(2*pi/n),sin(2*pi/n),cos(2*pi/n)-1,sin(2*pi/n),1]]
    for i in range(1,n):
        tab+=[[cos(i*2*pi/n),sin(i*2*pi/n),cos((i+1)*2*pi/n),sin((i+1)*2*pi/n),cos((i+1)*2*pi/n)-cos(i*2*pi/n),sin((i+1)*2*pi/n)-sin(i*2*pi/n),1]]
    return tab

def make_strip(L):
    tab=[[-L,0,L,0,1,0,1]]
    tab+=[[L,0,L,1,0,1,1]]
    tab+=[[L,1,-L,1,-1,0,1]]
    tab+=[[-L,1,-L,0,0,-1,1]]
    return tab

def make_triangle(Angle, L):
    tab=[[0,0,L*cos((np.pi-Angle)/2),L*sin(np.pi-(np.pi-Angle)/2),cos((np.pi-Angle)/2),sin(np.pi-(np.pi-Angle)/2),1]]
    tab+=[[L*cos((np.pi-Angle)/2),L*sin(np.pi-(np.pi-Angle)/2),L*cos(np.pi-(np.pi-Angle)/2),L*sin(np.pi-(np.pi-Angle)/2),-1,0,1]]
    tab+=[[L*cos(np.pi-(np.pi-Angle)/2),L*sin(np.pi-(np.pi-Angle)/2),0,0,-cos(np.pi-(np.pi-Angle)/2),-sin(np.pi-(np.pi-Angle)/2),1]]
    return tab


######################################
# STATES DEFINITION
######################################
# States has N+1 (pairs of, for rolling) rows, each row being the state after collision m (before and after rolling)
# Each state has eleven components
#[{0} rotational position (seldom used),
# {1} x pos in fixed frame,
# {2} y pos in fixed frame,
# {3} rotational velocity (not used for specular collisions),
# {4} x velocity,
# {5} y velocity,
# {6} t=time,
# {7} s=arclength,
# {8} phi=angle relative to outward normal,
# {9} wall number]

#####################################
# Main Routines
#####################################
# These will vary a lot, but often we will want to run through a
# series of different tables in a family, looking at trajectories,
# phase portraits, or Lyapunov exponents
#####################################

#######################################################################
############################ PARAMETERS ###############################
#######################################################################
epsilon= 0.01

@functions_framework.http
def billiard(request):
    if request.method == "OPTIONS":
        # Allows GET requests from any origin with the Content-Type
        # header and caches preflight response for an 3600s
        headers = {
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET",
            "Access-Control-Allow-Headers": "Content-Type",
            "Access-Control-Max-Age": "3600",
        }

        return ("", 204, headers)

    # Set CORS headers for the main request
    headers = {"Access-Control-Allow-Origin": "*"}

    n = int(request.args.get('n', 1))
    x = float(request.args.get('x', 0))
    y = float(request.args.get('y', 0))
    curve = float(request.args.get('curve', 0))
    eta = float(request.args.get('eta', 0))
    sides = int(request.args.get('sides', 3))
    alpha = float(request.args.get('alpha', 0.1))
    spin = float(request.args.get('spin', 0))

    # plt.figure(figsize=(20,20))
    plt.axis('equal')
    #######################################################################
    #######################################################################

    isStrip = False
    isWedge = False
    if sides == 1:
        isStrip = True
    elif sides == 2:
        isWedge = True

    if(not (isStrip or isWedge)):
        plt.axis('equal')

    if isStrip == True:
        table = make_strip(1000)
    elif isWedge:
        table = make_triangle(curve, 1000)
    elif curve != 0:
        table=make_crvngon(sides,math.radians(curve))
    else:
        table = make_ngon(sides)

    vx0=np.cos(alpha)
    vy0=np.sin(alpha)
    norm = (vx0**2+vy0**2+spin**2)**0.5
    state=[0,x,y,spin/norm,vx0/norm,vy0/norm]

    miX = 1000000000
    miY = 1000000000
    maX = -1000000000
    maY = -1000000000

    for _ in range(0,n-1):
        x0 = state[1]
        y0 = state[2]
        miX = min(miX,x0)
        maX = max(maX,x0)
        maY = max(maY,y0)
        miY = min(miY,y0)
        (runerror,cosphi)=BilliardMapIteration(table, state,eta)
        x1 = state[1]
        y1 = state[2]
        miX = min(miX,x1)
        maX= max(maX,x1)
        maY = max(maY,y1)
        miY = min(miY, y1)
        plt.plot([x0,x1],[y0,y1],color='black')


    draw_table(table)
    if(isStrip):
        miY = 0
        maY = 1
    elif(isWedge):
        aMa = max(abs(miX),abs(maX))
        miX = -aMa
        maX = aMa
    scaleX = maX-miX
    scaleY = maY-miY
    if((isStrip == True or isWedge== True) and scaleX<200 and scaleY<200 and scaleX>=0 and scaleY >= 0):
        plt.xlim(miX-0.2*scaleX, maX+0.2*scaleX)
        plt.ylim(miY-0.2*scaleY,maY + 0.2*scaleY)
    miScale = 2*min(scaleX,scaleY)
    if(miScale == 0):
        miScale = max(scaleX,scaleY)
    if(miScale>1 or not (isStrip or isWedge) or miScale<=0):
        miScale = 1
    plt.arrow(x,y,miScale*0.07*np.cos(alpha),miScale*0.07*np.sin(alpha),width=miScale*0.01, color='red', head_width=miScale*0.05, head_length=miScale*0.05,zorder=10)

    bIO = BytesIO()
    plt.savefig(bIO,dpi=200)
    bIO.seek(0)
    plt.clf()

    return (send_file(bIO, mimetype='image/png'), 200, headers)