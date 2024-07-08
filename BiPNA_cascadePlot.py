# now plot the second half as well:
# I want code to copy and print a number of profiles

import math
import matplotlib.pyplot as plt
import numpy as np

cos = math.cos
pi = math.pi
sin = math.sin
tan = math.tan


#==============================================
# number of elements to copy
no_cp = 5

no_cp2 = 3

no_cp3 = 5

# actuator length
a_len = 5
#actuator angle
a_ang = 10/2
#actuator height
a_hei = 10


a_len2 = 10#5
#actuator angle
a_ang2 = 24#10/2
#actuator height
a_hei2 = 10#10

# actuator length
a_len3 = 5
#actuator angle
a_ang3 = 10/2
#actuator height
a_hei3 = 10

#starting point
a_x = 0
a_y = 0
#==============================================

# Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to matrix multiply
def Mx(X, Y):
    T = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*Y)] for X_row in X]
    return T
#function to rotate
def Rx(X,deg_a):
    rot = -deg_a*2*pi/360
    #T = [[1,0,0,0],[0,cos(rot),-sin(rot),0],[0,sin(rot),cos(rot),0],[0,0,0,1]]
    #T = [[cos(rot),0,sin(rot),0],[0,1,0,0],[-sin(rot),0,cos(rot),0],[0,0,0,1]]
    T = [[round(cos(rot),7),round(-sin(rot),7),0,0],[round(sin(rot),7),round(cos(rot),7),0,0],[0,0,1,0],[0,0,0,1]]
    T2 = Mx(T,X)
    return T2
# function to translate
def Tx(X,len_a):
    tra_x = len_a
    tra_y = 0
    tra_z = 0
    #[1 0 0 Txyz(1); 0 1 0 Txyz(2); 0 0 1 Txyz(3); 0 0 0 1]
    T = [[1,0,0,tra_x],[0,1,0,tra_y],[0,0,1,tra_z],[0,0,0,1]]
    T2 = Mx(T,X)
    return T2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


p1x = a_x-a_len/2- round(a_hei/2*math.tan(a_ang*  (2*math.pi/360)),7)
p1y = a_hei/2

p2x = a_x+a_len/2+ round(a_hei/2*math.tan(a_ang*  (2*math.pi/360)),7)
p2y = a_hei/2

p3x = -a_len/2
p3y = 0

p4x = +a_len/2
p4y = 0

p5x = a_x-a_len/2+ round(a_hei/2*math.tan(a_ang*  (2*math.pi/360)),7)
p5y = -a_hei/2

p6x = a_x+a_len/2- round(a_hei/2*math.tan(a_ang*  (2*math.pi/360)),7)
p6y = -a_hei/2

Red_m = [[p1x,p2x,p3x,p4x,p5x,p6x],[p1y,p2y,p3y,p4y,p5y,p6y],[0,0,0,0,0,0],[1,1,1,1,1,1]]

# I need a reference matrix to transform cascades
Id = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
Id0 =[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

Tp = Rx(Tx(Red_m,a_len/2),a_ang)

xpoints = np.array(Tp[0])
ypoints = np.array(Tp[1])

xp_box = np.array([Tp[0][0],Tp[0][1],Tp[0][3],Tp[0][2],Tp[0][0],Tp[0][4],Tp[0][5],Tp[0][3]])
yp_box = np.array([Tp[1][0],Tp[1][1],Tp[1][3],Tp[1][2],Tp[1][0],Tp[1][4],Tp[1][5],Tp[1][3]])

for i in range(no_cp):

    
    plt.plot([0],[0],'o')
    plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'red')
    plt.plot(xpoints, ypoints, 'o',color='darkblue')
    
    if i == 0:
        #Tp = Rx(Tx(Red_m,a_len/2),a_ang)
        #Id = Rx(Tx(Id,a_len/2),a_ang)
        Tp = Mx(Id,Rx(Tx(Red_m,a_len/2),a_ang))
        Id = Mx(Id,Rx(Tx(Id0,a_len/2),a_ang))
    else:
        #Tp = Tx(Rx(Tx(Tp,a_len/2),a_ang*2),a_len/2)
        #Id = Tx(Rx(Tx(Id,a_len/2),a_ang*2),a_len/2)
        Tp = Mx(Id,Tx(Rx(Tx(Red_m,a_len/2),a_ang*2),a_len/2))
        Id = Mx(Id,Tx(Rx(Tx(Id0,a_len/2),a_ang*2),a_len/2))
    
    xpoints = np.array(Tp[0])
    ypoints = np.array(Tp[1])

    xp_box = np.array([Tp[0][0],Tp[0][1],Tp[0][3],Tp[0][2],Tp[0][0],Tp[0][4],Tp[0][5],Tp[0][3]])
    yp_box = np.array([Tp[1][0],Tp[1][1],Tp[1][3],Tp[1][2],Tp[1][0],Tp[1][4],Tp[1][5],Tp[1][3]])


plt.plot([0],[0],'o')
plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'red')
plt.plot(xpoints, ypoints, 'o',color='darkblue')        
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.plot([0],[0],'o')    


# plotting second half

p1x2 = a_x-a_len2/2- round(a_hei2/2*math.tan(a_ang2*  (2*math.pi/360)),7)
p1y2 = a_hei2/2

p2x2 = a_x+a_len2/2+ round(a_hei2/2*math.tan(a_ang2*  (2*math.pi/360)),7)
p2y2 = a_hei2/2

p3x2 = -a_len2/2
p3y2 = 0

p4x2 = +a_len2/2
p4y2 = 0

p5x2 = a_x-a_len2/2+ round(a_hei2/2*math.tan(a_ang2*  (2*math.pi/360)),7)
p5y2 = -a_hei2/2

p6x2 = a_x+a_len2/2- round(a_hei2/2*math.tan(a_ang2*  (2*math.pi/360)),7)
p6y2 = -a_hei2/2

Red_m = [[p1x2,p2x2,p3x2,p4x2,p5x2,p6x2],[p1y2,p2y2,p3y2,p4y2,p5y2,p6y2],[0,0,0,0,0,0],[1,1,1,1,1,1]]

Tp2 = Mx(Id,Red_m)

for i in range(no_cp2):
    
    print('i:',i)
    if i == 0:
        #Tp2 = Mx(Id,Tx(Rx(Tx(Red_m,a_len2/2),a_ang2+a_ang),a_len/2))
        #Id = Mx(Id,Tx(Rx(Tx(Id0,a_len2/2),a_ang2+a_ang),a_len/2))
        Tp2 = Mx(Id,Tx(Rx(Tx(Red_m,a_len2/2),a_ang2+a_ang),a_len/2))
        Id = Mx(Id,Tx(Rx(Tx(Id0,a_len2/2),a_ang2+a_ang),a_len/2))
    
    else:  
        #Tp2 = Mx(Id,Tx(Red_m,a_len2/2))
        Tp2 = Mx(Id,Tx(Rx(Tx(Red_m,a_len2/2),a_ang2*2),a_len2/2))
        Id = Mx(Id,Tx(Rx(Tx(Id0,a_len2/2),a_ang2*2),a_len2/2))
       
    xpoints = np.array(Tp2[0])
    ypoints = np.array(Tp2[1])

    xp_box = np.array([Tp2[0][0],Tp2[0][1],Tp2[0][3],Tp2[0][2],Tp2[0][0],Tp2[0][4],Tp2[0][5],Tp2[0][3]])
    yp_box = np.array([Tp2[1][0],Tp2[1][1],Tp2[1][3],Tp2[1][2],Tp2[1][0],Tp2[1][4],Tp2[1][5],Tp2[1][3]])


    plt.plot([0],[0],'o')
    plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'green')
    plt.plot(xpoints, ypoints, 'o',color='darkred')


plt.plot([0],[0],'o')
plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'green')
plt.plot(xpoints, ypoints, 'o',color='darkred')        
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.plot([0],[0],'o')    

#  Third section
p1x3 = a_x-a_len3/2- round(a_hei3/2*math.tan(a_ang3*  (2*math.pi/360)),7)
p1y3 = a_hei3/2

p2x3 = a_x+a_len3/2+ round(a_hei3/2*math.tan(a_ang3*  (2*math.pi/360)),7)
p2y3 = a_hei3/2

p3x3 = -a_len3/2
p3y3 = 0

p4x3 = +a_len3/2
p4y3 = 0

p5x3 = a_x-a_len3/2+ round(a_hei3/2*math.tan(a_ang3*  (2*math.pi/360)),7)
p5y3 = -a_hei3/2

p6x3 = a_x+a_len3/2- round(a_hei3/2*math.tan(a_ang3*  (2*math.pi/360)),7)
p6y3 = -a_hei3/2

Red_m = [[p1x3,p2x3,p3x3,p4x3,p5x3,p6x3],[p1y3,p2y3,p3y3,p4y3,p5y3,p6y3],[0,0,0,0,0,0],[1,1,1,1,1,1]]

Tp3 = Mx(Id,Red_m)

for i in range(no_cp2):
    
    print('i:',i)
    if i == 0:
        #Tp3 = Mx(Id,Tx(Rx(Tx(Red_m,a_len3/2),a_ang3+a_ang),a_len/2))
        #Id = Mx(Id,Tx(Rx(Tx(Id0,a_len3/2),a_ang3+a_ang),a_len/2))
        Tp3 = Mx(Id,Tx(Rx(Tx(Red_m,a_len3/2),a_ang3+a_ang2),a_len2/2))
        Id = Mx(Id,Tx(Rx(Tx(Id0,a_len3/2),a_ang3+a_ang2),a_len2/2))
    
    else:  
        #Tp3 = Mx(Id,Tx(Red_m,a_len3/2))
        Tp3 = Mx(Id,Tx(Rx(Tx(Red_m,a_len3/2),a_ang3*2),a_len3/2))
        Id = Mx(Id,Tx(Rx(Tx(Id0,a_len3/2),a_ang3*2),a_len3/2))
       
    xpoints = np.array(Tp3[0])
    ypoints = np.array(Tp3[1])

    xp_box = np.array([Tp3[0][0],Tp3[0][1],Tp3[0][3],Tp3[0][2],Tp3[0][0],Tp3[0][4],Tp3[0][5],Tp3[0][3]])
    yp_box = np.array([Tp3[1][0],Tp3[1][1],Tp3[1][3],Tp3[1][2],Tp3[1][0],Tp3[1][4],Tp3[1][5],Tp3[1][3]])


    plt.plot([0],[0],'o')
    plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'green')
    plt.plot(xpoints, ypoints, 'o',color='darkred')


plt.plot([0],[0],'o')
plt.plot(xp_box, yp_box, linestyle = 'solid',color = 'green')
plt.plot(xpoints, ypoints, 'o',color='darkred')        
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.plot([0],[0],'o')  
    
plt.show()


