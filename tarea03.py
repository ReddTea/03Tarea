#!/usr/bin/env python
###################################
#Se importan librerías importantes#
###################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode             #para parte 2
from mpl_toolkits.mplot3d import Axes3D     #para parte 2


#Algunos Parametros de la parte 1
mu=1.067 #RUT=19.077.067-2

#parámetros de la parte 2
sigma=10.0
beta=8.0/3.0
rho=28.0

#condiciones iniciales parte 1
y0_1=0.1
v0_1=0
y0_2=4
v0_2=0
n_steps=40000 #esto lo dejaremos acá
h=20*np.pi/n_steps #ancho del paso
tiempo=np.linspace(0,20*np.pi,n_steps) #tiempo

#condiciones iniciales de la parte 2
n_steps2 = 10000
tiempo2=np.linspace(0,100,n_steps2) #tiempo

###########################################
#Se definen funciones para resolver la edo#
###########################################

def vdp(y, v, mu=mu):
    '''
    Se le entrega y,v y el factor mu, nos devuelve la edo como arreglo
    '''
    return v, -y-mu*((y**2)-1)*v

def get_k1(y_n, v_n, h, vdp):
    '''
    Se le entrega y_n, v_n y devuelve el vector k1
    '''
    vdp_n = vdp(y_n, v_n)
    return h*vdp_n[0], h*vdp_n[1]

def get_k2(y_n, v_n, h, vdp):
    '''
    Se le entrega y_n, v_n y devuelve el vector k2
    '''
    k1 = get_k1(y_n, v_n, h, vdp)
    vdp_n= vdp(y_n+k1[0]/2., v_n+k1[1]/2.)
    return h*vdp_n[0], h*vdp_n[1]

def rk3_step(y_n, v_n, h, f):
    '''
    Vital, le entregamos y_n y nos devuelve el y_{n+1}
    '''
    k1 = get_k1(y_n, v_n, h, f)
    k2 = get_k2(y_n, v_n, h, f)
    k3 = get_k2(y_n, v_n, h, f)

    y_n1 = y_n+(k1[0]+4*k2[0]+k3[0])/6.0
    v_n1 = v_n+(k1[1]+4*k2[1]+k3[1])/6.0
    return y_n1, v_n1

def ecuacionlorentz(t, r):  #para la parte 2
    '''
    se le entrega el tiempo como arreglo (t) y los valores de x, y ,z
    en el vector desplazamiento r.
    Devuelve rr que son los valores de dx, dy, dz
    '''
    dx= (r[1]- r[0])*sigma          #dx
    dy= r[0]*(rho- r[2])- r[1]      #dy
    dz= r[0]*r[1]- beta*r[2]       #dz
    return [dx,dy,dz]

'''
######################################
#se setean las condiciones iniciales #
######################################
y=[y0_1]
v=[v0_1]

######################################
#Se echa a correr el RugneKutta      #
######################################

for i in range(n_steps-1):
    sig= rk3_step(y[i-1], v[i-1], h, vdp)
    y=np.append(y,sig[0])
    v= np.append(v,sig[1])

plt.figure(1)
plt.plot(tiempo,y)
plt.xlabel('$t[s]$', fontsize=20)
plt.ylabel('$y(t)$', fontsize=20)
plt.title('$y(t)$ vs t (con $\mu^*$= 1.067, $y(0)=0.1$, $v(0)=0$)', fontsize=24)
plt.savefig('y_vs_t_01.png')

plt.figure(2)
plt.plot(y,v)
plt.xlabel('$y(t)$', fontsize=20)
plt.ylabel('$v(t)$', fontsize=20)
plt.title('$y(t)$ vs v(t) (con $\mu^*= 1.067$, $y(0)=0.1$, $v(0)=0$)', fontsize=24)
plt.savefig('y_vs_v_01.png')

########################################
#se setean las condiciones iniciales 2 #
########################################
y=[y0_2]
v=[v0_2]

for i in range(n_steps-1):
    sig= rk3_step(y[i-1], v[i-1], h, vdp)
    y=np.append(y,sig[0])
    v= np.append(v,sig[1])

plt.figure(3)
plt.plot(tiempo,y)
plt.xlabel('$t[s]$', fontsize=20)
plt.ylabel('$y(t)$', fontsize=20)
plt.title('$y(t)$ vs t (con $\mu^*$= 1.067, $y(0)=4.0$, $v(0)=0$)', fontsize=24)
plt.savefig('y_vs_t_40.png')

plt.figure(4)
plt.plot(y,v)
plt.xlabel('$y(t)$', fontsize=20)
plt.ylabel('$v(t)$', fontsize=20)
plt.title('$y(t)$ vs v(t) (con $\mu^*= 1.067$, $y(0)=4.0$, $v(0)=0$)', fontsize=24)
plt.savefig('y_vs_v_40.png')

plt.show()
plt.draw()
'''
'''
convertir esto en funcion si queda tiempo ^
'''
#######################################
################PARTE 2################
#######################################

r0= [1, 2, 1] #condicion inicial
#donde meteré los x,y,z
x=[]
y=[]
z=[]

#se setea el ode solver
marcapagina= ode(ecuacionlorentz).set_integrator('dopri5') #se usa RK4

marcapagina.set_initial_value(r0)
#se echa a correr
while marcapagina.successful() and marcapagina.t<100: #100=tf
    marcapagina.integrate(marcapagina.t+0.01) #0.01=dt
    x=np.append(x,marcapagina.y[0])
    y=np.append(y,marcapagina.y[1])
    z=np.append(z,marcapagina.y[2])

#se plotea
fig=plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('auto')
ax.plot(x, y, z, 'g')
plt.title('Atractor de Lorenz')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.grid()
plt.show()
plt.draw()
