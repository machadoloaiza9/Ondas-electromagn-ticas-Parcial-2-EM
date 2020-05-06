#Importamos librerias necesarias para el código

import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import math

#Discretizamos el espacio 
Npts = 500

#Nombramos las variables de la permitividad elÃ©ctrica y magnÃ©tica en el vacÃ­o  
epsi_0 = constants.epsilon_0
mu_0 = constants.mu_0
velocidadpropagacion = constants.c #Propagación en el vacio

#Campo eléctrico, campo magnético y vector asociado al cambio medio
Ez = np.zeros((Npts), dtype=float)
Hy = np.zeros((Npts), dtype=float)
eps = np.zeros(Npts)

#Inicializamos el campo eléctrico en cero y frecuencia entre 10-10.5GHz
E0= 0
Hertz= 10

#Dimensiones físicas del problema, dimesiones para el dieléctrico
bordeinferior_dielectrico = 200
bordesuperior_dielectrico = 400

##Condiciones de frontera en los bordes, condicion Dirichlet
limite_der = [0, 0]
limite_izq = [0, 0]

#Nuevo valor para permitividad elÃ©ctrica en el dielecticro
epsilon_dielectrico = 12

# VariaciÃ³n de la permitividad por influencia del dielÃ©ctrico
eps[0:bordeinferior_dielectrico] = epsi_0
eps[bordeinferior_dielectrico:bordesuperior_dielectrico] = epsilon_dielectrico
eps[bordesuperior_dielectrico:Npts] = epsi_0

# Condiciones de malla
dt = 0.018e-9                          #Diferencial de tiempo
dx = 2.1*dt*velocidadpropagacion   #Diferencial de espacio basado en teorema de 
                                    #Nyquist 

# Impedancias intrinsecasde cada medio
n1 = (mu_0/epsi_0)**(1/2)
n2 = (mu_0/(epsilon_dielectrico*epsi_0))**(1/2)

# Coeficiente de transmision
t12 =(2*n1)/(n1+n2) 
t21 =(2*n2)/(n1+n2)

#Coeficientes de reflexion
r12 =(n2-n1)/(n1+n2)
r21 =(n1-n2)/(n1+n2)

input('Ingrese el número de pasos para el programa')
Numero_pasos = int(input())


##############################################################################
#Utilizacion FDTD Method, tomado del libro ELECTROMAGNETIC SIMULATION USING 
#THE FDTD METHOD WITH PYTHON, Sullivan

#ex [k] = ex [k] + 0.5 ∗ (hy [k − 1] − hy [k]) 
#hy [k] = hy [k] + 0.5 ∗ (ex [k] − ex [k + 1] 
##############################################################################


for i in range(0, Numero_pasos):     #Update the electric field
    
    for j in range(1, Npts):
        Ez[j] = Ez[j] + (dt/(eps[j]*dx))*(Hy[j] - Hy[j-1])
    
    # Pulso senoidal
    if i < 200:
        pulse = math.sin(2 * math.pi * Hertz * dt * i)
    Ez[0] = pulse
    
    if i > 201:
        Ez[0] = limite_izq.pop(0)
        limite_izq.append(Ez[1])
        Ez[Npts - 1] = limite_der.pop(0)
        limite_der.append(Ez[Npts - 2])
     
    for k in range(Npts - 1):   #Update the magnetic field
        Hy[k] = Hy[k] + (dt/(mu_0*dx))*(Ez[k+1] - Ez[k])

#Grafica  
plt.plot(Ez)
plt.plot([bordeinferior_dielectrico, bordeinferior_dielectrico],[-(E0+0.6), E0+0.6])
plt.plot([bordesuperior_dielectrico, bordesuperior_dielectrico],[-(E0+0.6), E0+0.6])
plt.ylabel('Campo eléctrico [V/m]')
plt.xlabel('Desplazamiento')
plt.show()