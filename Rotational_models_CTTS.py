import concurrent.futures
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d
import random
from scipy import optimize

def Run(Mass, Prot, Macc, Tdisk, Bfield, betta, gamma, APSW):

  #####################################################
  Msun=1.98e+33 #g
  Rsun=6.96e+10 #cm
  Omegasun=2.8e-06 #s-1
  Isun=6.4e+53 #g*cm2
  Jsun=1.8e+48 #g*cm2/s
  Kw=9e48#2.7e47#*3.154e7*3.154e7 #g*cm2/yr2
  tg=35e+6 #Escala de tiempo gravitacional en yr
  ta=2.1e+6 #2.3e6 yr fedele 2010
  G=6.67259e-8 #cgs
  yr2sec=3.154e+7 #s
  #####################################################

  #Lectura del archivo de baraffe

  baraffe=np.genfromtxt('./Baraffe/'+'Baraffe_'+str(Mass)+'Msun.dat', comments='#', dtype='S')

  #Datos del modelo solar de Baraffe 2015

  RADIUS=baraffe[:,5].astype(float)
  AGE=pow(10,baraffe[:,1].astype(float))
  TEFF=baraffe[:,2].astype(float)
  LUM=baraffe[:,3].astype(float)
  K2ENV=baraffe[:,11].astype(float)
  K2CORE=baraffe[:,12].astype(float)

  R=interp1d(AGE/tg,RADIUS,kind='linear')
  C=interp1d(AGE/tg,K2ENV,kind='linear')
  D=interp1d(AGE/tg,K2CORE,kind='linear')

  S=interp1d(AGE/tg,LUM,kind='linear')
  T=interp1d(AGE/tg,TEFF,kind='linear')

  ##############################################################################################

  #Tiempo normalizado a la escala de tiempo gravitacional tg
  time=np.arange(5e5,Tdisk,1e5)/tg #5e9

  radius=R(time)
  K2env=C(time)
  K2core=D(time)
  K2total=(K2env*K2env)+(K2core*K2core)

  #Tasa de acrecion
  Madot=(Macc)*np.exp(-time*tg/ta)

  #Masa estelar (Solucion Analitica)
  Mstar=-(ta)*(Macc)*np.exp(-time*tg/ta)*(1-APSW)+ Mass


  #Momento de inercia de la estrella
  Istar=interp1d(time,K2total*radius*radius*Mstar, kind='linear')
  I_total=Istar(time)

  #Derivada del Radio Estelar (En desuso)

  def R_dot(radius,time):
    "Derivada del radio estelar"
    rdot=[]
    for i in range(len(time)):
      if (i==0):
        rdot.append((radius[i+1]-radius[i])/(time[i+1]-time[i]))
      if (i==len(time)-1):
        return np.array(rdot)
      else:
        rdot.append((radius[i+1]-radius[i])/(time[i+1]-time[i]))

  #Interpolacion de la derivada del radio
  Rdot=interp1d(time,R_dot(radius,time),kind='linear')

  #Derivada del momento de inercia

  def Idot(I_total,time):
    "Derivada del momento de inercia"
    I_dot=[]
    for i in range(len(time)):
      if (i==0):
        I_dot.append((I_total[i+1]-I_total[i])/(time[i+1]-time[i]))
      if (i==len(time)-1):
        return np.array(I_dot)
      else:
        I_dot.append((I_total[i+1]-I_total[i])/(time[i+1]-time[i]))


  #Interpolacion de la derivada del momento de inercia

  J=interp1d(time,Idot(I_total,time),kind='linear')
  I_dot=J(time)

  ##################################################################################################

  #Ecuacion diferencial a resolver

  def Omega_dot(Omega, time, tdisk, B, betta, gamma, APSW):

    #Tasa de acrecion
    Madot=(Macc)*np.exp(-time*tg/ta)

    #Masa de la estrella
    Mstar=-(ta)*(Macc)*np.exp(-time*tg/ta)*(1-APSW)+ Mass

    #Fase en donde existe un disco protoplanetario
    if (time<tdisk):
      #Velocidad angular normalizada a la velocidad de ruptura
      f=Omega*Omegasun*np.sqrt(pow(R(time)*Rsun,3.0)/(G*Mstar*Msun))

      #Limite para la velocidad angular normalizada
      if (f>1):
        f=1
        Omega=f/(Omegasun*np.sqrt(pow(R(time)*Rsun,3.0)/(G*Mstar*Msun)))

      #Radio de Corrotacion
      Rco=pow((G*Mstar*Msun/(pow(Omega*Omegasun,2.0))),1.0/3.0)/Rsun

      #Radio Externo a la region de corrotacion
      Rout=Rco*pow(1+betta*gamma, 2.0/3.0)

      #Radio Interno a la region de corrotacion
      Rin=Rco*pow(1-betta*gamma, 2.0/3.0)

      #Dipolo magnetico
      dipole=B*pow(R(time)*Rsun,3.0)

      psi=2*pow(dipole,2.0)*pow(Madot*Msun/yr2sec,-1.0)*pow(np.sqrt(G*Mstar*Msun),-1.0)*pow(R(time)*Rsun,-3.5)

      #Ecuacion 15 (MP05) a resolver a traves del metodo Newton Raphson (Proposito: Encontrar el valor de x=Rt/Rco)
      y=[lambda x: pow(x, -3.5) - pow(x, -2.0) - (betta*pow(f,-2.33)/psi), lambda x: 2.0*pow(x,-3.0) - 3.5 * pow(x,-4.5), lambda x: -6*pow(x,-4.0) + 15.75 * pow(x,-5.5)]

      #Condicion especifica para el campo magnetico
      if (B==0):
        Rt=R(time)
        torquemag=0

      # (Estado 1) Lineas de campo magnetico abiertas mas alla del radio de truncamiento
      if (f<(1-betta*gamma)*pow(gamma*psi,(-3.0/7.0))):
        Rt = pow(gamma*psi, 2.0/7.0)*R(time)
        if (Rt<R(time)):
          Rt=R(time)
        torquemag=0

      #(Estado 2) Disk-Locking Phase
      else:
        zero=optimize.newton(y[0],0.5, fprime=y[1], fprime2=y[2]) #zero is Rt/Rco
        if (zero<R(time)/Rco):
          zero=R(time)/Rco
        Rt=zero*Rco
        torquemag=(pow(dipole,2.0)/(3.0*betta*pow(Rco*Rsun,3.0))) * (2.0*pow(Rco/Rout,1.5)-pow(Rco/Rout,3.0)-2.0*pow(Rco/Rt,1.5)+pow(Rco/Rt,3.0))

      torqueacc=(Madot*Msun/yr2sec)*np.sqrt(G*Mstar*Msun*Rt*Rsun)

      rA=2.11*pow(B*B*R(time)*Rsun*R(time)*Rsun/(APSW*(Madot*Msun/(yr2sec))*np.sqrt(2*G*Mstar*Msun/(R(time)*Rsun))),0.223)

      torquewind=-APSW*(Madot*Msun/(yr2sec))*Omega*Omegasun*rA*rA*R(time)*R(time)*Rsun*Rsun

      torque=(torqueacc+torquemag+torquewind)

      ode=((tg*yr2sec)*torque/(Istar(time)*Msun*Rsun*Rsun*Omegasun))-(Omega*(J(time)/Istar(time)))

      return ode, torqueacc, torquemag, Rt, Rco, R(time), Mstar, f, torquewind

  ##################################################################################################

  def solution(omegainicial, tiempoinicial, tdisk, B, betta, gamma, APSW):
	  # Set initial conditions.
	  t = tiempoinicial
	  x = omegainicial

	  # Set initial step size.
	  dt = 1e4/tg

	  # Set minimal step size.
	  dt_min = 1e3/tg

	  # Set relative change tolerances.
	  dx_max = 1e-1 # Enables faster speed.
	  dx_min = 1e-3 # Controls accuracy.
	  x_tol = 1e-3

	  a=[]
	  b=[]

	  tormag=[]
	  toracc=[]
	  torwind=[]
	  Rcorr=[]
	  Rtr=[]
	  Rast=[]
	  Mast=[]
	  timer=[]
	  FF=[]

	  a.append(x)
	  b.append(t)


	  while (t < max(time)-dt*2):
		  #print("age=%i, Omega=%.3f\n" %(t*tg, x))
		  # Calculate partial steps.
		  k1 = Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[0]
		  k2 = Omega_dot(x+dt*k1/2, t+dt/2, tdisk, B, betta, gamma, APSW)[0]
		  k3 = Omega_dot(x+dt*k2/2, t+dt/2, tdisk, B, betta, gamma, APSW)[0]
		  k4 = Omega_dot(x+dt*k3, t+dt, tdisk, B, betta, gamma, APSW)[0]
		  # Combine partial steps.
		  step_x = x + dt/6*(k1+2*k2+2*k3+k4)

		  # Calculate partial steps.
		  k2 = Omega_dot(x+dt*k1/4, t+dt/4, tdisk, B, betta, gamma, APSW)[0]
		  k3 = Omega_dot(x+dt*k2/4, t+dt/4, tdisk, B, betta, gamma, APSW)[0]
		  k4 = Omega_dot(x+dt*k3/2, t+dt/2, tdisk, B, betta, gamma, APSW)[0]
		  # Combine partial steps.
		  half_step_x = x + dt/12*(k1+2*k2+2*k3+k4)

		  # Calculate partial steps.
		  k2 = Omega_dot(x+dt*k1, t+dt, tdisk, B, betta, gamma, APSW)[0]
		  k3 = Omega_dot(x+dt*k2, t+dt, tdisk, B, betta, gamma, APSW)[0]
		  k4 = Omega_dot(x+2*dt*k3, t+2*dt, tdisk, B, betta, gamma, APSW)[0]
		  # Combine partial steps.
		  dble_step_x = x + dt/3*(k1+2*k2+2*k3+k4)

		  if (abs(step_x) < x_tol): # Use a fixed step size for small values of x.
			  if (dt != dt_min):
				  #print("New step size",dt_min)
				  dt = dt_min
			  new_x = step_x
		  else:
			  if (abs(step_x) > x_tol and abs(step_x-half_step_x)/abs(step_x) > dx_max):
				  dt = dt/2 # Error is too large; decrease step size.
				  #print("New step size",dt)
				  new_x = half_step_x
			  elif (abs(step_x) > x_tol and abs(step_x-dble_step_x)/abs(step_x) < dx_min):
				  dt = dt*2 # Larger error is acceptable; increase step size.
				  #print("New step size",dt)
				  new_x = dble_step_x
			  else:
				  new_x = step_x # This step size is just right.

		  x = new_x
		  a.append(x)
		  t = t + dt
		  b.append(t)

		  tormag.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[2])
		  toracc.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[1])
		  Rcorr.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[4])
		  Rtr.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[3])
		  Rast.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[5])
		  Mast.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[6])
		  FF.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[7])
		  torwind.append(Omega_dot(x, t, tdisk, B, betta, gamma, APSW)[8])
		  timer.append(t)
	  return np.array(a), np.array(b), np.array(toracc), np.array(tormag), np.array(Rcorr), np.array(Rtr), np.array(Rast), np.array(Mast), np.array(FF), np.array(timer), np.array(torwind)


  ###################################################################################################
  tdisk=Tdisk/tg
  Omegasat=10
  Kconstant=2.5e-4

  omega=2.0*np.pi/(Prot*Omegasun*86400)

  Omega_star, times, torA, torM, RCO, RTR, RAST, MAST, ff, TIMER, torW = solution(omega, time[0], tdisk, Bfield, betta, gamma, APSW) #Omega,time0,disktime,Bfield,betta,gamma

  period=2.0*np.pi/(Omega_star[:-1]*Omegasun*86400)
  vrot=Omega_star[:-1]*Omegasun*RAST*Rsun/1e5

  #################################################################################################

  return TIMER*tg, vrot, period
