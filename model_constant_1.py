import concurrent.futures
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d
import random
from scipy import optimize
import smtplib, ssl
from email.message import EmailMessage
from scipy.stats import ks_2samp
from scipy.stats import norm

def email_sender():
    port = 465  # For SSL
    smtp_server = "smtp.gmail.com"
    sender_email = "tessextractor.app@gmail.com"  # Enter your address
    receiver_email = "jserna@astro.unam.mx"
    password = "vdvkjdfkfuqyqwvs"
    subject = f'Confirmation Job Done'
    body = """
    Job done!
    """
    em = EmailMessage()
    em['From'] = sender_email
    em['To'] = receiver_email
    em['Subject'] = subject
    em.set_content(body)
    context = ssl.create_default_context()
    with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, em.as_string())


def REFUGEE(Mass, Prot, Macc, Tdisk, Bfield, betta, gamma, APSW):

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

  #Masa del disco en unidades solares
  #Mdisk=0.01
  #Tasa de acrecion
  Madot=(Macc)*np.exp(-time*tg/ta)

  #Masa estelar (Solucion Analitica)
  #Mstar=-(ta)*(Macc)*np.exp(-time*tg/ta)*(1-APSW)+1
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

    #Campo Magnetico
    #B=2500 #Unidad de Gauss
    #global betta, gamma
    #betta=0.01
    #gamma=1

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
        #print("state1:" "Rco=%f,Rt=%f,Rin=%f,Rout=%f,f=%f\n" %(Rco, Rt, Rin, Rout, f))
        if (Rt<R(time)):
          Rt=R(time)
        torquemag=0

      #(Estado 2) Disk-Locking Phase
      else:
        zero=optimize.newton(y[0],0.5, fprime=y[1], fprime2=y[2]) #zero is Rt/Rco
        if (zero<R(time)/Rco):
          zero=R(time)/Rco
        Rt=zero*Rco
        #print("state2:" "Rco=%f,Rt=%f,Rin=%f,Rout=%f,f=%f\n" %(Rco, Rt, Rin, Rout, f))
        #torquemag=(pow(dipole,2.0)/(3.0*betta*pow(Rco*Rsun,3.0))) * (2.0*pow(1+betta*gamma,-1.0)-pow(1+betta*gamma,-2.0)-2.0*pow(Rco/Rt,1.5)+pow(Rco/Rt,3.0))
        torquemag=(pow(dipole,2.0)/(3.0*betta*pow(Rco*Rsun,3.0))) * (2.0*pow(Rco/Rout,1.5)-pow(Rco/Rout,3.0)-2.0*pow(Rco/Rt,1.5)+pow(Rco/Rt,3.0))

      #print("Rco=%f,Rt=%f,R=%f,f=%f,Mstar=%f\n" %(Rco, Rt, R(time), f, Mstar))

      torqueacc=(Madot*Msun/yr2sec)*np.sqrt(G*Mstar*Msun*Rt*Rsun)

      rA=2.11*pow(B*B*R(time)*Rsun*R(time)*Rsun/(APSW*(Madot*Msun/(yr2sec))*np.sqrt(2*G*Mstar*Msun/(R(time)*Rsun))),0.223)

      torquewind=-APSW*(Madot*Msun/(yr2sec))*Omega*Omegasun*rA*rA*R(time)*R(time)*Rsun*Rsun

      torque=(torqueacc+torquemag+torquewind)

      #print(Omega*(J(time)/Istar(time)),((tg*yr2sec)*torque/(Istar(time)*Msun*Rsun*Rsun*Omegasun)))

      #print(torqueacc,torquemag)

      #print("state=%i, f=%f, psi=%f, cond=%f\n" %(state,f,psi,(1-betta)*pow(psi,-3.0/7.0)))

      ode=((tg*yr2sec)*torque/(Istar(time)*Msun*Rsun*Rsun*Omegasun))-(Omega*(J(time)/Istar(time)))
      #ode=((tg*yr2sec)*torque/(Istar(time)*Msun*Rsun*Rsun*Omegasun))-(Omega*((Madot/Mstar)+(2*Rdot(time)/R(time))))

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
  #Prot=8 #days
  tdisk=Tdisk/tg
  Omegasat=10
  Kconstant=2.5e-4

  omega=2.0*np.pi/(Prot*Omegasun*86400)

  Omega_star, times, torA, torM, RCO, RTR, RAST, MAST, ff, TIMER, torW = solution(omega, time[0], tdisk, Bfield, betta, gamma, APSW) #Omega,time0,disktime,Bfield,betta,gamma

  period=2.0*np.pi/(Omega_star*Omegasun*86400)
  vrot=Omega_star[:-1]*Omegasun*RAST*Rsun/1e5

  #################################################################################################


  return TIMER*tg, vrot


""" """
# Create a function to sample random x values from the distribution
def sample_random_x(cdf_values):
    # Generate a random value between 0 and 1
    random_value = random.uniform(0, 1)

    # Find the bin where the random value falls
    for i in range(len(cdf_values)):
        if random_value <= cdf_values[i]:
            return i

"""**Model Generator**"""

def compute(params):
    global tbin
    m, Protin, Maccin, B, A = params
    Models = REFUGEE(Mass=m, Prot=Protin, Macc=Maccin, Tdisk=15e6, Bfield=B, betta=0.01, gamma=1, APSW=A)
    time = Models[0]
    vsini = Models[1] * np.pi / 4
    function = interp1d(time, vsini, kind='linear')
    if tbin==0.5:
        value = function(5.25e5)
    if tbin==1.5:
        value = function(1.5e6)
    if tbin==2.5:
        value = function(2.5e6)
    if tbin==8:
        value = function(8e6)
    return Protin, Maccin, B, value


def combine_histograms(min_bin_edges, max_bin_edges, bin_number, repetitions, histos):
    """ plot as plt.bar(output1, output2, width=output3, color='red', alpha=0.2, yerr=output4, capsize=5, align='edge', label='Models')
    """
    # Calculate histograms for each repetition with a consistent minimum value
    histograms = [np.histogram(np.array(histos[i]), bins=np.linspace(min_bin_edges, max_bin_edges, bin_number)) for i in range(repetitions)]

    # Determine the maximum bin count based on the new histograms
    max_bins = len(histograms[0][1])

    # Extract bin edges for each new histogram
    bin_edges_a = [histogram[1] for histogram in histograms]

    # Zero-pad histograms to the maximum bin count
    padded_histograms = [np.pad(histogram[0], (0, max_bins - len(histogram[0])), 'constant') for histogram in histograms]

    # Calculate the median histogram
    median_histogram = np.median(padded_histograms, axis=0)
    std_histogram = np.std(padded_histograms, axis=0)

    # Calculate the width of each bar
    bar_width = (bin_edges_a[0][-1] - bin_edges_a[0][0]) / (max_bins - 1)

    # Zero-pad all bin_edges_ to the maximum bin count
    padded_bin_edges = [np.pad(bin_edges, (0, max_bins - len(bin_edges)), 'constant') for bin_edges in bin_edges_a]

    # Calculate the mean of the padded bin edges
    padded_bin_centers = np.mean(padded_bin_edges, axis=0)

    return padded_bin_centers, median_histogram, bar_width, std_histogram


# Parallelize the computation
if __name__ == '__main__':

    ages=[0.5, 1.5, 2.5, 8]
    num_processes = 60  # Adjust the number of processes as needed
    repetitions = 100
    for tbin in ages:
        Bfield=np.linspace(500, 3500, 300)
        result = [[] for _ in range(repetitions)]
        vhis = [[] for _ in range(repetitions)]
        if tbin==0.5:
            num_samples = 62  #[62 68 40 30] # Change this to the number of samples you want
            cdf_values = norm.cdf(Bfield, loc=1000, scale=200)
            Bin1_Macc=np.genfromtxt('mdot_bin1_logmdot', comments='#', dtype='S')
            filename_1="./Results/bin1_vsini_result.txt"
            plot_="./Plots/bin1_vsini_model_constant_1.png"
            plot_2="./Plots/bin1_B_model_constant_1.png"
            plot_3="./Plots/bin1_Prot_model_constant_1.png"
            plot_4="./Plots/bin1_Macc_model_constant_1.png"
            filename_2="./Results/bin1_parameters_result.txt"
            databin = np.loadtxt('./Plots/vsini_bin1.txt')

        if tbin==1.5:
            num_samples = 68  #[62 68 40 30] # Change this to the number of samples you want
            cdf_values = norm.cdf(Bfield, loc=1000, scale=200)
            Bin1_Macc=np.genfromtxt('mdot_bin1_logmdot', comments='#', dtype='S')
            filename_1="./Results/bin2_vsini_result.txt"
            plot_="./Plots/bin2_vsini_model_constant_1.png"
            plot_2="./Plots/bin2_B_model_constant_1.png"
            plot_3="./Plots/bin2_Prot_model_constant_1.png"
            plot_4="./Plots/bin2_Macc_model_constant_1.png"
            filename_2="./Results/bin2_parameters_result.txt"
            databin = np.loadtxt('./Plots/vsini_bin2.txt')

        if tbin==2.5:
            num_samples = 40  #[62 68 40 30] # Change this to the number of samples you want
            cdf_values = norm.cdf(Bfield, loc=1000, scale=200)
            Bin1_Macc=np.genfromtxt('mdot_bin1_logmdot', comments='#', dtype='S')
            filename_1="./Results/bin3_vsini_result.txt"
            plot_="./Plots/bin3_vsini_model_constant_1.png"
            plot_2="./Plots/bin3_B_model_constant_1.png"
            plot_3="./Plots/bin3_Prot_model_constant_1.png"
            plot_4="./Plots/bin3_Macc_model_constant_1.png"
            filename_2="./Results/bin3_parameters_result.txt"
            databin = np.loadtxt('./Plots/vsini_bin3.txt')

        if tbin==8:
            num_samples = 30  #[62 68 40 30] # Change this to the number of samples you want
            cdf_values = norm.cdf(Bfield, loc=1000, scale=200)
            Bin1_Macc=np.genfromtxt('mdot_bin1_logmdot', comments='#', dtype='S')
            filename_1="./Results/bin4_vsini_result.txt"
            plot_="./Plots/bin4_vsini_model_constant_1.png"
            plot_2="./Plots/bin4_B_model_constant_1.png"
            plot_3="./Plots/bin4_Prot_model_constant_1.png"
            plot_4="./Plots/bin4_Macc_model_constant_1.png"
            filename_2="./Results/bin4_parameters_result.txt"
            databin = np.loadtxt('./Plots/vsini_bin4.txt')

        bin_=Bin1_Macc.astype(float)
        central_values=bin_#+ np.log10(np.exp((tbin-0.5)/2.1))
        bin_heights, bin_edges = np.histogram(central_values, bins=1000)
        Macc_=np.linspace(min(bin_edges[1:]),max(bin_edges[1:]),100)

        for k in range(repetitions):
            APSW=[0.3]
            Pin=[1,2,3,4,5,6,7,8]
            Mass=[0.5]
            #############################################################################
            random_B_samples = [sample_random_x(cdf_values) for _ in range(num_samples)]
            random_values_B = [Bfield[i] for i in random_B_samples if i is not None]
            #############################################################################
            cumulative_sum = np.cumsum(bin_heights)
            total_height = cumulative_sum[-1]  # Total area under the histogram
            cdf = cumulative_sum / total_height
            cdf_interpol = interp1d(bin_edges[1:],cdf,kind='linear')
            cdf_values_2 = cdf_interpol(Macc_)
            Macc = pow(10,Macc_)
            #################################################################
            random_Macc_samples = [sample_random_x(cdf_values_2) for _ in range(num_samples)]
            random_values_Macc = [Macc[i] for i in random_Macc_samples if i is not None]
            #############################################################################
            params_list = [(Mass[0], random.choice(Pin), j, i, APSW[0]) for i, j in zip(random_values_B,random_values_Macc)]

            with Pool(num_processes) as pool:
                results = pool.map(compute, params_list)

            # Separate the results
            for item in results:
                Protin, Maccin, B, value = item
                result[k].append((Protin, Maccin, B))
                vhis[k].append(value)

        result = np.array(result)
        per_histograms = [result[i][:, 0] for i in range(repetitions)]
        per = np.concatenate(per_histograms)
        B_histograms = [result[i][:, 2] for i in range(repetitions)]
        B = np.concatenate(B_histograms)
        Ma_histograms = [result[i][:, 1] for i in range(repetitions)]
        Ma = np.concatenate(Ma_histograms)

    ########################################################################################################

        comb_1, comb_2, comb_3, comb_4 = combine_histograms(-3.5,100,16, repetitions, vhis)
        comb_1_B, comb_2_B, comb_3_B, comb_4_B = combine_histograms(500,3500,14, repetitions, B_histograms)
        comb_1_Ma, comb_2_Ma, comb_3_Ma, comb_4_Ma = combine_histograms(-10,-6, 10, repetitions, np.log10(Ma_histograms))
        comb_1_P, comb_2_P, comb_3_P, comb_4_P = combine_histograms(0, 8, 9, repetitions, per_histograms)

        plt.figure(1)
        if tbin==0.5:
            plt.title(f'Bin 1')
        if tbin==1.5:
            plt.title(f'Bin 2')
        if tbin==2.5:
            plt.title(f'Bin 3')
        if tbin==8:
            plt.title(f'Bin 4')
        plt.bar(comb_1_B, comb_2_B, width=comb_3_B, color='blue', alpha=0.2, yerr=comb_4_B, capsize=5, align='edge', label='Models')
        plt.ylim(0,)
        plt.xlabel(r"$B_{\ast}\,\,\,(G)$")
        plt.savefig(plot_2,dpi=300)
        plt.cla()
        plt.close()


        plt.figure(2)
        if tbin==0.5:
            plt.title(f'Bin 1')
        if tbin==1.5:
            plt.title(f'Bin 2')
        if tbin==2.5:
            plt.title(f'Bin 3')
        if tbin==8:
            plt.title(f'Bin 4')
        plt.bar(comb_1_P, comb_2_P, width=comb_3_P, color='red', alpha=0.2, yerr=comb_4_P, capsize=5, align='edge', label='Models')
        plt.ylim(0,)
        plt.xlabel(r"$P^{in}_{rot}\,\,\,(d)$")
        plt.savefig(plot_3,dpi=300)
        plt.cla()
        plt.close()


        plt.figure(3)
        if tbin==0.5:
            plt.title(f'Bin 1')
        if tbin==1.5:
            plt.title(f'Bin 2')
        if tbin==2.5:
            plt.title(f'Bin 3')
        if tbin==8:
            plt.title(f'Bin 4')
        plt.hist(bin_, bins=np.linspace(-10,-6, 10), cumulative=False, density=False, alpha=0.2, label='Observations')
        plt.bar(comb_1_Ma, comb_2_Ma, width=comb_3_Ma, color='orange', alpha=0.2, yerr=comb_4_Ma, capsize=5, align='edge', label='Models')
        plt.ylim(0,)
        plt.xlabel(r"$\log(\dot{M}_{acc})\,\,\,(\frac{M_{\ast}}{yr})$")
        plt.legend()
        plt.savefig(plot_4,dpi=300)
        plt.cla()
        plt.close()

        # Find the minimum and maximum value of bin_edges among all histograms
        min_bin_edges = -3.5#min(min(bin_edges) for bin_edges in bin_edges_)
        max_bin_edges = 100#max(max(bin_edges) for bin_edges in bin_edges_)
        bin_step = 16

        # Create the bar plot with error bars
        kk = np.histogram(databin, bins=np.linspace(min_bin_edges, max_bin_edges, bin_step))
        ks_statistic2, ks_pvalue2 = ks_2samp(kk[0], comb_2[:-1])

        plt.figure(4)
        if tbin==0.5:
            plt.title(f'Bin 1 (sample:{num_samples:.0f}, runs:{repetitions:.0f})')
        if tbin==1.5:
            plt.title(f'Bin 2 (sample:{num_samples:.0f}, runs:{repetitions:.0f})')
        if tbin==2.5:
            plt.title(f'Bin 3 (sample:{num_samples:.0f}, runs:{repetitions:.0f})')
        if tbin==8:
            plt.title(f'Bin 4 (sample:{num_samples:.0f}, runs:{repetitions:.0f})')
        plt.bar(comb_1, comb_2, width=comb_3, color='red', alpha=0.2, yerr=comb_4, capsize=5, align='edge', label='Models')
        #plt.bar(padded_bin_centers, median_histogram, width=bar_width, color='red', alpha=0.2, yerr=std_histogram, capsize=5, align='edge', label='Models')
        y2, _ , _ = plt.hist(databin, bins=np.linspace(min_bin_edges, max_bin_edges, bin_step), alpha=0.2, label='Observations',color='darkgreen')
        plt.text(50, max(y2)/2, f'K-S: {ks_statistic2:.3f}', fontsize=8)
        plt.text(50, max(y2)/2.3, f'P-Value: {ks_pvalue2:.3f}', fontsize=8)
        plt.ylim(0,)
        plt.xlim(-3,100)
        plt.xlabel(r"$v\sin(i)\,\,\,[\frac{km}{s}]$")
        plt.legend()
        plt.savefig(plot_,dpi=300)
        plt.cla()
        plt.close()
    #plt.show()

email_sender()
#np.savetxt(filename_1, vhis, header="vsini")
#np.savetxt(filename_2, result, header="Protin Maccin Bfield")
