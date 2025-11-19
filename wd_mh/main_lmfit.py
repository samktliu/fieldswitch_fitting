import numpy as np
import lmfit
import matplotlib.pyplot as plt
import corner
import scienceplots

plt.style.use(['science','nature'])

#### flags #####

mumax = True
gen_emcees = False

#### define constants #####

if mumax:
    samprate = int(1E9)
else:
    samprate = int(2E3)
deltat = 1/samprate 
kb = 1.380649e-16 # erg/K
Temp = 300
kbT = kb*Temp
mu0 = 4*np.pi*1E-7
diam = 65e-7 # cm
Ms = 1495 # emu/cm^3
t_fl = 1.5e-7 # cm

#### define fitting function #####

def ph_fit(h,Hk,fatt,Aex,delta):

    # edw calculation
    edw = np.sqrt(8*Ms*Hk*Aex)

    # q1 calculation
    eps = edw/(Ms*np.abs(h)*diam)
    q1 = 1 + eps - np.sqrt(1+eps*eps)
    q1_p = q1 + delta
    q1_n = q1 - delta

    # thetaq calculation
    thetaq1 = np.atan(q1*(1-q1/2)/(1-q1))
    thetaq1_p = np.atan(q1_p*(1-q1_p/2)/(1-q1_p))
    thetaq1_n = np.atan(q1_n*(1-q1_n/2)/(1-q1_n))

    thetaq1 = np.where(q1>1,np.pi-np.atan(q1*(1-q1/2)/(q1-1)),np.atan(q1*(1-q1/2)/(1-q1)))
    thetaq1_p = np.where(q1_p>1,np.pi-np.atan(q1_p*(1-q1_p/2)/(q1_p-1)),np.atan(q1_p*(1-q1_p/2)/(1-q1_p)))
    thetaq1_n = np.where(q1_n>1,np.pi-np.atan(q1_n*(1-q1_n/2)/(q1_n-1)),np.atan(q1_n*(1-q1_n/2)/(1-q1_n)))

    # Ad calculation
    Adq1_p = 0.25*diam*diam*(thetaq1_p-np.tan(thetaq1_p)+(np.pi/2-thetaq1_p)*np.tan(thetaq1_p)*np.tan(thetaq1_p))
    Adq1_n = 0.25*diam*diam*(thetaq1_n-np.tan(thetaq1_n)+(np.pi/2-thetaq1_n)*np.tan(thetaq1_n)*np.tan(thetaq1_n))
    
    # Ldw calculation
    Ldq1 = diam*(np.pi/2-thetaq1)*np.tan(thetaq1)
    
    # Eq calculation
    Eq1 = Ldq1*edw*t_fl+np.abs(h)*Ms*t_fl*(((diam*diam*np.pi)/4)-Adq1_p-Adq1_n)
    Eb = Eq1-((np.pi/4)*diam*diam*Ms*h*t_fl)
    
    # Calculate P(H)
    ph = 1-np.exp(-fatt*deltat*np.exp(-(Eb/kbT)))
    
    return ph

#### MAIN #####

# load model
ph_model = lmfit.Model(ph_fit)

# initialize fitting parameters
params = lmfit.Parameters()
params.add('Aex',value=1.5E-6,min=1E-6)                   # erg/cm
params.add('delta',value=0.2,min=0.05)     
if mumax:              
    params.add('fatt',value=5E9,min=0.5E9,max=1E10)         # Hz
    params.add('Hk',value=2.8E3,min=1.0E3,max=5E3)        # Oe
else:
    params.add('fatt',value=1E9,vary=False)         # Hz
    params.add('Hk',value=1.81E3,vary=False)        # Oe

# load data
data = np.load('pmm_A1.776.npy')
Hz = data[:,0]
ph_p = data[:,1]

# run fit
result = ph_model.fit(data=ph_p,h=Hz,params=params,method='powell')
if result.success:
    print(result.fit_report())
else:
    print(f'Fit failed')

# generate emcees
if gen_emcees:
    emcee_kws = dict(steps=5000,burn=500,thin=20,is_weighted=False,progress=False)
    emcee_params = result.params.copy()
    emcee_params.add('__lnsigma',value=np.log(0.1),min=np.log(0.001),max=np.log(2.0))
    result_emcee = ph_model.fit(data=ph_p,h=Hz,params=emcee_params,method='emcee',nan_policy='omit',fit_kws=emcee_kws)
    print(result_emcee.fit_report())

# generate fitting data
hfit = np.linspace(1,np.max(Hz),1000)
pfit = ph_model.eval(result.params,h=hfit)

# plot fit and data
fig,ax = plt.subplots(1,1,figsize=(3,2))
ax.plot(Hz,ph_p,'^',color='black',mfc='none')
ax.plot(hfit,pfit,color='tab:red')
fig.savefig('fitting_lmfit.svg')

# plot emcee corner plot
if gen_emcees:
    plt.figure()
    emcee_corner = corner.corner(result_emcee.flatchain,labels=result_emcee.var_names,
                                truths=list(result_emcee.params.valuesdict().values()))
    plt.savefig('fitting_emcee.svg')