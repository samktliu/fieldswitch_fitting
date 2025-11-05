# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 16:22:40 2025

@author: pufall
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May  8 22:35:36 2025
DW energy barrier function for fitting probability  
@author: pufall
"""


import numpy as np
import matplotlib.pyplot as plt
import math

#### define constants #####
kb = 1.380649E-23  #J/K
Temp = 300
kbT = kb*Temp
fatt = 1E9
mu0 = 4*np.pi*1E-7
class Ebclass:
    # set of equations to define energy barrier due to domain wall formation and switching
    # from Appl. Phys. Lett. 117, 242404 (2020); doi: 10.1063/5.0023852
    # This is supposed to be in SI units. The paper's units are cgs, so a 4pi and/or a mu0 might be missing. I think Ms should be in A/m, and a mu0 should come into either H or Ms... 
    # this class has functions to give external access to sub-parts of the calculation. 
    def __init__(self,H,edw,wdw,t,D,Ms):
        self.H = H
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D  = D
        self.Ms = Ms
        self.delta = wdw/D
        self.eps = edw/(mu0*Ms*np.abs(H)*D) # eqn 2g
        self.q1 = 1 + self.eps - np.sqrt(1+self.eps*self.eps) # eqn 2f. key approx to solve this thing
        
    def thetaqf(self, q):
        # Eqn 2a
        theta1 = np.atan(q*(1-q/2)/(1-q))
        theta2 = np.atan(q*(1-q/2)/(q-1))
        # if q <1:
        #     retVal = thetaq
        # elif q>=1 and q < 2:
        #     retVal = np.pi-thetaq
        # else:
        #     print(f'q = {q}!!')
        retVal = np.where(q>1,np.pi-theta2, theta1)
        return retVal
    
    def Adtheta(self,thet):
        # Eqn 2b
        # need to trap for q- delta less than zero
        theVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        return theVal

    def Adq(self, q):
        # Eqn 2b
        # in terms of q. need to trap for q <0
        thet = self.thetaqf(q)
        AdqVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        theVal = np.where(q>0,AdqVal, 0)
#        theRetVal = np.where(AdqVal>=0,AdqVal,0)
        #theRetVal = np.where(theVal>=0,theVal,0)
        return theVal
    
    def Ldw(self, q):
        # Eqn 2        
        thet = self.thetaqf(q)
        Ldw = self.D*(np.pi/2 - thet)*np.tan(thet)
        return Ldw        

    def Eq(self):   
# calc areas as func of various q+/-delta. Adq() calculates theta(q) internally
        Adthetap = self.Adq(self.q1+self.delta)  # Areas calculated at q +/- delta = wdw/D
        Adthetam = self.Adq(self.q1-self.delta)
        devArea = np.pi*(self.D/2)**2
        #Ldwtheta = D*(np.pi/2 - thetaq)*np.tan(thetaq) # Eval at theta(q). 
        Ldwtheta = self.Ldw(self.q1) #evaluated at q1, not at q1 +/- delta. 
        EqVal = Ldwtheta*self.edw*self.t + mu0*np.abs(self.H)*self.Ms*self.t*((np.pi*self.D*self.D)/4 - Adthetap - Adthetam)
        return EqVal
    
    def Eb(self):
    #form energy function, and Eb by subtracting unreversed energy 
    # (recall H sign is defined as positive being AP to initial state)
        Ebdw = self.Eq() -mu0*(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

class EbclassH:
    # set of equations to define energy barrier due to domain wall formation and switching
    # from Appl. Phys. Lett. 117, 242404 (2020); doi: 10.1063/5.0023852
    # This is supposed to be in SI units. The paper's units are cgs, so a 4pi and/or a mu0 might be missing. I think Ms should be in A/m, and a mu0 should come into either H or Ms... 
    # this class has functions to give external access to sub-parts of the calculation.
    # remake this to make it an explicit function of H, not fully calculated liek teh above class upon definition.
    # maybe this only involves one additional function setH? try that
    def __init__(self,edw,wdw,t,D,Ms):
        #self.H = H
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D  = D
        self.Ms = Ms
        self.delta = wdw/D
        
    def setH(self, H):
        # calcs new things based on H
        self.H = H
        self.eps = self.edw/(mu0*self.Ms*np.abs(H)*self.D) # eqn 2g
        self.q1 = 1 + self.eps - np.sqrt(1+self.eps*self.eps) # eqn 2f. key approx to solve this thing
        
    def thetaqf(self, q):
        # Eqn 2a
        theta1 = np.atan(q*(1-q/2)/(1-q))
        theta2 = np.atan(q*(1-q/2)/(q-1))
        # if q <1:
        #     retVal = thetaq
        # elif q>=1 and q < 2:
        #     retVal = np.pi-thetaq
        # else:
        #     print(f'q = {q}!!')
        retVal = np.where(q>1,np.pi-theta2, theta1)
        return retVal
    
    def Adtheta(self,thet):
        # Eqn 2b
        # need to trap for q- delta less than zero
        theVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        return theVal

    def Adq(self, q):
        # Eqn 2b
        # in terms of q. need to trap for q <0
        thet = self.thetaqf(q)
        AdqVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        theVal = np.where(q>0,AdqVal, 0)
#        theRetVal = np.where(AdqVal>=0,AdqVal,0)
        #theRetVal = np.where(theVal>=0,theVal,0)
        return theVal
    
    def Ldw(self, q):
        # Eqn 2        
        thet = self.thetaqf(q)
        Ldw = self.D*(np.pi/2 - thet)*np.tan(thet)
        return Ldw        

    def Eq(self):   
# calc areas as func of various q+/-delta. Adq() calculates theta(q) internally
        Adthetap = self.Adq(self.q1+self.delta)  # Areas calculated at q +/- delta = wdw/D
        Adthetam = self.Adq(self.q1-self.delta)
        devArea = np.pi*(self.D/2)**2
        #Ldwtheta = D*(np.pi/2 - thetaq)*np.tan(thetaq) # Eval at theta(q). 
        Ldwtheta = self.Ldw(self.q1) #evaluated at q1, not at q1 +/- delta. 
        EqVal = Ldwtheta*self.edw*self.t + mu0*np.abs(self.H)*self.Ms*self.t*((np.pi*self.D*self.D)/4 - Adthetap - Adthetam)
        return EqVal
    
    def Eb(self, H):
    #form energy function, and Eb by subtracting unreversed energy 
    # (recall H sign is defined as positive being AP to initial state)
        self.setH(H)
        Ebdw = self.Eq() -mu0*(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

    def Ebfunc(self,H,edw,wdw,t,D,Ms):
        # function for fitting.   
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D  = D
        self.Ms = Ms
        self.delta = wdw/D
        self.setH(H)
        Ebdw = self.Eq() -mu0*(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

class Ebclassfunc:
    # set of equations to define energy barrier due to domain wall formation and switching
    # from Appl. Phys. Lett. 117, 242404 (2020); doi: 10.1063/5.0023852
    # This is supposed to be in SI units. The paper's units are cgs, so a 4pi and/or a mu0 might be missing. I think Ms should be in A/m, and a mu0 should come into either H or Ms... 
    # this class has functions to give external access to sub-parts of the calculation.
    # remake this to make it an explicit function of H, not fully calculated liek teh above class upon definition.
    # maybe this only involves one additional function setH? try that
    # def __init__(self,edw,wdw,t,D,Ms):
    #     #self.H = H
    #     self.edw = edw
    #     self.wdw = wdw
    #     self.t = t
    #     self.D  = D
    #     self.Ms = Ms
    #     self.delta = wdw/D
        
    def setH(self, H):
        # calcs new things based on H
        self.H = H
        self.eps = self.edw/(mu0*self.Ms*np.abs(H)*self.D) # eqn 2g
        self.q1 = 1 + self.eps - np.sqrt(1+self.eps*self.eps) # eqn 2f. key approx to solve this thing
        
    def thetaqf(self, q):
        # Eqn 2a
        theta1 = np.atan(q*(1-q/2)/(1-q))
        theta2 = np.atan(q*(1-q/2)/(q-1))
        # if q <1:
        #     retVal = thetaq
        # elif q>=1 and q < 2:
        #     retVal = np.pi-thetaq
        # else:
        #     print(f'q = {q}!!')
        retVal = np.where(q>1,np.pi-theta2, theta1)
        return retVal
    
    def Adtheta(self,thet):
        # Eqn 2b
        # need to trap for q- delta less than zero
        theVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        return theVal

    def Adq(self, q):
        # Eqn 2b
        # in terms of q. need to trap for q <0
        thet = self.thetaqf(q)
        AdqVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        theVal = np.where(q>0,AdqVal, 0)
#        theRetVal = np.where(AdqVal>=0,AdqVal,0)
        #theRetVal = np.where(theVal>=0,theVal,0)
        return theVal
    
    def Ldw(self, q):
        # Eqn 2        
        thet = self.thetaqf(q)
        Ldw = self.D*(np.pi/2 - thet)*np.tan(thet)
        return Ldw        

    def Eq(self):   
# calc areas as func of various q+/-delta. Adq() calculates theta(q) internally
        Adthetap = self.Adq(self.q1+self.delta)  # Areas calculated at q +/- delta = wdw/D
        Adthetam = self.Adq(self.q1-self.delta)
        devArea = np.pi*(self.D/2)**2
        #Ldwtheta = D*(np.pi/2 - thetaq)*np.tan(thetaq) # Eval at theta(q). 
        Ldwtheta = self.Ldw(self.q1) #evaluated at q1, not at q1 +/- delta. 
        EqVal = Ldwtheta*self.edw*self.t + mu0*np.abs(self.H)*self.Ms*self.t*((np.pi*self.D*self.D)/4 - Adthetap - Adthetam)
        return EqVal
    
    def Eb(self, H):
    #form energy function, and Eb by subtracting unreversed energy 
    # (recall H sign is defined as positive being AP to initial state)
        self.setH(H)
        Ebdw = self.Eq() -mu0*(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

    def Ebfunc(self,H,edw,wdw,t,D,Ms):
        # function for fitting.   
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D  = D
        self.Ms = Ms
        self.delta = wdw/D
        self.setH(H)
        Ebdw = self.Eq() -mu0*(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

class basicEb:
    # basic energy barrier function with quadratic field dependence
    def __init__(self,Hk,D,t,Ms):
        self.meow = 1
        self.Hk = Hk
        self.t = t
        self.D = D
        self.Ms = Ms
        self.Keff = mu0*Hk*Ms/2 
        self.Vol = self.t*np.pi*(self.D/2)**2
        
    def Ebfunc(self,H,Hk,D,t,Ms):
        #Eb = self.Keff*self.Vol*(1-H/self.Hk)**2
        self.Keff = mu0*Hk*Ms/2 
        self.Vol = self.t*np.pi*(self.D/2)**2
        Eb = np.where(H>Hk, 0, self.Keff*self.Vol*(1-H/self.Hk)**2)
        return Eb
    
class EbclassHCGS:
    # set of equations to define energy barrier due to domain wall formation and switching
    # from Appl. Phys. Lett. 117, 242404 (2020); doi: 10.1063/5.0023852
    # this class has functions to give external access to sub-parts of the calculation.
    # remake this to make it an explicit function of H, not fully calculated liek teh above class upon definition.
    def __init__(self,edw,wdw,t,D,Ms):
        #self.H = H
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D  = D
        self.Ms = Ms
        self.delta = wdw/D
        
    def setH(self, H):
        # calcs new things based on H
        self.H = H
        self.eps = self.edw/(self.Ms*np.abs(H)*self.D) # eqn 2g
        self.q1 = 1 + self.eps - np.sqrt(1+self.eps*self.eps) # eqn 2f. key approx to solve this thing
        
    def thetaqf(self, q):
        # Eqn 2a
        theta1 = np.atan(q*(1-q/2)/(1-q))
        theta2 = np.atan(q*(1-q/2)/(q-1))
        # if q <1:
        #     retVal = thetaq
        # elif q>=1 and q < 2:
        #     retVal = np.pi-thetaq
        # else:
        #     print(f'q = {q}!!')
        retVal = np.where(q>1,np.pi-theta2, theta1)
        return retVal
    
    def Adtheta(self,thet):
        # Eqn 2b
        # need to trap for q- delta less than zero
        theVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        return theVal

    def Adq(self, q):
        # Eqn 2b
        # in terms of q. need to trap for q <0
        thet = self.thetaqf(q)
        AdqVal = (thet - np.tan(thet)+(np.pi/2 - thet)*np.tan(thet)*np.tan(thet))*self.D*self.D/4
        theVal = np.where(q>0,AdqVal, 0)
        # theRetVal = np.where(AdqVal>=0,AdqVal,0)
        # theRetVal = np.where(theVal>=0,theVal,0)
        return theVal
    
    def Ldw(self, q):
        # Eqn 2        
        thet = self.thetaqf(q)
        Ldw = self.D*(np.pi/2 - thet)*np.tan(thet)
        return Ldw        

    def Eq(self):   
# calc areas as func of various q+/-delta. Adq() calculates theta(q) internally
        Adthetap = self.Adq(self.q1+self.delta)  # Areas calculated at q +/- delta = wdw/D
        Adthetam = self.Adq(self.q1-self.delta)
        devArea = np.pi*(self.D/2)**2
        #Ldwtheta = D*(np.pi/2 - thetaq)*np.tan(thetaq) # Eval at theta(q). 
        Ldwtheta = self.Ldw(self.q1) #evaluated at q1, not at q1 +/- delta. 
        EqVal = Ldwtheta*self.edw*self.t + np.abs(self.H)*self.Ms*self.t*((np.pi*self.D*self.D)/4 - Adthetap - Adthetam)
        return EqVal
    
    def Eb(self, H):
    #form energy function, and Eb by subtracting unreversed energy 
    # (recall H sign is defined as positive being AP to initial state)
        self.setH(H)
        Ebdw = self.Eq()-(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw

    def Ebfunc(self,H,edw,wdw,t,D,Ms):
        # function for fitting.   
        self.edw = edw
        self.wdw = wdw
        self.t = t
        self.D = D
        self.Ms = Ms
        self.delta = wdw/D
        self.setH(H)
        Ebdw = self.Eq()-(np.pi/4)*self.D*self.D*self.Ms*self.H*self.t
        return Ebdw
