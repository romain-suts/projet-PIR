import math
from math import *
import numpy as np
import matplotlib.pyplot as plt
def densite(R,rayon_noyau,etat):
    masse_vol_glace = 917
    masse_vol_eau = 1000
    masse_vol_silicate = 3700
    if (R<=rayon_noyau):
        return masse_vol_silicate
    else:
        if (etat==0):
            return masse_vol_glace
        else:
            return masse_vol_eau
        
def S(r,rayon_noyau):
    if(r<=rayon_noyau):
        return 2E-12
    else:
        return 0

def lambda_(R,rayon_noyau):
    conduc_ther_glace = 2.1
    conduc_ther_silicate = 0.22
    if (R<=rayon_noyau):
        return conduc_ther_silicate
    else:
        return conduc_ther_glace

def calc_I(rayon_lune,rayon_noyau,masse_lune,pas):
    r=0
    I=0
    while (r<=rayon_lune):
        r+=pas
        I+=8/3*pi*r**4*densite(r,rayon_noyau,etat)
    I=I/(np.power(rayon,2)*masse_lune)
    return I

def affichage(tour,g,pression,chaleur,T,I,etat,rayon,rayon_noyau):
    x = np.linspace(rayon,0,tour+1) 
    plt.plot (x,g)
    plt.xlabel("rayon")
    plt.ylabel("g")
    plt.show()
    plt.plot(x,pression)
    plt.xlabel("Rayon")
    plt.ylabel("pression")
    plt.show()
    plt.plot(x,chaleur)
    plt.xlabel("Rayon")
    plt.ylabel("chaleur")
    plt.show()
    plt.plot(x,T)
    plt.xlabel("Rayon")
    plt.ylabel("température")
    plt.show()
    plt.plot(x,I)
    plt.xlabel("Rayon")
    plt.ylabel("moment d'inertie")
    plt.show()
    plt.plot(x,etat)
    plt.xlabel("Rayon")
    plt.ylabel("etat")
    plt.show()

def calc_forme (masse,masse_lune,rayon_noyau):
    masse_vol_glace = 1000
    masse_vol_silicate = 3700
    delta_masse = masse-masse_lune
    delta_volume = (masse_vol_silicate-masse_vol_glace)*delta_masse
    nv_rayon = math.cbrt( np.power(rayon_noyau,3)+3/(4*pi)*delta_volume)
    return nv_rayon

def calc_forme_init(masse,rayon): 
    
    densite_glace = 0.917 #glace pure
    densite_silicate = 3.5#valeur a peut etre redeterminée suivant le noyau que l'on cherche a avoir 
    volume =4/3*pi*np.power(rayon,3)
    densite_lune = masse/(1000*volume)
    print("la densité de la lune est " + str(densite_lune))
    volume_noyau = volume*(densite_lune-densite_glace)/(densite_silicate-densite_glace)
    rayon_noyau = math.cbrt(volume_noyau*3/(4*pi))#!! si noyau metallique present on peut avoir un rayon de silicate > rayon de la lune ex io
    print("le rayon du noyau de la lune est " + str(rayon_noyau))
    return rayon_noyau,densite_lune
def calc_g(rayon,rayon_noyau,G,etat,pas):
    r=0
    masse=0
    while (r<=rayon):
        masse += pas*4*pi*np.power(r,2)*densite(r,rayon_noyau,etat)
        r += pas
    return (G*masse/np.power(rayon,2))

def etat_l(P,T):
    delta_V = -0.0000885
    chaleur_latente = 333.55
    P_atm = 101325
    if (math. log(T)>math. log(273.15)+(delta_V/chaleur_latente)*(P-P_atm)):
        return 1
    else:
        return 0
    
masse_lune = 1.4819E23
rayon_lune = 2631.2E3
moment_inertie = 0.3115
pression_surface = 0. 
temperature_surface = 110.
chaleur_surface = 0.002
pas = 200.
G = 6.6743E-11
g_lune = masse_lune*G/np.power(rayon_lune,2)
print(g_lune)
tour_global = 0
masse = 0.
rayon_noyau,densite_lune = calc_forme_init(masse_lune,rayon_lune)
nb_iteration_max = 40

while (tour_global<nb_iteration_max):
    tour_global += 1
    masse = 0
    rayon = rayon_lune
    tour = 0
    etat = [0]
    g = [g_lune]
    pression = [pression_surface]
    T = [temperature_surface]
    chaleur = [chaleur_surface]
    I = [0]
    calibration = 0
    
    while (rayon>2*pas):#probleme de g a regler
        tour += 1
        rayon -= pas
        masse += pas*4*pi*np.power(rayon,2)*densite(rayon,rayon_noyau,etat[-1])
        if ((masse_lune-masse)>0 and rayon>=0.1*rayon_lune):
            g.append((masse_lune-masse)*G/(rayon**2))
        elif(rayon<0.2*rayon_lune):
            g.append(calc_g(rayon,rayon_noyau,G,0,pas))
        else:
            g.append(0)
        pression.append(pression[-1]+pas*g[-1]*densite(rayon,rayon_noyau,etat[-1]))
        chaleur.append(chaleur[-1]+pas*(densite(rayon,rayon_noyau,etat)*S(rayon,rayon_noyau) - 2*chaleur[-1]/rayon))
        T.append(T[-1]+pas*(chaleur[-1]/lambda_(rayon,rayon_noyau)))
        I.append(I[-1]-(pas*8/3*np.pi*rayon**4*densite(rayon,rayon_noyau,etat[-1])/(masse_lune*rayon_lune**2)))
        
        if (rayon>rayon_noyau):
            etat.append(etat_l(pression[-1],T[-1]))
        else:
            etat.append(2)
            
    while (rayon>0):
        rayon -= pas
        masse += pas*4*pi*np.power(rayon,2)*densite(rayon,rayon_noyau,etat[-1])
        
    for i in range(len(I)):
        
        I[i] -= I[-1] #redressement moment d'inertie 
        
    masse_excedentaire = masse-masse_lune
    calibration = ((rayon_noyau**3)+(3*masse_excedentaire)/(4*pi*(densite(rayon_lune,rayon_noyau,1)-densite(0,rayon_noyau,2))))**(1/3)
    print(calibration)
    rayon_noyau = calibration

affichage(tour,g,pression,chaleur,T,I,etat,rayon_lune,rayon_noyau)