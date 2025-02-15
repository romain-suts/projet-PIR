import math
from math import *
import numpy as np
import matplotlib.pyplot as plt

def densite(R,rayon_noyau,etat,P,T):#fonction qui retourne la masse volumique de la couche 
    masse_vol_glace_I = 917#etablissement des constantes a des pressions et temperature de reference 
    masse_vol_glace_III = 1160
    masse_vol_glace_V = 1240
    masse_vol_glace_VI = 1310
    masse_vol_glace_VII = 1650
    masse_vol_eau = 1000
    masse_vol_silicate = 3700
    compression_isotherme_eau = 2.2*10**9
    P_atm = 101325
    if (R<=rayon_noyau):#si on est dans le noyau => silicate
        return masse_vol_silicate
    else: #si on est dans le manteau glace => revoie la masse volumique de la glace/eau de la couche
        if (etat==1):
            return masse_vol_glace_I
        elif (etat==3):
            return masse_vol_glace_III
        elif (etat==5):
            return masse_vol_glace_V
        elif (etat==6):
            return masse_vol_glace_VI
        elif (etat==7):
            return masse_vol_glace_VII
        else:
            return masse_vol_eau*(1+(P-P_atm)/compression_isotherme_eau)*(1+2.07*10**(-4)*(T-293.15))
        
def S(r,rayon_noyau):#renvoi la puissance généré par kg de la couche
    if(r<=rayon_noyau):#puissance emise pas le silicate
        return 2E-12
    else:#la glace n'emet pas de puissance radiogénique 
        return 0

def lambda_(R,rayon_noyau):#renvoi la conduction thermique de la couche 
    conduc_ther_glace = 2.1
    conduc_ther_silicate = 0.22
    if (R<=rayon_noyau):
        return conduc_ther_silicate
    else:
        return conduc_ther_glace

def affichage(tour,g,pression,chaleur,T,I,etat,masse_vol,rayon):#fonction qui affiche les graphiques
    x = np.linspace(rayon,0,tour+1) 
    plt.grid()
    plt.plot (x,g)
    plt.xlabel("rayon")
    plt.ylabel("g")
    plt.show()
    plt.grid()
    plt.plot(x,pression)
    plt.xlabel("Rayon")
    plt.ylabel("pression")
    plt.show()
    plt.grid()
    plt.plot(x,chaleur)
    plt.xlabel("Rayon")
    plt.ylabel("chaleur")
    plt.show()
    plt.grid()
    plt.plot(x,T)
    plt.xlabel("Rayon")
    plt.ylabel("température")
    plt.show()
    plt.grid()
    plt.plot(x,I)
    plt.xlabel("Rayon")
    plt.ylabel("moment d'inertie")
    plt.show()
    plt.grid()
    plt.plot(x,etat)
    plt.xlabel("Rayon")
    plt.ylabel("etat")
    plt.show()
    plt.grid()
    plt.plot(T,pression)
    plt.xlabel("T")
    plt.ylabel("P")
    plt.yscale("log") 
    plt.show()
    plt.plot(x,masse_vol)
    plt.xlabel("Rayon")
    plt.ylabel("masse volumique")
    plt.show()

def calc_forme_init(masse,rayon): #fonction qui intialise la premiere forme de la lune pour des masses volumique constante (on cherche a retrouver la masse de la lune)
    
    densite_glace = 0.917 #glace pure
    densite_silicate = 3.7#valeur a peut etre redeterminée suivant le noyau que l'on cherche a avoir !! penser a changer les valeurs pour les autres fonctions 
    volume =4/3*pi*np.power(rayon,3)
    densite_lune = masse/(1000*volume)
    print("la densité de la lune est " + str(densite_lune))
    volume_noyau = volume*(densite_lune-densite_glace)/(densite_silicate-densite_glace)
    rayon_noyau = math.cbrt(volume_noyau*3/(4*pi))#!! si noyau metallique present on peut avoir un rayon de silicate > rayon de la lune ex io
    print("le rayon du noyau de la lune est " + str(rayon_noyau))
    return rayon_noyau,densite_lune

def calc_g(rayon,G,pas,masse_vol):#calcule la valeur de g vers le centre du noyau grace a la formule theorique
    r=0
    masse=0
    while (r<=rayon):
        masse += pas*4*pi*np.power(r,2)*masse_vol
        r += pas
    return (G*masse/np.power(rayon,2))

def etat_l(P,T):#determine l'etat de la couche 
    delta_V_ice_I = -0.0000885
    chaleur_latente_I = 333.55
    P_atm = 101325
    if ((math. log(T)<math. log(273.15)+(delta_V_ice_I/chaleur_latente_I)*(P-P_atm)) and (P<209.9*10**6)):#glace 1 
        return 1
    elif ((P*10**(-6)>209.5+101.1*((T/251.15)**(42.86)-1)) and (P<350.1*10**6) and (P>=209.9*10**6)):#glace 3
        return 3
    elif((P*10**(-6)>355+373.6*((T/256.43)**(8.66)-1)) and (P<632.4*10**6) and (P>=350.1*10**6)):#glace 5
        return 5
    elif((P*10**(-6)>618.4+661.4*((T/272.73)**(4.69)-1)) and (P<2.216*10**9) and (P>=632.4*10**6) ):#glace 6
        return 6
    elif ((P*10**(-9)>2.67348*10**(-4)*T**(1.55299)-0.22933) and (P>=2.216*10**9)): #glace 7
        return 7
    else:
        return 0 #eau liquide
    
masse_vol_glace = 917#intialisation des constantes
masse_vol_silicate = 3700
masse_lune = 1.4819E23
rayon_lune = 2631.2E3
moment_inertie = 0.3115
pression_surface = 0. 
temperature_surface = 110.
chaleur_surface = 0.002
pas = 1000.
G = 6.6743E-11
g_lune = masse_lune*G/np.power(rayon_lune,2)
print(g_lune)
tour_global = 0
masse = 0.
rayon_noyau,densite_lune = calc_forme_init(masse_lune,rayon_lune)
nb_iteration_max = 10

while (tour_global<nb_iteration_max):#boucle de calcule
    
    tour_global += 1#initialisation des constantes pour une nouvelle iteration
    masse = 0
    rayon = rayon_lune
    tour = 0
    etat = [1]
    g = [g_lune]
    pression = [pression_surface]
    T = [temperature_surface]
    chaleur = [chaleur_surface]
    I = [0]
    calibration = 0
    lim_g_estimation = 0.2*rayon_noyau #limite du basculement de methode pour calculer g
    masse_vol = [densite(rayon,rayon_noyau,etat[-1],pression[-1],T[-1])]
    
    while (rayon>pas):
        tour += 1
        rayon -= pas
        masse += pas*4*pi*np.power(rayon,2)*masse_vol[-1]
        #iteration de g
        if ((masse_lune-masse)>0 and rayon>=lim_g_estimation):#methode direct
            g.append((masse_lune-masse)*G/(rayon**2))
        elif(rayon<lim_g_estimation):#zone dans laquelle la premiere methode semble invalide du aux erreurs dans le calcul de la masse
            g.append(calc_g(rayon,G,pas,masse_vol[-1]))
        else:
            g.append(0)#si probleme renvoie la valeur 0
        pression.append(pression[-1]+pas*g[-1]*masse_vol[-1])
        chaleur.append(chaleur[-1]+pas*(masse_vol[-1]*S(rayon,rayon_noyau) - 2*chaleur[-1]/rayon))
        #evolution de la temperature
        if (len(etat)>1 and ((etat[-1] != 0 and etat[-2]==0) or (etat[-1] == 0 and etat[-2]!=0))):#on cherche des interface liquide => loi de newton
            T.append(T[-1]+chaleur[-1]/10)#1000 a redefinir!!!
        elif(len(etat)>1 and (etat[-1] == 0 and (etat[-2] == 0))):#cherche a savoir si on est dans de l'eau loin des parois => convection
            T.append(T[-1])
        else:#on est dans une zone solide => diffusion
            T.append(T[-1]+pas*(chaleur[-1]/lambda_(rayon,rayon_noyau)))

        I.append(I[-1]-(pas*8/3*np.pi*rayon**4*masse_vol[-1]/(masse_lune*rayon_lune**2))) #incrementation du moment d'inertie de la couche
        
        if (rayon>rayon_noyau):#etablissement de l'etat de la couche
            etat.append(etat_l(pression[-1],T[-1]))
        else:
            etat.append(-1)#on revoie -1 comme valeur pour le silicate
            
        masse_vol.append(densite(rayon,rayon_noyau,etat[-1],pression[-1],T[-1]))
    for i in range(len(I)):
        
        I[i] -= I[-1] #redressement moment d'inertie 
        
    masse_excedentaire = masse-masse_lune
    calibration = ((rayon_noyau**3)+(3*masse_excedentaire)/(4*pi*(masse_vol_glace-masse_vol_silicate)))**(1/3)
    print(masse_excedentaire)
    rayon_noyau = calibration

affichage(tour,g,pression,chaleur,T,I,etat,masse_vol,rayon_lune)