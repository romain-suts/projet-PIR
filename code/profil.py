import math
from math import *
import numpy as np
import matplotlib.pyplot as plt

def affichage(tour,g,pression,rayon):
    x = np.linspace(rayon,0,tour) 
    plt.plot (x,g)
    plt.yscale('log')
    plt.show()
    plt.plot(x,pression)
    plt.yscale('log')
    plt.show()
    
def calc_forme(masse,rayon):#modele densite constante 
    densite_glace = 0.917 #glace pure
    densite_silicate = 3#valeur a peut etre redeterminée suivant le noyau que l'on cherche a avoir 
    volume =4/3*pi*np.power(rayon,3)
    densite_lune = masse/(1000*volume)
    print("la densité de la lune est " + str(densite_lune))
    volume_noyau = volume*(densite_lune-densite_glace)/(densite_silicate-densite_glace)
    rayon_noyau = math.cbrt(volume_noyau*3/(4*pi))#!! si noyau metallique present on peut avoir un rayon de silicate > rayon de la lune ex io
    print("le rayon du noyau de la lune est " + str(rayon_noyau))
    return rayon_noyau

masse = float(input("rentrer la masse de la lune en kg"))
rayon_lune = float(input("rentrer le rayon de la lune en m")) 
pas = 0.1
G = 6.6743E-11
g = [masse*G/np.power(rayon_lune,2)]
masse_vol_glace = 917
masse_vol_silicate = 3000
pression = [0] #peut changer si atmosphere
rayon_noyau = calc_forme(masse,rayon_lune)
rayon = rayon_lune
tour = 1
fichier = open("data.txt", "a")
while (rayon>rayon_noyau):#1ere couche
    tour +=1
    rayon-=pas
    g.append(g[-1]-pas*(4*pi*G*masse_vol_glace-2*g[-1]/rayon))
    pression.append(pression[-1]+pas*g[-1]*masse_vol_glace)
    fichier.write("{}\t{}\t{}\n".format(rayon,pression[-1],g[-1]))
while (rayon>=pas):#2eme couche
    tour +=1
    rayon=rayon-pas
    g.append(g[-1]-pas*(4*pi*G*masse_vol_silicate-2*g[-1]/rayon))
    pression.append(pression[-1]+pas*g[-1]*masse_vol_silicate)
    fichier.write("{}\t{}\t{}\n".format(rayon,pression[-1],g[-1]))
affichage(tour,g,pression,rayon_lune)
fichier.close()