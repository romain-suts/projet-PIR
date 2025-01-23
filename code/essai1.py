# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:02:50 2025

@author: louau
"""

#Importation des modules

from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# DEFINITIONS DES CONSTANTES
densite_glace = 917 # kg/m^3
densite_silicate = 3000 # kg/m^3
G = 6.674*10**(-11) #SI

# PARAMETRES DE L'OBJET ETUDIE : Europa
pourcentage_glace = 40/100
pourcentage_silicate = 60/100
rayon = 1510.8*10**3 # m
masse = 4.8 *10**22 # kg

def calcul_masse():
    densite_corps = pourcentage_glace*densite_glace + pourcentage_silicate*densite_silicate
    volume_corps = 4/3*np.pi*rayon**3 
    masse_corps = densite_corps*volume_corps
    return masse_corps , densite_corps

masse_corps = calcul_masse()[0]
densite_corps = calcul_masse()[1]

print("La masse d'Europa sous nos hypothèses est de {} kg.".format(masse_corps))
print("Ce qui représente un écart de {} %.".format(abs((masse-masse_corps)/masse_corps)))

# HYPOTHESE : LES COUCHES SONT DIFFERENCIEES

def calcul_rayon_silicate():
    masse_silicate = masse_corps*pourcentage_silicate
    
    def drdm(m, r):
        return 1 / (4 * np.pi * r**2 * densite_silicate)
    

    m_span = (0, masse_silicate)   #intervalle
    r0 = [0.001]  #condition initiale
    m_eval = np.linspace(0, masse_silicate, 100)
    
    # Résolution de l'équation différentielle
    solution = solve_ivp(drdm, m_span, r0, t_eval=m_eval)
    
    return solution.t, solution.y[0]


distribution_masse_silicate, distribution_rayon_silicate = calcul_rayon_silicate()

print("Le rayon du noyau de silicate vaut {} km.".format(distribution_rayon_silicate[99]*10**(-3)))

def calcul_coquille_glace():
    masse_glace = masse_corps*densite_glace
    
    def drdm(m,r):
        return 1 / (4 * np.pi * r**2 * densite_glace)

    m_span = (distribution_masse_silicate[99], masse_corps)   #intervalle
    r0 = [distribution_rayon_silicate[99]]  #condition initiale
    m_eval = np.linspace(distribution_masse_silicate[99], masse_corps, 100)
    
    # Résolution de l'équation différentielle
    solution = solve_ivp(drdm, m_span, r0, t_eval=m_eval)
    
    return solution.t, solution.y[0]

distribution_masse_glace, distribution_rayon_glace = calcul_coquille_glace()

# Concaténation des distributions de masse
distribution_masse_totale = np.concatenate((distribution_masse_silicate, distribution_masse_glace))

# Concaténation des distributions de rayon correspondantes
distribution_rayon_totale = np.concatenate((distribution_rayon_silicate, distribution_rayon_glace))


plt.figure(1)
plt.plot(distribution_rayon_totale*10**(-3),distribution_masse_totale)
plt.plot([distribution_rayon_silicate[99]*10**(-3),distribution_rayon_silicate[99]*10**(-3)],[0,masse_corps],color="red")
plt.xlabel("Rayon (en km)")
plt.ylabel("Masse (en kg)")
plt.title("Evolution de la masse au sein du corps en fonction du rayon")
plt.grid()

# CALCUL DE LA GRAVITE EN TOUT POINT 
# def calcul_gravite():
#     g=np.zeros(200)
#     for i in range(200):
#         g[i]=G*distribution_masse_totale[i]/distribution_rayon_totale[i]**2
#     return g

# g = calcul_gravite()

# plt.figure(2)
# plt.plot(distribution_rayon_totale*10**(-3),g)
# plt.grid()

#retourne R
# def ech(tab):
#     # Parcours du tableau avec un pas de 2
#     for i in range(0, len(tab) - 1, 1):
#         # Échange de tab[i] et tab[99-i]
#         tab[i], tab[99-i] = tab[99-i], tab[i]
#     return tab

# # EQUATION QUI REGIE LA PRESSION EN FONCTION DU RAYON

# def calcul_pression():
#     def dPdr(P,r):
#         return - densite_silicate*g
    
#     masse_inverse = ech(distribution_masse_silicate)
#     rayon_inverse = ech(distribution_rayon_silicate)
    
#     r_span = (rayon_inverse[0], rayon_inverse[99])   #intervalle
#     P0 = [0.001] #condition initiale
#     r_eval = np.linspace(rayon_inverse[0], rayon_inverse[99], 100)
    
#     # Résolution de l'équation différentielle
#     solution = solve_ivp(dPdr, r_span, P0, t_eval=r_eval)
    
#     return solution.t, solution.y[0]
    

