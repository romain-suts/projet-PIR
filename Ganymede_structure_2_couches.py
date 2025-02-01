# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 20:41:52 2025

@author: louau
"""

# IMPORTATION DES MODULES
from math import *
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# CONSTANTES
densite_glace = 917  # kg/m^3
densite_silicate = 3500  # kg/m^3
G = 6.674 * 10**(-11)  # SI

def densite(r):
    if r <= rayon_silicate:
        return densite_silicate
    else:
        return densite_glace

# PARAMÈTRES DE GANYMEDE
pourcentage_glace = 40 / 100
pourcentage_silicate = 60 / 100
rayon = 2612.1 * 10**3  # m
masse = 1.48 * 10**23  # kg


# CALCUL DE LA MASSE DU CORPS
def calcul_masse():
    densite_corps = pourcentage_glace * densite_glace + pourcentage_silicate * densite_silicate
    volume_corps = 4 / 3 * np.pi * rayon**3 
    masse_corps = densite_corps * volume_corps
    return masse_corps, densite_corps

masse_corps, densite_corps = calcul_masse()
print("La masse de Ganymède sous nos hypothèses est de {:.2e} kg.".format(masse_corps))
print("Écart avec la réalité : {:.2f} %.".format(abs((masse - masse_corps) / masse_corps) * 100))


# CALCUL DES RAYONS DES COUCHES
def calcul_rayon():
    # SILICATE
    masse_silicate = masse_corps * pourcentage_silicate
    rayon_silicate = ((3 * masse_silicate) / (4 * np.pi * densite_silicate))**(1/3)
    # GLACE
    masse_glace = masse_corps * pourcentage_glace
    rayon_glace = rayon - rayon_silicate
    return rayon_silicate, masse_silicate, rayon_glace, masse_glace

rayon_silicate, masse_silicate, rayon_glace, masse_glace = calcul_rayon()
print("Rayon du noyau de silicate : {:.2f} km, soit {:.2f} % du rayon total.".format(rayon_silicate * 10**(-3), (rayon_silicate / rayon) * 100))
print("Masse du noyau de silicate : {:.2e} kg.".format(masse_silicate))

#%%

# DEFINITION METHODE D'INTEGRATION
def euler_method(f, y0, t0, t_end, h):
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0

    while t < t_end:
        y += h * f(t, y)
        t += h
        t_values.append(t)
        y_values.append(y)
    
    return t_values, y_values

#%%

# PROFIL DE MASSE
distr_rayon = np.linspace(0, rayon, 1000)

def dmdr(r, m):
    return 4 * np.pi * densite(r) * r**2

distr_masse = euler_method(lambda r, m: dmdr(r, m), 0, 0, rayon, rayon / 1000)[1]
distr_masse = distr_masse[:1000] # la résolution ajout 1 élement en trop dans ce tableau

plt.figure(1)
plt.plot(distr_rayon*10**(-3), distr_masse, label="Distribution de masse")
plt.xlabel("Rayon (km)")
plt.ylabel("Masse (kg)")
plt.title("Distribution de masse dans Ganymède")
plt.grid()

#%% 

# PROFIL DE DENSITE
distr_densite = np.zeros(1000)
for i in range(len(distr_rayon)):
    distr_densite[i] = densite(distr_rayon[i])
    
plt.figure(2)
plt.plot(distr_rayon*10**(-3), distr_densite)
plt.xlabel("Rayon (km)")
plt.ylabel("Masse volumique (kg/m^3)")
plt.title("Profil de densité dans Ganymède")
plt.grid()

#%% 

# PROFIL DE g 
def calcul_gravite_analytique():
    g_analytique = np.zeros(1000)
    for i in range(1000):
        g_analytique[i] = G*distr_masse[i]/distr_rayon[i]**2
    return g_analytique

g_analytique = calcul_gravite_analytique()

# CALCUL DE g(r) AVEC MÉTHODE D'EULER
g = np.zeros(1000)
for i in range(1, len(distr_rayon)):
    r = distr_rayon[i]
    m = distr_masse[i]
    densite = distr_densite[i-1]
    h = rayon/1000
    g_i = g[i - 1]
    
    # Calcul de dg/dr à l'étape i
    dgdr = 4 * np.pi * G * densite  - (2 * g_i) / r
    
    # Mise à jour de g(r) en utilisant la méthode d'Euler
    g[i] = g[i - 1] + h * dgdr


plt.figure(3)
plt.plot(distr_rayon*10**(-3),g_analytique,color='red',label="solution analytique")
plt.plot(distr_rayon*10**(-3),g,color='blue',label="solution numérique")
plt.xlabel("Rayon (en km)")
plt.ylabel("g (m/s²)")
plt.title("Variation de l'accélération de la pesanteur dans Ganymède")
plt.legend()
plt.grid()

#%%

# PROFIL DE PRESSION dP/dr
def dPdr(r, P):
    densite_r = distr_densite[int(r / (rayon / len(distr_densite)))]  # Extraction de la densité pour le rayon r
    g_r = g[int(r / (rayon / len(g)))]  # Extraction de la gravité pour le rayon r
    return -densite_r * g_r

# Condition initiale : pression centrale
P0 = 7.57e9  #Pa

# Exemple d'utilisation de la méthode d'Euler
distr_pression = euler_method(lambda r, P: dPdr(r, P), P0, 0, rayon, rayon / 1000)[1]
distr_pression = distr_pression[1:]


plt.figure(4)
plt.plot(distr_rayon*10**(-3), distr_pression)
plt.xlabel("Rayon (km)")
plt.ylabel("Pression (Pa)")
plt.title("Profil de la pression dans Ganymède")
plt.grid()
plt.show()

#%%

# CALCUL DU MOMENT D'INERTIE

def dIdr(r):
    def densite(r):
        if r <= rayon_silicate:
            return densite_silicate
        else:
            return densite_glace
    rho = densite(r)
    return 8/3*np.pi*r**4*rho

I = np.zeros(len(distr_rayon))

for i in range (len(distr_rayon)):
    I[i] = I[i-1] + distr_rayon[1]*dIdr(distr_rayon[i])

alpha = I[len(I)-1]/(masse*rayon**2)

plt.figure(5)
plt.plot(distr_rayon, I)
















