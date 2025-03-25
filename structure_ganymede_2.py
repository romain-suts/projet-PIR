"""
Created on Fri Jan 24 19:47:44 2025

Programme Python 3 permettant de déterminer la structure interne d'un corps (de type Ganymède par ex) 
composé de deux couches de densités constantes.

@author: lou, Romain, Alexandre
"""

# Modules nécessaires
import numpy as np
import matplotlib.pyplot as plt 
import CoolProp.CoolProp as CP

plt.close('all')


# CONSTANTES
densite_glace = 917  # kg/m^3
densite_silicate = 3500  # kg/m^3
G = 6.674E-11  # SI

# PARAMÈTRES DE GANYMEDE
pourcentage_glace = 40 / 100
pourcentage_silicate = 60 / 100
rayon = 2612.1E3  # m
masse = 1.48E23  # kg

# ECHANTILLONAGE DU RAYON
distr_rayon = np.linspace(0, rayon, 1000)

# DEFINITION DU PAS D'INTEGRATION
h = distr_rayon[1]


                                     #############################################
                                     #      STRUCTURE INTERNE DE GANYMEDE       #
                                     #############################################

# DEFINITION DU PROFIL DE DENSITE DE GANYMEDE
def densite(r) :
    if r <= rayon_silicate :
        return densite_silicate
    else:
        return densite_glace


# CALCUL DE LA MASSE DU CORPS ET DU RAYON/MASSE DES COUCHES
def structure_corps() :
    
    # Masse du corps
    densite_corps = pourcentage_glace * densite_glace + pourcentage_silicate * densite_silicate
    volume_corps = 4 / 3 * np.pi * rayon**3 
    masse_corps = densite_corps * volume_corps
    
    # Rayon et masse de la couche de silicate 
    masse_silicate = masse_corps * pourcentage_silicate
    rayon_silicate = ((3 * masse_silicate) / (4 * np.pi * densite_silicate))**(1/3)
    
    # Rayon et masse de la couche de glace 
    masse_glace = masse_corps * pourcentage_glace
    epaisseur_glace = rayon - rayon_silicate
    
    return masse_corps, densite_corps, rayon_silicate, masse_silicate, epaisseur_glace, masse_glace

masse_corps, densite_corps, rayon_silicate, masse_silicate, epaisseur_glace, masse_glace = structure_corps()

print("La masse de Ganymède sous nos hypothèses est de {:.2e} kg.".format(masse_corps))
print("Écart avec la réalité : {:.2f} %.".format(abs((masse - masse_corps) / masse_corps) * 100))

print("Rayon du noyau de silicate : {:.2f} km, soit {:.2f} % du rayon total.".format(rayon_silicate * 1E-3, (rayon_silicate / rayon) * 100))
print("Masse du noyau de silicate : {:.2e} kg.".format(masse_silicate))



                                   #############################################
                                   #  GRAVITE / PRESSION / FACTEUR D'INERTIE   #
                                   #############################################

# CALCUL DU PROFIL DE MASSE, GRAVITE, PRESSION
def calcul_profil_gravite() :
    
    # PROFIL DE MASSE 
    def dmdr(r) :
        return 4 * np.pi * densite(r) * r**2
    
    profil_masse = np.zeros(len(distr_rayon))
    for i in range (len(distr_rayon)) :
        profil_masse[i] = profil_masse[i-1] + h * dmdr(distr_rayon[i])
        
    # PROFIL DE GRAVITE OBTENU NUMERIQUEMENT 
    def dgdr(r, g) :
        return 4 * np.pi * G * densite(r)  - (2 * g) / r
    
    g = np.zeros(len(distr_rayon))
    for i in range(1, len(distr_rayon)) :
        g[i] = g[i - 1] + h * dgdr(distr_rayon[i], g[i - 1])
       
    # PROFIL DE GRAVITE OBTENU ANALYTIQUEMENT (pour comparer avec g)
    g_analytique = np.zeros(len(distr_rayon))
    for i in range(1, len(distr_rayon)) :
        g_analytique[i] = G * profil_masse[i] / distr_rayon[i]**2
        
    # PROFIL DE PRESSION
    def dPdr(r, grav) :
        return - densite(r) * grav
    
    P = np.zeros(len(distr_rayon))
    for i in range (len(distr_rayon)) :
        P[i] = P[i-1] + h * dPdr(distr_rayon[i], g[i])
   
    # Réajustement du profil grâce à la température de surface
    for i in range(len(distr_rayon)) :
        P[i] = P[i] - P[len(distr_rayon)-1]
        
    return profil_masse, g, g_analytique, P

profil_masse, profil_gravite, profil_gravite_analytique, profil_pression = calcul_profil_gravite()


# TRACE DU PROFIL DE GRAVITE DANS GANYMEDE
plt.figure(3)
plt.plot(distr_rayon*1E-3, profil_gravite_analytique, color='red', label="solution analytique")
plt.plot(distr_rayon*1E-3, profil_gravite, color='blue', label="solution numérique")
plt.xlabel("Rayon (en km)")
plt.ylabel("g (m/s²)")
plt.title("Variation de l'accélération de la pesanteur dans Ganymède")
plt.legend()
plt.grid()

# TRACE DU PROFIL DE PRESSION
plt.figure(4)
plt.plot(distr_rayon*1E-3, profil_pression)
plt.xlabel("Rayon (en km)")
plt.ylabel("Pression (Pa)")
plt.title("Variation de la pression dans Ganymède")
plt.grid()


def calcul_facteur_inertie() :
    
    # CALCUL DU MOMENT D'INERTIE
    def dIdr(r) :
        return 8/3 * np.pi * r**4 * densite(r)
    
    moment_inertie = np.zeros(len(distr_rayon))
    
    for i in range (len(distr_rayon)) :
        moment_inertie[i] = moment_inertie[i-1] + distr_rayon[1] * dIdr(distr_rayon[i])
    
    # FACTEUR D'INERTIE
    alpha = moment_inertie[len(moment_inertie)-1] / (masse * rayon**2)
    
    return alpha

facteur_inertie = calcul_facteur_inertie()


                                 #############################################
                                 #                TEMPERATURE               #
                                 #############################################

# CALCUL DU FLUX DE TEMPERATURE
prod_radio = 2E-12 # W/kg
k_silicate = 0.22 # W/m/kg
k_glace = 2.1 # W/m/kg
flux_surf = 0.002 # W/m²
T_surf = 110 # K


def source_chaleur(r) :
    if r <= rayon_silicate :
        return prod_radio
    else :
        return 0

def conductivite_th(r) :
    if r <= rayon_silicate :
        return k_silicate
    else :
        return k_glace
"""
a mettre au dessus 
"""
    
def calcul_profil_temperature() :
    
    # FLUX DE CHALEUR 
    def dqdr(q, r) :
        return densite(r) * source_chaleur(r) - 2 * q /r
    
    flux = np.zeros(len(distr_rayon))
    
    # Résolution du flux de chaleur grâce à la méthode d'Euler
    for i in range (1, len(distr_rayon)) :
        flux[i] = flux[i-1] + h * dqdr(flux[i-1], distr_rayon[i])
    
    # PROFIL DE TEMPERATURE
    def dTdr(q, r) :
        return - q / conductivite_th(r)

    T = np.zeros(len(distr_rayon))
    
    # Résolution du profil de température grâce à la méthode d'Euler
    for i in range(len(distr_rayon)) :
        T[i] = T[i-1] + h * dTdr(flux[i], distr_rayon[i])
    
    # Réajustement du profil grâce à la température de surface
    for i in range(len(distr_rayon)) :
        T[i] = T[i] - T[len(distr_rayon)-1] + T_surf
        
    
    # SOLUTION ANALYTIQUE DU PROFIL DE TEMPERATURE (pour comparaison)
    T_analytique = np.zeros(len(distr_rayon))
    # for i in range(len(distr_rayon)):
    #     if distr_rayon[i] <= rayon_silicate:
    #         T_analytique[i] = T_surf + source_chaleur / (6*k_silicate) * (rayon**2 - distr_rayon[i]**2)
    #     else:
    #         T_analytique[i] = T_surf + rayon * flux_surf / k_glace *(rayon/distr_rayon[i] - 1)

    return T, T_analytique

profil_temperature, profil_temperature_analytique = calcul_profil_temperature()

plt.figure(5)
# plt.plot(distr_rayon*10**(-3),profil_temperature_analytique, color='red', label="solution analytique")
plt.plot(distr_rayon*1E-3, profil_temperature, color='blue', label="solution numerique")
plt.xlabel("Rayon (en km)")
plt.ylabel("Température (K)")
plt.title("Variation de la température dans Ganymède")
plt.legend()
plt.grid()


# def phase_eau(r, P, T) :
#     if r <= rayon_silicate :
#         return np.nan
#     else :
#         return cp.PhaseSI('P', P, 'T', T, 'Water')

# def etat_couches() :
#     etat = []
    
#     # trad : solide=0 liq=1
#     trad = []
#     for i in range(len(distr_rayon)) : 
#         etat.append(phase_eau(distr_rayon[i], profil_pression[i], profil_temperature[i]))
#         if etat[i] == 'supercritical' :
#             trad.append(0)
#         elif etat[i] == 'supercritical_liquid' :
#             trad.append(1)
#         else : 
#             trad.append(-1)
        
            
#     return etat, trad

# etat, cbe_etat = etat_couches()

# plt.figure(6)
# plt.plot(distr_rayon*1E-3, cbe_etat, color='blue')
# plt.xlabel("Rayon (en km)")
# plt.ylabel("Etat")
# plt.title("Variation de la phase de l'eau dans Ganymède")
# plt.legend()
# plt.grid()




# FLUX DE CHALEUR 
def dqdr(q, r) :
    return densite(r) * source_chaleur(r) - 2 * q /r

# PROFIL DE TEMPERATURE
def dTdr(q, r) :
    return - q / conductivite_th(r)

flux = np.zeros(len(distr_rayon))
T = np.zeros(len(distr_rayon))
T[0] = 200000

Water = CP.AbstractState("HEOS", "Water")
state = []

T_test = [1]
T_test[0] = 300000
T_test.append(T[0])
j = 1

while abs(T_test[j]-T_test[j-1]) >= 100 :
    state = []
    for i in range (0, len(distr_rayon)-1) :
        
        # Dans la couche de glace/eau 
        if distr_rayon[i] >= rayon_silicate :
            
            # Si dans la glace : calcul normal du flux
            if T[i] < Water.melting_line(CP.iT, CP.iP, profil_pression[i]):
                flux[i+1] = flux[i] + h * dqdr(flux[i], distr_rayon[i+1])
                T[i+1] = T[i] + h * dTdr(flux[i+1], distr_rayon[i+1])
                
                state.append(1) #solide
                
            # Si dans l'eau : T=cte
            else : 
                T[i+1] = T[i]
                
                state.append(0) #liquide
                
        # Dans silicate : calcul normal de T
        else :
            flux[i+1] = flux[i] + h * dqdr(flux[i], distr_rayon[i+1])
            T[i+1] = T[i] + h * dTdr(flux[i+1], distr_rayon[i+1])
            
            state.append(-1)
            
    # Réajustement du profil grâce à la température de surface
    for i in range(len(distr_rayon)) :
        T[i] = T[i] - T[len(distr_rayon)-1] + T_surf
    
    state.append(1) # glace à la surface 
    T_test.append(T[0])
    j += 1
    print(j)
    
    
        



plt.figure(6)
plt.plot(distr_rayon*1E-3, state, color='blue')

plt.xlabel("Rayon (en km)")
plt.ylabel("Etat : silicate=-1 ; eau_liq=0 ; glace=1")
plt.title("Variation de la phase de l'eau dans Ganymède")
plt.grid()






