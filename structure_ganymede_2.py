"""
Created on Fri Jan 24 19:47:44 2025

Programme Python 3 permettant de déterminer la structure interne d'un corps (de type Ganymède par ex) 
composé de trois couches de densités constantes.

n = noyau
m = manteau
e = couche externe

@author: lou, Romain, Alexandre
"""

# Modules nécessaires
import numpy as np
import matplotlib.pyplot as plt 
import CoolProp.CoolProp as CP

plt.close('all')

# CONSTANTES
G = 6.674E-11  # SI

###############################################################################
### PARAMETRES A INITIALISER ###

# Rayon du corps (m) 
rayon = 2612.1E3  

# Masse du corps (kg) 
masse = 1.48E23

# Rayon du noyau (m)
rayon_n = 850E3

# Fraction du rayon occupée par le noyau
Rn_R = rayon_n / rayon

# Densité du noyau (kg/m³)  
densite_n = 8000

# Densité de la couche externe (kg/m³)
densite_e = 1000  

# Densité moyenne du corps (kg/m³)
densite_moy = 1940  

# Facteur d'inertie du corps 
inertie = 0.3115

###############################################################################
### PARAMETRES DE CALCUL ###


# Fraction du rayon du manteau sur le rayon total
Rm_R = np.linspace(Rn_R, 1.0, 1001)

# Initialisation des tableaux de résultats
densite_m_densite_moy = np.zeros(1000)
alpha = np.zeros(1000)

###############################################################################
### CALCUL ###
Rm_R = np.delete(Rm_R, 0)

for i, Rm in enumerate(Rm_R):
    
    numerateur = (1 - (densite_n / densite_moy) * Rn_R**3 - (densite_e / densite_moy) * (1 - Rm**3))
    denominateur = Rm**3 - Rn_R**3
    densite_m_densite_moy[i] = numerateur/denominateur
    
    # Calcul du facteur d'inertie C/MR²
    alpha[i] = (2/5) * (densite_e * (1 - Rm**5) + densite_m_densite_moy[i] * densite_moy * (Rm**5 - Rn_R**5) + densite_n * Rn_R**5) / densite_moy


###############################################################################
### CALCUL GAMME DE DONNEES RECEVABLE ###

deltaI = np.zeros(len(Rm_R))

for i in range (len(Rm_R)) :
    deltaI[i] = np.fabs(inertie - alpha[i])

indice = np.argmin(deltaI)
print("L'indice de alpha qui vérifie au plus proche le facteur d'inertie est :", indice, "\n")
    
rayon_m = Rm_R[indice] * rayon
densite_m = densite_m_densite_moy[indice] * densite_moy

print("Avec les paramètres renseignés, le rayon du manteau doit être de {:.2f} km. Sa densité est de {:.2f} kg.m^-3. \n".format(rayon_m*1E-3, densite_m))



# PARAMÈTRES DE GANYMEDE
pourcentage_glace = 40 / 100
pourcentage_silicate = 60 / 100

# ECHANTILLONAGE DU RAYON
ech_rayon = np.linspace(0, rayon, 1000)

# DEFINITION DU PAS D'INTEGRATION
h = ech_rayon[1]


                                     #############################################
                                     #        STRUCTURE INTERNE DU CORPS         #
                                     #############################################

# DEFINITION DU PROFIL DE DENSITE DU CORPS
def densite(r) :
    if r <= rayon_n :
        return densite_n
    elif rayon_n < r <= rayon_m :
        return densite_m
    else :
        return densite_e


# CALCUL DE LA MASSE DU CORPS ET DU RAYON/MASSE DES COUCHES
def structure_corps() :
    
    # Masse du corps
    volume_corps = 4 / 3 * np.pi * rayon**3 
    masse_corps = densite_moy * volume_corps
    
    # Masse du manteau
    masse_m = densite_m * 4 * np.pi / 3 * (rayon_m**3 - rayon_n**3)
    
    # Rayon et masse de la couche de glace 
    masse_e = densite_e * 4 * np.pi / 3 * (rayon**3 - rayon_m**3)
    epaisseur_e = rayon - rayon_n
    
    return masse_corps, masse_m, epaisseur_e, masse_e

masse_corps, masse_m, epaisseur_e, masse_e = structure_corps()

print("La masse du corps sous nos hypothèses est de {:.2e} kg.".format(masse_corps))
print("Ecart avec la réalité : {:.2f} %.".format(abs((masse - masse_corps) / masse_corps) * 100), "\n")



                                   #############################################
                                   #  GRAVITE / PRESSION / FACTEUR D'INERTIE   #
                                   #############################################

# CALCUL DU PROFIL DE MASSE, GRAVITE, PRESSION
def calcul_profils() :
    
    # PROFIL DE MASSE 
    def dmdr(r) :
        return 4 * np.pi * densite(r) * r**2
    
    profil_masse = np.zeros(len(ech_rayon))
    for i in range (len(ech_rayon)) :
        profil_masse[i] = profil_masse[i-1] + h * dmdr(ech_rayon[i])
        
    # PROFIL DE GRAVITE OBTENU NUMERIQUEMENT 
    def dgdr(r, g) :
        return 4 * np.pi * G * densite(r)  - (2 * g) / r
    
    g = np.zeros(len(ech_rayon))
    for i in range(1, len(ech_rayon)) :
        g[i] = g[i - 1] + h * dgdr(ech_rayon[i], g[i - 1])
       
    # PROFIL DE GRAVITE OBTENU ANALYTIQUEMENT (pour comparer avec le profil numérique)
    g_analytique = np.zeros(len(ech_rayon))
    for i in range(1, len(ech_rayon)) :
        g_analytique[i] = G * profil_masse[i] / ech_rayon[i]**2
        
    # PROFIL DE PRESSION
    def dPdr(r, g) :
        return - densite(r) * g
    
    P = np.zeros(len(ech_rayon))
    for i in range (len(ech_rayon)) :
        P[i] = P[i-1] + h * dPdr(ech_rayon[i], g[i])
   
    # Réajustement du profil grâce à la pression de surface
    for i in range(len(ech_rayon)) :
        P[i] = P[i] - P[len(ech_rayon)-1]
        
    return profil_masse, g, g_analytique, P

profil_masse, profil_gravite, profil_gravite_analytique, profil_pression = calcul_profils()


# TRACE DU PROFIL DE GRAVITE
plt.figure(3)
plt.plot(ech_rayon*1E-3, profil_gravite_analytique, color='red', label="solution analytique")
plt.plot(ech_rayon*1E-3, profil_gravite, color='blue', label="solution numérique")
plt.xlabel("Rayon (km)")
plt.ylabel("g (m/s²)")
plt.title("Variation de l'accélération de la pesanteur dans Ganymède")
plt.legend()
plt.grid()

# TRACE DU PROFIL DE PRESSION
plt.figure(4)
plt.plot(ech_rayon*1E-3, profil_pression)
plt.xlabel("Rayon (km)")
plt.ylabel("Pression (Pa)")
plt.title("Variation de la pression dans Ganymède")
plt.grid()


def calcul_facteur_inertie() :
    
    # CALCUL DU MOMENT D'INERTIE
    def dIdr(r) :
        return 8/3 * np.pi * r**4 * densite(r)
    
    moment_inertie = np.zeros(len(ech_rayon))
    
    for i in range (len(ech_rayon)) :
        moment_inertie[i] = moment_inertie[i-1] + ech_rayon[1] * dIdr(ech_rayon[i])
    
    # FACTEUR D'INERTIE
    alpha = moment_inertie[len(moment_inertie)-1] / (masse * rayon**2)
    
    return alpha

facteur_inertie = calcul_facteur_inertie()

print("Le facteur d'inertie sous nos hypothèses est de {:.5}.".format(facteur_inertie))
print("Écart avec la réalité : {:.2f} %.".format(abs((inertie - facteur_inertie) / facteur_inertie) * 100), "\n")

                                 #############################################
                                 #                TEMPERATURE               #
                                 #############################################

# # CALCUL DU FLUX DE TEMPERATURE
# prod_radio = 2E-12 # W/kg
# k_silicate = 0.22 # W/m/kg
# k_glace = 2.1 # W/m/kg
# flux_surf = 0.002 # W/m²
# T_surf = 110 # K


# def source_chaleur(r) :
#     if r <= rayon_silicate :
#         return prod_radio
#     else :
#         return 0

# def conductivite_th(r) :
#     if r <= rayon_silicate :
#         return k_silicate
#     else :
#         return k_glace
# """
# a mettre au dessus 
# """
    
# def calcul_profil_temperature() :
    
#     # FLUX DE CHALEUR 
#     def dqdr(q, r) :
#         return densite(r) * source_chaleur(r) - 2 * q /r
    
#     flux = np.zeros(len(ech_rayon))
    
#     # Résolution du flux de chaleur grâce à la méthode d'Euler
#     for i in range (1, len(ech_rayon)) :
#         flux[i] = flux[i-1] + h * dqdr(flux[i-1], ech_rayon[i])
    
#     # PROFIL DE TEMPERATURE
#     def dTdr(q, r) :
#         return - q / conductivite_th(r)

#     T = np.zeros(len(ech_rayon))
    
#     # Résolution du profil de température grâce à la méthode d'Euler
#     for i in range(len(ech_rayon)) :
#         T[i] = T[i-1] + h * dTdr(flux[i], ech_rayon[i])
    
#     # Réajustement du profil grâce à la température de surface
#     for i in range(len(ech_rayon)) :
#         T[i] = T[i] - T[len(ech_rayon)-1] + T_surf
        
    
#     # SOLUTION ANALYTIQUE DU PROFIL DE TEMPERATURE (pour comparaison)
#     T_analytique = np.zeros(len(ech_rayon))
#     # for i in range(len(ech_rayon)):
#     #     if ech_rayon[i] <= rayon_silicate:
#     #         T_analytique[i] = T_surf + source_chaleur / (6*k_silicate) * (rayon**2 - ech_rayon[i]**2)
#     #     else:
#     #         T_analytique[i] = T_surf + rayon * flux_surf / k_glace *(rayon/ech_rayon[i] - 1)

#     return flux, T, T_analytique

# flux, profil_temperature, profil_temperature_analytique = calcul_profil_temperature()

# plt.figure(5)
# # plt.plot(ech_rayon*10**(-3),profil_temperature_analytique, color='red', label="solution analytique")
# plt.plot(ech_rayon*1E-3, profil_temperature, color='blue', label="solution numerique")
# plt.xlabel("Rayon (en km)")
# plt.ylabel("Température (K)")
# plt.title("Variation de la température dans Ganymède")
# plt.legend()
# plt.grid()

# plt.figure(55)
# # plt.plot(ech_rayon*10**(-3),profil_temperature_analytique, color='red', label="solution analytique")
# plt.plot(ech_rayon*1E-3, flux, color='blue', label="solution numerique")
# plt.xlabel("Rayon (en km)")
# plt.ylabel("Flux (W/m²)")
# plt.title("Flux de température dans Ganymède")
# plt.legend()
# plt.grid()




# # FLUX DE CHALEUR 
# def dqdr(q, r) :
#     return densite(r) * source_chaleur(r) - 2 * q /r

# # PROFIL DE TEMPERATURE
# def dTdr(q, r) :
#     return - q / conductivite_th(r)

# flux = np.zeros(len(ech_rayon))
# T = np.zeros(len(ech_rayon))
# T[0] = 200000

# Water = CP.AbstractState("HEOS", "Water")
# state = []

# T_test = [1]
# T_test[0] = 300000
# T_test.append(T[0])
# j = 1

# while abs(T_test[j]-T_test[j-1]) >= 100 :
#     state = []
#     for i in range (0, len(ech_rayon)-1) :
        
#         # Dans la couche de glace/eau 
#         if ech_rayon[i] >=   :
            
#             # Si dans la glace : calcul normal du flux
#             if T[i] < Water.melting_line(CP.iT, CP.iP, profil_pression[i]):
#                 flux[i+1] = flux[i] + h * dqdr(flux[i], ech_rayon[i+1])
#                 T[i+1] = T[i] + h * dTdr(flux[i+1], ech_rayon[i+1])
                
#                 state.append(1) #solide
                
#             # Si dans l'eau : T=cte
#             else : 
#                 T[i+1] = T[i]
                
#                 state.append(0) #liquide
                
#         # Dans silicate : calcul normal de T
#         else :
#             flux[i+1] = flux[i] + h * dqdr(flux[i], ech_rayon[i+1])
#             T[i+1] = T[i] + h * dTdr(flux[i+1], ech_rayon[i+1])
            
#             state.append(-1)
            
#     # Réajustement du profil grâce à la température de surface
#     for i in range(len(ech_rayon)) :
#         T[i] = T[i] - T[len(ech_rayon)-1] + T_surf
    
#     state.append(1) # glace à la surface 
#     T_test.append(T[0])
#     j += 1
#     print(j)
    
    
        



# plt.figure(6)
# plt.plot(ech_rayon*1E-3, state, color='blue')

# plt.xlabel("Rayon (en km)")
# plt.ylabel("Etat : silicate=-1 ; eau_liq=0 ; glace=1")
# plt.title("Variation de la phase de l'eau dans Ganymède")
# plt.grid()






