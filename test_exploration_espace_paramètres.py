# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 19:13:27 2025

@author: alexa
"""

# Modules nécessaires
import numpy as np
import matplotlib.pyplot as plt 

plt.close('all')

# CONSTANTES
G = 6.674E-11  # SI

###############################################################################
### PARAMETRES OBSERVER ###

# Rayon du corps (m) 
rayon = 2612.3E3  

# Masse du corps (kg) 
masse = 1.4819E23

# Densité moyenne du corps (kg/m³)
densite_moy = masse / (4/3 * np.pi * rayon**3)  
print("La densité moyenne du corps est : ", densite_moy,"\n")

# Facteur d'inertie du corps 
inertie = 0.3549


###############################################################################
# DEFINITION DES FONCTIONS 

# Echantillonnage du rayon 
ech_rayon = np.linspace(0, rayon, 1000)

# Pas d'intégration
h = ech_rayon[1]



# DEFINITION DU PROFIL DE DENSITE DU CORPS
def densite(r) :
    if r <= rayon_n :
        return densite_n
    elif rayon_n < r <= rayon_m :
        return densite_m
    else :
        return densite_e



# CALCUL DE LA MASSE DU CORPS ET DU RAYON/MASSE DES COUCHES
def calcul_masse() :
    
    def dmdr(r) :
        return 4 * np.pi * r**2 *densite(r)
    
    masse_corps = np.zeros(len(ech_rayon))
    for i in range (len(ech_rayon)):
        masse_corps[i] = masse_corps[i-1] + h * dmdr(ech_rayon[i])
    
    return masse_corps[len(ech_rayon)-1]




# CALCUL DU FACTEUR D'INERTIE
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



###############################################################################
# VARIABLES POUR L'EXPLORATION DES PARAMETRES

# Densité du noyau (kg/m³)  
densite_n = 0

# Rayon du noyau (m)
rayon_n = 0

# Densité de la couche externe (kg/m³)
densite_e = 1000


nb_point = 2
resultat = np.zeros((nb_point**2, 6))
tour = 0 

for densite_n in np.linspace(4000, 8000, nb_point) :
    for rayon_n in np.linspace(0, 800E3, nb_point) :
        
        # Fraction du rayon occupée par le noyau
        Rn_R = rayon_n / rayon
        
        # Fraction du rayon du manteau sur le rayon total
        Rm_R = np.linspace(Rn_R, 1.0, 1001)
        Rm_R = np.delete(Rm_R, 0)
        
        # Initialisation des tableaux de résultats
        densite_m_densite_moy = np.zeros(1000)
        alpha = np.zeros(1000)

        for i, Rm in enumerate(Rm_R) :
            
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
            
        # indice qui vérifie au plus proche le facteur d'inertie
        indice = np.argmin(deltaI)
        
        
        rayon_m = Rm_R[indice] * rayon
        densite_m = densite_m_densite_moy[indice] * densite_moy
        masse_corps = calcul_masse()
        
        facteur_inertie = calcul_facteur_inertie()
        
        resultat[tour, 2] = densite_m
        resultat[tour, 3] = rayon_m / rayon
        resultat[tour, 4] = ( masse_corps - masse) / masse  
        resultat[tour, 5] = (facteur_inertie - inertie) / inertie
        
        resultat[tour, 0] = densite_n
        
        # if 2500 < densite_m < 4000 :
        resultat[tour, 1] = rayon_n / rayon
        
        tour += 1







### FIGURE ###

# fig, ax1 = plt.subplots()

# # Tracé de C/MR² (courbe verte)
# ax1.plot(Rm_R, alpha, color='green', label=r'$C / MR^2$')
# ax1.axhline(inertie, color="red", label=fr"Moment inertie = {inertie}")
# ax1.set_xlabel(r'$R_m / R$', fontsize = 20)
# ax1.set_ylabel(r'$C / MR^2$', color='green', fontsize = 20)

# ax1.tick_params(axis='y', labelcolor='green')
# ax1.grid()

# # Ajout du second axe Y pour rho_s / rho_moyen
# ax2 = ax1.twinx()
# ax2.plot(Rm_R, densite_m_densite_moy, color='blue', label=r'ps_p')
# ax2.set_ylabel(r'$\rho_m / \rho_{\text{moyen}}$', color='blue', fontsize = 20)


# ax2.legend(loc='lower right',fontsize = 15)
# ax1.legend(loc='upper left', fontsize = 15)
# ax1.plot([Rm_R[indice], Rm_R[indice]],[0.2, 0.4], color='black', linestyle='--')
# plt.title('Contraintes pour la couche de silicates', fontsize = 20)

# plt.show()



                                     #############################################
                                     #        STRUCTURE INTERNE DU CORPS         #
                                     #############################################




print("Ce qui nous donne :\nRayon noyau = ", rayon_n*1E-3, " km \nRayon manteau = ", rayon_m*1E-3, "km \nEpaisseur de de la couche de glace = ", rayon*1E-3-rayon_m*1E-3, "km \n")


print("La masse du corps sous nos hypothèses est de {:.2e} kg.".format(masse_corps))
print("Ecart avec la réalité : {:.2f} %.".format(abs((masse - masse_corps) / masse_corps) * 100), "\n")



                                   #############################################
                                   #  GRAVITE / PRESSION / FACTEUR D'INERTIE   #
                                   #############################################





