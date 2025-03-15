# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 14:03:00 2025

@author: louau

Modèle à 2 couches - Moment inertie et choix densité
"""

# Modules nécessaires
import numpy as np
import matplotlib.pyplot as plt 
plt.close('all')

Ri_Re = np.linspace(0.,1.,100)
# pi_pe = [1., 2., 3.4, 4., 5., 6.]
pi_pe = [3.4]

# color = 'tab:red'

fig, ax1 = plt.subplots()

for rho in pi_pe :
    alpha = 2 / 5 * (((rho - 1) * Ri_Re**5 + 1) / (((rho - 1) * Ri_Re**3 + 1)))
    ax1.plot(Ri_Re, alpha, label=r'$\rho_i/\rho_e = 3.4$', color='green')

ax1.set_xlabel('$Ri / Re$')
ax1.set_ylabel('$C / MR²$')
ax1.set_xlim(0, 1)
ax1.set_ylim(0.2, 0.41)
ax1.legend(title=r"Facteur d'inertie", loc='upper left')


ax2 = ax1.twinx()  

for rho in pi_pe :
    pm_pe = (rho - 1) * Ri_Re**3 + 1
    ax2.plot(Ri_Re, pm_pe, linestyle='dashed', color='green', label=r'$\rho_i/\rho_e = 3.4$')


ax2.set_ylabel('$pm/pe$') 
ax2.set_ylim(1, 6)
ax2.legend(title=r'$\rho_m/\rho_e$', loc='lower left')

fig.tight_layout()
ax1.plot([0,1], [0.3115, 0.3115], color='tab:pink')
ax1.plot([0.60943, 0.60943], [0.2, 0.3115], color='black', linestyle='-.')
ax1.plot([0.74622, 0.74622], [0.2, 0.3115], color='black', linestyle='-.')
ax2.plot([0.60943, 1], [1.543332, 1.54332], color='black', linestyle='-.')
ax2.plot([0.74622, 1], [1.997327274625, 1.997327274625], color='black', linestyle='-.')

plt.show()

#%%

inertie = 0.3115

pe_pm = np.linspace(0., 1., 100)

Ri_Re = ((5/2 * inertie - pe_pm) / (1 - pe_pm))**(1/2)

pi_pm = (1 - pe_pm)*(Ri_Re)**(-3) + pe_pm 


# Création de la figure
fig1, ax1 = plt.subplots()

# Tracé de pi/pm en fonction de pe/pm
ax1.plot(pe_pm, pi_pm, label='$p_i / p_m$', color='tab:blue')
ax1.set_xlabel('$p_e / p_m$')
ax1.set_ylabel('$p_i / p_m$', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')
ax1.set_ylim(1,3)
ax1.legend(loc='lower left')

# Deuxième axe pour Ri/Re
ax2 = ax1.twinx()  
ax2.plot(pe_pm, Ri_Re, linestyle='dashed', color='tab:red', label='$R_i / R_e$')
ax2.set_ylabel('$R_i / R_e$', color='tab:red') 
ax2.tick_params(axis='y', labelcolor='tab:red')
ax2.legend(loc='lower right')
ax2.set_xlim(0,1)
# ax2.plot([0.501,0.501],[0.1,0.9], color='black', linestyle='-.')
# ax1.plot([0,0.501],[1.703,1.703], color='black', linestyle='-.')
# ax2.plot([0.501,1],[0.746,0.746], color='black', linestyle='-.')


# Affichage de la figure
fig1.tight_layout()
plt.show()





