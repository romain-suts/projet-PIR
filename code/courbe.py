import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data.txt')

masse_silicate = data[:,0]
rayon_metallique = data[:,1]
rayon_silicate = data[:,2]
z = np.polyfit(masse_silicate,rayon_metallique, 3)
modele_metal = [z[3] + z[2] * val + z[1] * val**2 + z[0]*val**3 for val in masse_silicate]

fit = np.polyfit(masse_silicate,rayon_silicate, 3)
modele_silicate = [fit[3] + fit[2] * val + fit[1] * val**2 + fit[0]*val**3 for val in masse_silicate]

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot(masse_silicate, modele_silicate, color="red")
ax2.plot(masse_silicate, modele_metal, color="blue")

ax1.set_xlabel("masse volumique du silicate (kg/m3)")
ax1.set_ylabel("rayon du noyau silicaté (m)",color = "red")
ax1.tick_params(axis="y", labelcolor="red")

ax2.set_ylabel("rayon du noyau métallique (m)",color="blue")
ax2.tick_params(axis="y",labelcolor = "blue")

fig.suptitle("rayon des differents noyaux en fonction de la masse volumique du silicate")
ax1.grid()
fig.tight_layout()
plt.show()