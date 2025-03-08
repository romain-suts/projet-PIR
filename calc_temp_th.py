import math
from math import *
import numpy as np
import matplotlib.pyplot as plt
rayon_lune = 2631.2E3
temperature_surface = 110.
chaleur_surface = 0.002
pas = 1000
conduc_ther_glace = 2.1
R = [rayon_lune]
T = [temperature_surface]
fichier = open("temp_th.txt", "w")
rayon = rayon_lune
while (R[-1]>2*pas):
    fichier.write("{} {}\n".format(R[-1],T[-1]))
    R.append(R[-1]-pas)
    T.append(temperature_surface+(1/R[-1]-1/rayon_lune)*chaleur_surface*(rayon_lune**2)/conduc_ther_glace)
fichier.close()
plt.plot (R,T)
plt.xlabel("rayon")
plt.ylabel("T")
plt.show()