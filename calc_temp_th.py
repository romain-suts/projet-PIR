import math
from math import *
import numpy as np

rayon_lune = 2631.2E3
temperature_surface = 110.
chaleur_surface = 0.002
pas = 1000
conduc_ther_glace = 2.1
R = rayon_lune
fichier = open("temp_th.txt", "w")
rayon = rayon_lune
while (R>2*pas):
    fichier.write("{} {}\n".format(R,temperature_surface+(1/R-1/rayon_lune)*chaleur_surface*(rayon_lune**2)/conduc_ther_glace))
    R-=pas
fichier.close()