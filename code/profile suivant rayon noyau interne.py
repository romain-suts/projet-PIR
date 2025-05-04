# Modules nécessaires
import math
from math import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#fonction qui renvoie 1 si l'element est solide, 0 si il est liquide
def solide_(etat):
    if etat>=1:
        return 1
    elif (etat==-1):
        return 1
    elif (etat==-3):
        return 1 
    return 0
#fonction qui retourne la viscosité de la glace 
def ice_viscosity(T):
	#tau is the stress in MPa and T the temperature in K
	Tref=255
	R=8.314
	Ea=where(T<263,5.9e4,1.9e5)
	eta_ref=1.e15
	return eta_ref*np.exp(Ea/(R*Tref)*(Tref/T-1))
#fonction qui renvoie la chaleur produite par les condrite en W/kg 
def chaleur_radiactive (SSage):


	chondrite='CM'

	MeV=1.602e-13        #Joule
	nspy=365.25*24*3600. #number of seconds per a year
	Na=6.02214129e23     #Avogadro constant
	nX    =np.array([' 26Al',  ' 40K',  ' 60Fe', '232Th',  '235U',  '238U'])
	t_demi=np.array([ 7.17e5, 1.265e9, 2.620e6, 1.405e10, 7.038e8, 4.468e9])  #years
	E     =np.array([   3.12,  0.6087,    2.71,     39.0,    42.7,    46.1])  #Mev
	m_mole=np.array([  26.98,    39.1,   55.85,   323.04,  238.03,  238.03])  #g
	X_CI  =np.array([ 8.6e-3,  5.6e-4, 1.82e-1,   2.9e-8,  8.2e-9,  8.2e-9])/0.82  
	X_CM  =np.array([11.8e-3,  4.0e-4, 2.10e-1,   4.0e-8, 11.0e-9, 11.0e-9])/0.874  
	nXsXpd=np.array([     0., 1.18e-4,      0.,       1.,  7.2e-3,  0.9928])
	nXsXi =np.array([  5.e-5,      0.,   1.e-8,       0.,      0.,      0.])
	tau   =t_demi/log(2.)
	Xpd   =np.ones(6)
	Xi    =np.zeros(6)
	nXi   =np.zeros(6)
	nXpd  =np.zeros(6)
	Qpd   =np.zeros(6)
	Hpd   =np.zeros(6)
	Qi    =np.zeros(6)
	Hi    =np.zeros(6)



	for i in range(0,6):
		if nXsXpd[i]>0:
			nXpd[i]  = nXsXpd[i]*Xpd[i]
			nXi[i]   = nXpd[i]*exp(SSage/(tau[i]))
			Xi[i]    = Xpd[i]-nXpd[i]+nXi[i]
			nXsXi[i] = nXi[i]/Xi[i]
		else:
			Xi[i]    = Xpd[i]
			nXi[i]   = nXsXi[i]*Xi[i]


	Xi[4:6]=nXi[4]+nXi[5]
	nXsXi[4:6] = nXi[4:6]/Xi[4:6]


	print('{:<10}{:^20}{:^20}{:^20}{:^20}'.format(' nX', '[nX]/[X]pd','[X]pd', '[nX]/[X]i','[X]i'))
	for i in range(0,6):
		print( '{:<10}   {:<20.3e} {:<20.3e} {:<20.3e} {:<20.3e}'.format(nX[i], nXsXpd[i], Xpd[i], nXsXi[i], Xi[i]))

	print( '')
	print( '{:<10}{:^10}{:^30}{:^10}{:^30}{:^10}'.format(' nX','t_1/2', 'Qpd','Hpd', 'Qi','Hi'))
	for i in range(0,6):
		Qpd[i]=nXpd[i]*E[i]*MeV*Na/(m_mole[i]*1.e-3)
		Hpd[i]=Qpd[i]/(tau[i]*nspy)
		Qi[i]=nXi[i]*E[i]*MeV*Na/(m_mole[i]*1.e-3)
		Hi[i]=Qi[i]/(tau[i]*nspy)
		print( '{:<10}{:<20.3e}{:<20.3e}{:<20.3e}{:<20.3e}{:<20.3e}'.format(nX[i], t_demi[i], Qpd[i], Hpd[i], Qi[i], Hi[i]))

	nt=300

	time=SSage*exp(-linspace(13,0,nt))

	Xchond=X_CI
	if (chondrite == 'CM'):
		Xchond=X_CM


	H=np.zeros((7,nt))

	H[0,:]=Hi[0]*exp(-(time[:])/tau[0])*Xchond[0]
	for i in range(0,6):
		H[i,:]=H[i-1,:]+Hi[i]*Xchond[i]*exp(-(time[:])/tau[i])
		H[i,:]=         Hi[i]*Xchond[i]*exp(-(time[:])/tau[i])
		H[6,:]=H[6,:]+H[i,:]

	print("\n Initial total heat production in CM chondrites    : {:.3e} W/kg".format(H[-1,0]))
	print(" Present day total heat production in CM chondrites: {:.3e} W/kg\n".format(H[-1,-1]))

	return (H[-1,-1])

#fonction qui renvoie la capacite thermique de chaque element 
def capa_ther(etat):
    capa_ther_glace = 2060 #J/kg/K
    capa_ther_silicate = 44.4 #J/kg/K
    capa_ther_metal = 0.449 #J/kg/K
    if (etat >= 1):
        return capa_ther_glace
    elif (etat <= 3):
        return capa_ther_metal 
    return capa_ther_silicate
#fonction qui retourne la masse volumique de la couche 
def densite(R,rayon_noyau,rayon_metallique,etat,P,T,fraction_silicate_glace,masse_vol_silicate,masse_vol_metallique):
    #etablissement des constantes a des pressions et temperature de reference 
    
    #diferents masses volumiques de reference
    
    masse_vol_glace_I = 917 # kg/m^3
    masse_vol_glace_III = 1160 # kg/m^3
    masse_vol_glace_V = 1240 # kg/m^3
    masse_vol_glace_VI = 1310 # kg/m^3
    masse_vol_glace_VII = 1650 # kg/m^3
    masse_vol_eau = 1000 # kg/m^3
    
    #diferents compressibilité isotherme de reference
    
    compressibilite_isotherme_eau = 2.2*10**9 #Pa^-1
    compressibilite_isotherme_silicate = 2*10**(-11) #Pa^-1
    compressibilite_isotherme_glace_I = 12*10**(-12) #Pa^-1
    compressibilite_isotherme_glace_III = 5.5*10**(-12) #Pa^-1
    compressibilite_isotherme_glace_V = 4*10**(-12) #Pa^-1
    compressibilite_isotherme_glace_VI = 3.2*10**(-12) #Pa^-1
    compressibilite_isotherme_glace_VII =2.5*10**(-12) #Pa^-1
    
    #diferents expansions thermiques de reference
    
    expansion_thermique_glace_I = 9*10**(-5) #K^-1
    expansion_thermique_glace_III = 2*10**(-5) #K^-1
    expansion_thermique_glace_V = 3*10**(-5) #K^-1
    expansion_thermique_glace_VI = 1.5*10**(-5) #K^-1
    expansion_thermique_glace_VII = 2.5*10**(-6) #K^-1
    
    #si on est dans le noyau => silicate
    if (R<=rayon_metallique):
        
        return masse_vol_metallique
    
    #si on est dans le noyau => silicate
    
    if (R<=rayon_noyau):
        
        return masse_vol_silicate*(1+(P-101325)*compressibilite_isotherme_silicate+3*10**(-5)*T)
    
    #si on est dans le manteau glace => revoie la masse volumique de la glace/eau de la couche
    #les nombres renvoies au numéro de glace 
    else:
        
        masse_silicate_glace = fraction_silicate_glace*masse_vol_silicate*(1+(P-101325)*compressibilite_isotherme_silicate+3*10**(-5)*T) #masse de la fraction de silicate dans la glace
        
        if (etat==1):#faire attention aux temperatures de references pour les valeurs de compressibilite de d'expansion
            
            return (1-fraction_silicate_glace)*masse_vol_glace_I*(1+(P-101325)*compressibilite_isotherme_glace_I+expansion_thermique_glace_I*(T-273.15))+masse_silicate_glace
        
        elif (etat==3):
            
            return (1-fraction_silicate_glace)*masse_vol_glace_III*(1+(P-350*10**6)*compressibilite_isotherme_glace_III+expansion_thermique_glace_III*(T-273.15))+masse_silicate_glace
        
        elif (etat==5):
        
            return (1-fraction_silicate_glace)*masse_vol_glace_V*(1+(P-350*10**6)*compressibilite_isotherme_glace_V+expansion_thermique_glace_V*(T-273.15))+masse_silicate_glace
        
        elif (etat==6):
        
            return (1-fraction_silicate_glace)*masse_vol_glace_VI*(1+(P-0.6*10**9)*compressibilite_isotherme_glace_VI+expansion_thermique_glace_VI*(T-273.15))+masse_silicate_glace

        elif (etat==7):
        
            return (1-fraction_silicate_glace)*masse_vol_glace_VII*(1+(P-2.5*10**9)*compressibilite_isotherme_glace_VII+expansion_thermique_glace_VII*(T-(273.15+25)))+masse_silicate_glace
        
        #les conditions ne remplicent aucun type de glace donc on a de l'eau
        else:
        
            return masse_vol_eau*(1+(P-101325)/compressibilite_isotherme_eau+2.07*10**(-4)*(T-293.15))
        
#fonction qui renvoie l'emission radiogénique de la couche W/kg
def Source(r,rayon_noyau,emmission_radiogénique,effet_maree,masse_vol):#renvoi la puissance générée par kg de la couche
    
    #puissance emise par le silicate
    if(r<=rayon_noyau):
        
        return masse_vol*emmission_radiogénique+effet_maree
    
    #les autres elements ne retournent pas de chaleur de source radiogénique
    else: 
        
        return effet_maree

#renvoi la conduction thermique de la couche en fonction du materiau et de la temperature
def lambda_(R,rayon_noyau,rayon_metal,T,etat):
    
    #Constantes 
    conduc_ther_glace = 2.1
    conduc_ther_eau = 0.6
    conduc_ther_silicate = 0.22
    conduc_ther_metal = 80
    
    #si on se trouve dans le métal 
    if (R<=rayon_metal):
        
        return conduc_ther_metal
    
    #si on est dans le noyau => silicate
    if (R<=rayon_noyau):
    
        return conduc_ther_silicate
    
    
    else: 
    
        #si on est dans le manteau glace => possibilite d'avoir de l'eau
        if (etat == 0):
    
            return conduc_ther_eau
    
        else:
            
        #si on est dans le manteau glace => revoie la conductivité de la glace de la couche car ce n'est pas de l'eau
            return conduc_ther_glace-0.012*(T-273.15)

#fonction qui permet d'afficher les resultats
def affichage(tour,g,pression,chaleur,T,I,etat,masse_vol,rayon):#fonction qui affiche les graphiques
    x = np.linspace(rayon,0,tour+1) 
    plt.grid()
    plt.plot (x,g)
    plt.title("Acceleration de pesanteur en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("g(m/s²)")
    plt.show()
    plt.grid()
    plt.plot(x,pression)
    plt.title("Pression en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("Pression (Pa)")
    plt.show()
    plt.grid()
    plt.plot(x,T)
    plt.title("Température en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("Température (K)")
    plt.show()
    plt.grid()
    plt.plot(x,I)
    plt.title("Facteur d'inertie en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("moment d'inertie(kg.m²)")
    plt.show()
    plt.grid()
    plt.plot(x,etat)
    plt.title("Etat en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("etat")
    plt.show()
    plt.grid()
    plt.plot(T,pression)
    plt.title("Diagramme de Clapeyron")
    plt.xlabel("T (K)")
    plt.ylabel("P (Pa)")
    plt.yscale("log") 
    plt.show()
    plt.plot(x,masse_vol)
    plt.title("Masse volumique en fonction du rayon")
    plt.xlabel("Rayon(m)")
    plt.ylabel("masse volumique (kg/m^3)")
    plt.show()  
#fonction qui permet de calculer une premiere de forme de lune qu'avec de la glace en surfa ce 
def calc_forme_init(masse,rayon,rayon_metallique,densite_glace,densite_silicate,densite_metal): #fonction qui intialise la premiere forme de la lune pour des masses volumique constante (on cherche a retrouver la masse de la lune)
    
    #constante 

    volume =4/3*pi*np.power(rayon,3) #m^3
    densite_lune = masse/(volume) # kg/m^3
    #calcule du rayon du noyau par conservation de la masse
    rayon_noyau = ((rayon**3*(densite_lune-densite_glace)+rayon_metallique**3*(densite_silicate-densite_metal))/(densite_silicate-densite_glace))**(1/3)# m !! si noyau metallique pas bien calibre on peut avoir un rayon de silicate > rayon de la lune ex io
    return rayon_noyau,densite_lune 

#calcule la valeur de g vers le centre du noyau grace a la formule theorique
def calc_g(rayon,G,pas,masse_vol):
    #initialisation des constantes
    r=0 #km
    masse=0 #kg
    
    # calcule de la masse a un certain rayon
    while (r<=rayon):
        masse += pas*4*pi*np.power(r,2)*masse_vol
        r += pas
        
    #retourne la valeur theorique a partir de la masse trouvee
    return (G*masse/np.power(rayon,2))

#determine l'etat de la couche 
def etat_l(rayon,rayon_noyau,rayon_metallique,P,T):
    
    #constantes
    delta_V_ice_I = 1/densite(rayon_noyau+1,rayon_noyau,rayon_metallique,1,P,T,0,0,0)-1/densite(rayon_noyau+1,rayon_noyau,rayon_metallique,0,P,T,0,0,0) #m^3
    chaleur_latente_I = 333.55*1000 #J/kg
    P_atm = 101325 #Pa
    temp_fusion_metal =  1553+273.15 #K
    if (rayon<=rayon_metallique):
        
        if (T<temp_fusion_metal):#si le metal est encore solide
            return -3
        
        else:#sinon il est liquide
            return -4
    
    #si on est dans le manteau glacé on cherche le type de glace que l'on a 
    elif (rayon>rayon_noyau):
        
        if (T<(math. log(T)<math. log(273.15)+(delta_V_ice_I/chaleur_latente_I)*(P-P_atm)) and (P<209.9*10**6)):#glace 1 (math. log(T)<math. log(273.15)+(delta_V_ice_I/chaleur_latente_I)*(P-P_atm))
            return 1
        
        elif ((P*10**(-6)>209.5+101.1*((T/251.15)**(42.86)-1)) and (P<350.1*10**6) and (P>=209.9*10**6)):#glace 3
            return 3
        
        elif((P*10**(-6)>355+373.6*((T/256.43)**(8.66)-1)) and (P<632.4*10**6) and (P>=350.1*10**6)):#glace 5
            return 5
        
        elif((P*10**(-6)>618.4+661.4*((T/272.73)**(4.69)-1)) and (P<2.216*10**9) and (P>=632.4*10**6) ):#glace 6
            return 6
        
        elif ((P*10**(-9)>2.67348*10**(-4)*T**(1.55299)-0.22933) and (P>=2.216*10**9)): #glace 7
            return 7
        
        #eau liquide car ne correspond a aucun type de glace
        else:
            return 0 
        
    #si on est dans le  coeur silicate solide ou liquide ?
    
    else:
        if (T>1800):
            return (-2)
        
        else :
            return (-1)

#fonction de l'itération de la température avec un developement en ordre
def calc_T(T,Q,k,rayon,Source,pas):
    
    ordre_0 = T
    ordre_1 = -1.*Q/k
    ordre_2 = -Source/k + 2*Q/(rayon*k)
    ordre_3 = 2*((ordre_2*rayon*k)-Q*k)/(rayon**2*k)
    return ((ordre_0-pas*ordre_1)-(pas**2)*(ordre_2/2))#-((pas**3)/6)*ordre_3

#intialisation des constantes
#age du systeme etudié
SSage=4.5673e9       #years

#constante de masse
masse_vol_glace_pur = 917. #kg/m^3
masse_vol_glace = 1000. #kg/m^3
masse_vol_metal = 9045.2
masse_lune = 0.47998*10**23 #kg
masse = 0. #kg

#constante géometrique
rayon_lune = 1560.8*10**3 #m
rayon_metallique = 456.25*10**3

#constante lie a la gravitation
G = 6.6743E-11 #SI
g_lune = masse_lune*G/np.power(rayon_lune,2) #ms^-2

#constante mecanique
moment_inertie = 0.346 #I/MR^2 
pression_surface = 0. #Pa

#constante thermodynamique
temperature_surface = 110. #K
chaleur_surface = 0.002 #W/m^2
emmission_radiogénique = chaleur_radiactive(SSage) # W/kg
tidal_heating = 0 #W/m^3
fichier = open("data.txt", "w")
I = [0]
var = []
while (rayon_metallique>456*10**3):
    print("rayon metallique {}".format(rayon_metallique))
    masse_vol_silicate = 2000. #kg/m^3
    I[0] = 0
    m = [100,100]
    while(masse_vol_silicate<4000 and m[0]*m[1]>0 and abs(I[0]-moment_inertie)>=(10**(-4))):
        #constante lie aux iterations
        m[1]= m[0]
        masse_vol_silicate+=10
        print("masse silicate {}".format(masse_vol_silicate))
        print("ecart {}".format((I[0]-moment_inertie)))
        pas = 100. #m
        tour_global = 0
        rayon_noyau,densite_lune = calc_forme_init(masse_lune,rayon_lune,rayon_metallique,masse_vol_glace,masse_vol_silicate,masse_vol_metal)
        nb_iteration_max = 20
        fraction_silicate_glace = 1-masse_vol_glace_pur/masse_vol_glace
        calibration = rayon_noyau
        delta_rayon = 2*pas #m
        #constantes pour eviter des erreurs de calcul 
        lim_g_estimation = 0.1*rayon_noyau #limite du basculement de methode pour calculer g
        #debut de la boucle d'iteration globale
        tour_global=0
        while (tour_global<nb_iteration_max and delta_rayon>pas):
            
            tour_global += 1 #compte le nombre de tour effectue
            tour = 0
            
            #initialisation des constantes physique pour une nouvelle iteration
            masse = 0
            rayon = rayon_lune
            etat = [1] #la surface est observé comme etant de la glace
            g = [g_lune]
            pression = [pression_surface]
            T = [temperature_surface]
            chaleur = [chaleur_surface]
            I = [0]
            masse_vol = [densite(rayon,rayon_noyau,rayon_metallique,etat[-1],pression[-1],T[-1],fraction_silicate_glace,masse_vol_silicate,masse_vol_metal)]
            emmission_condrite_totale = 0
            while (rayon>pas):#la borne inferieur peut etre modifier pour eviter les erreurs de calcul numérique et approximation
                
                tour += 1
                
                rayon -= pas
                masse += pas*4*pi*np.power(rayon,2)*masse_vol[-1]
                
                #iteration de g
                #methode direct
                if ((masse_lune-masse)>0 and rayon>=lim_g_estimation):
                    g.append((masse_lune-masse)*G/(rayon**2))
                
                #zone dans laquelle la premiere methode semble invalide du aux erreurs dans le calcul de la masse
                elif(rayon<lim_g_estimation):
                    g.append(calc_g(rayon,G,pas,masse_vol[-1]))#calcule la masse a partir du centre de la lune 
                
                else:
                    g.append(0)#si probleme renvoie la valeur 0
                    
                pression.append(pression[-1]+pas*g[-1]*masse_vol[-1])
                
                emmission_condrite_totale = 4*pi*rayon*Source(rayon,rayon_noyau,emmission_radiogénique,tidal_heating,masse_vol[-1])-tidal_heating
                
                chaleur.append(chaleur[-1]-pas*(tidal_heating+Source(rayon,rayon_noyau,emmission_radiogénique,tidal_heating,masse_vol[-1]) - 2*chaleur[-1]/rayon))
                
                #correction ou la chaleur devient negative impossible par l'équilibre thermodynamique
                if (chaleur[-1]<0):
                    chaleur[-1]=0
                    
                #evolution de la temperature
                #cherche a savoir si on est dans un liquide => convection
                #on est dans une zone solide => diffusion
                if(solide_(etat[-1])):
                    
                    T.append(calc_T(T[-1],chaleur[-1],lambda_(rayon,rayon_noyau,rayon_metallique,T[-1],etat[-1]),rayon,Source(rayon,rayon_noyau,emmission_radiogénique,tidal_heating,masse_vol[-1]),pas))
                
                else:
                    T.append(T[-1])
                    
                #incrementation du moment d'inertie de la couche
                I.append(I[-1]-(pas*8/3*np.pi*rayon**4*masse_vol[-1]/(masse_lune*rayon_lune**2)))
                
                #etablissement de l'etat de la couche
                etat.append(etat_l(rayon,rayon_noyau,rayon_metallique,pression[-1],T[-1]))
                masse_vol.append(densite(rayon,rayon_noyau,rayon_metallique,etat[-1],pression[-1],T[-1],fraction_silicate_glace,masse_vol_silicate,masse_vol_metal))
                #manque de valeurs sur la viscosite ce qui entraine de grande erreurs interprétation surement impossible en profondeur(endroit ou on a plus de galce ih)
            #fin de la boucle d'iteration pour une iteration
            
            
            #redressement moment d'inertie car incrementé de maniere negative
            for i in range(len(I)):
                
                I[i] -= I[-1]
                
            #calcule du nouveau rayon de l'interface glace silicate pour avoir la bonne masse
            masse_excedentaire = masse-masse_lune
            calibration = ((rayon_noyau**3)+(3*masse_excedentaire)/(4*pi*(masse_vol_glace-masse_vol_silicate)))**(1/3)
            delta_rayon = rayon_noyau-calibration
            rayon_noyau = calibration
            chaleur_surface=((1/3)*masse_vol_silicate*(rayon_noyau**3-rayon_metallique**3))*emmission_radiogénique/rayon_lune**2
        m[0] = I[0]-moment_inertie
    fichier.write("{}\t {}\t {}\t {}\n".format(masse_vol_silicate,rayon_metallique,rayon_noyau,abs(I[0]-moment_inertie)))
    rayon_metallique-=250
#fin des iteration globales 


#sauvegarde dans un fichier de la derniere lune iteree
fichier.close()

affichage(tour,g,pression,chaleur,T,I,etat,masse_vol,rayon_lune)