import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt

plt.close('all')

Water = CP.AbstractState("HEOS", "Water")
pc = Water.keyed_output(CP.iP_critical)
Tc = Water.keyed_output(CP.iT_critical)
Tmin = 200
Tmax = 1000
pmax = Water.keyed_output(CP.iP_max)
pt = 611.657
Tt = 273.16
fillcolor = 'g'

fig = plt.figure()
ax = fig.add_subplot(111)

# --------------
# Melting curve
# --------------

TT = []
PP = list(np.logspace(np.log10(pt), np.log10(pmax),1000))
for p in PP:
    TT.append(Water.melting_line(CP.iT, CP.iP, p))

# #Zone VI
# for T in np.linspace(max(TT), 355):
#     TT.append(T)
#     theta = T/273.31
#     pi = 1-1.07476*(1-theta**4.6)
#     p = pi*632.4e6
#     PP.append(p)

plt.plot(np.array(TT),PP,'darkblue', label='Courbe de fusion')
# ----------------
# Saturation curve
# ----------------
Ts = np.linspace(273.16, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,'Water')

# ------
# Labels
# ------

plt.plot(Ts,ps,'orange', label='Courbe de saturation')

# Critical lines
plt.axvline(Tc, dashes = [2, 2])
plt.axhline(pc, dashes = [2, 2])

# Labels
plt.text(850, 1e8, 'supercritical',ha= 'center')
plt.text(350, 3e6, 'liquid', rotation = 45)
plt.text(450, 5e4, 'gas', rotation = 45)


plt.gca().set_yscale('log')
plt.gca().set_xlim(0, 1000)
plt.ylabel(u'Pression [Pa]')
plt.xlabel(u'Température [K]')
plt.tight_layout()
plt.legend()
plt.grid()
plt.show()

state = []

for i in range (len(T)):
    if T[i]<Water.melting_line(CP.iT, CP.iP, P[i]):
        state.append(1)
    else : 
        state.append(0)