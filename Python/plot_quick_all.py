from numpy import array, append
import matplotlib.pylab as plt



f = open("./H2_RHF_multi1/H2_plot.dat","r")
f.readline()
N_RHF             = array([])
B1_RHF            = array([])
SCF_RHF           = array([])
for ligne in f:
    words           = ligne.replace(",",".").strip().split()
    if (words[0] == "----"):
        continue
    N_RHF               = append(N_RHF,float(words[0]))
    B1_RHF              = append(B1_RHF,float(words[1]))
    SCF_RHF             = append(SCF_RHF,float(words[2]))


x_RHF = B1_RHF
y_RHF = SCF_RHF
print(x_RHF)
print(y_RHF)


f = open("./H2_UHF_multi1/H2_plot.dat","r")
f.readline()
N_UHF             = array([])
B1_UHF            = array([])
SCF_UHF           = array([])
for ligne in f:
    words           = ligne.replace(",",".").strip().split()
    if (words[0] == "----"):
        continue
    N_UHF               = append(N_UHF,float(words[0]))
    B1_UHF              = append(B1_UHF,float(words[1]))
    SCF_UHF             = append(SCF_UHF,float(words[2]))


x_UHF = B1_UHF
y_UHF = SCF_UHF
print(x_UHF)
print(y_UHF)

f = open("./H2_UHF_multi3/H2_plot.dat","r")
f.readline()
N_UHF3             = array([])
B1_UHF3            = array([])
SCF_UHF3           = array([])
for ligne in f:
    words           = ligne.replace(",",".").strip().split()
    if (words[0] == "----"):
        continue
    N_UHF3               = append(N_UHF3,float(words[0]))
    B1_UHF3              = append(B1_UHF3,float(words[1]))
    SCF_UHF3             = append(SCF_UHF3,float(words[2]))


x_UHF3 = B1_UHF3
y_UHF3 = SCF_UHF3
print(x_UHF3)
print(y_UHF3)


plt.figure(0)


plt.plot(x_RHF,y_RHF,linestyle='-',color='orange', label='RHF')
plt.plot(x_UHF,y_UHF,linestyle='-',color='blue', label='UHF multiplicity 1')
plt.plot(x_UHF3,y_UHF3,linestyle='-',color='green', label='UHF multiplicity 3')
plt.ylabel("Energy (eV)")
plt.xlabel("$Distance (Bhor)$")
plt.title("$H_2$ dissociation with 3 different methods")
plt.legend(loc=4,fontsize=8)

plt.savefig("H2_dissociation.png",format="png")
plt.show()
plt.close(0)
