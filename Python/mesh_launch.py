import os
import subprocess
#import numpy as np


# Functions

def create_input(matrix_size,alpha):

    """The size of the matrix is the number of coefficients. The number of 'alpha' sent must be equal to the maximum of the matrix size.
    Be careful! You need to be in the "./Fortran" directory.
    """

    working_directory = "../Python/input_mesh/input.txt"
    working_directory_ = "../Python/input.txt"

    if (not os.path.isfile(working_directory)):
        HF_in = open(working_directory_,"x")
    else:
        HF_in = open(working_directory,"w")

    HF_in.write("###################################\n")
    HF_in.write("##           Input file          ##\n")
    HF_in.write("###################################\n")
    HF_in.write("\n")
    HF_in.write("##-------------------------------##\n")
    HF_in.write("##M is the dimension of the basis##\n")
    HF_in.write("###################################\n")
    HF_in.write("##We need as many Slater         ##\n")
    HF_in.write("# coefficients as the dimension  ##\n")
    HF_in.write("# of the basis set.              ##\n")
    HF_in.write("##-------------------------------##\n")
    HF_in.write("\n")
    HF_in.write("M\n")
    HF_in.write(str(matrix_size)+"\n")
    HF_in.write("--Begin--\n")

    for i in range(matrix_size):
        HF_in.write(str(alpha[i])+"\n")


    HF_in.write("--End--\n")

    HF_in.write("\n")
    HF_in.write("--SCF_Begin--\n")
    HF_in.write("p =\n")
    HF_in.write("1\n")
    HF_in.write("cp1 =\n")
    HF_in.write("1\n")
    HF_in.write("0\n")
    HF_in.write("thr_SCF =\n")
    HF_in.write("1e-8\n")
    HF_in.write("--SCF_End--\n")

    HF_in.close()

def read_ouput():
    #print("Work in progress/read_input()")
    Energy = ""
    Stop_read = " ============Converged_energy==============".strip().split()
    Output = open("./output/output.txt",'r')
    Start_read = False
    for ligne in Output:
        words           = ligne.replace(",",".").strip().split()
        # Verification to skip the read
        if Start_read:
            Energy = float(words[0])
            break
        if (len(words) == 0):
            continue
        for i in range(len(Stop_read)):
            if (words[i] == Stop_read[i]):
                Start_read = True
    return Energy

def print_energy(Energy,Coeff1,Coeff2):
    #print("Work in progress/print_energy()")

    working_directory = "./output/output_mesh.txt"

    if (not os.path.isfile(working_directory)):
        HF_out = open(working_directory,"x")
        HF_out.write("    E       Coeff1      Coeff2\n")
    else:
        HF_out = open(working_directory,"a")

    HF_out.write(str(Energy)+"    "+str(Coeff1)+"   "+str(Coeff2)+"   "+"\n")

    HF_out.close()

# Alpha values

min_value   = 0.005 # If it's 0 we have nothing.
max_value   = 2.74  # After this value program can't access the result

alpha_central_value = [round(1.45,3),round(2.90,3)]

Matrix_size         = len(alpha_central_value)
Step_alpha          = 0.005
Max_step            = max_value
Time                = round((6*60+31.49)/117,3) #s => One simulation.
Size_alpha          = int((Max_step)/(Step_alpha))


# In case it's needed
first_value_1       = min_value
end_value_1         = max_value
#first_value_2       = alpha_central_value[1]-Max_step
#end_value_2         = alpha_central_value[1]+Max_step

# Recode of linspace due to technical problems with numpy.
Coeff1 = [round(first_value_1,3)]
Coeff2 = [round(first_value_1,3)]

for i in range(Size_alpha-1):
    Coeff1.append(round(Coeff1[i]+Step_alpha,3))

for i in range(Size_alpha-1):
    Coeff2.append(round(Coeff2[i]+Step_alpha,3))
#######


#print(Coeff1, Coeff2)
print(Size_alpha)
print("Boundary of the coeff:",round(min(Coeff1),3),round(max(Coeff1),3))
#print("For coeff 2:",round(first_value_2,3),round(end_value_2,3))
print(len(Coeff1),len(Coeff2))
print("Time expected for the simulation")
print(Time*len(Coeff2)*len(Coeff1),"s = ",Time*len(Coeff2)*len(Coeff1)*1/(60*60),"h")
print("The time expected is higher than the real time!")

# Start loop

print("Should we start now? Their is no going back.")
print("Verify that you are inside the ./Fortran folder.")
Start = input("-> ")

if Start == "0" or Start == "No" or Start == "N":
    exit()

loop = 0
for i in range(len(Coeff1)):
    for j in range(len(Coeff2)):
        create_input(Matrix_size,[Coeff1[i],Coeff2[j]])
        # Launch the fortran code
        os.system("cp ../Python/input_mesh/input.txt ./input/")
        subprocess.call(['sh','./compile.sh'])

        # Retrieve the datas
        Energy = read_ouput()
        print("The Energy entreposed:", Energy, "for Coeff1 =",Coeff1[i],"and Coeff2 =",Coeff2[j])
        # Dump data inside a new file.
        print_energy(Energy=Energy,Coeff1=Coeff1[i],Coeff2=Coeff2[j])
        
        print("Security =", j, "Loop = ", loop)
        loop += 1
        # No loop security!!
        #if i == 1:
        #    HF_out = open("./output/output_mesh.txt","a")
        #    HF_out.write("STOP\n")
        #    HF_out.close()
        #    print("Emergency exit")
        #    exit()


HF_out = open("./output/output_mesh.txt","a")
HF_out.write("STOP\n")
HF_out.close()


print("My work here is done.")


