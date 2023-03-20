import os

print("What is the size of the matrix?")
matrix_size = int(input("-> "))

alpha = [1.0]
alpha_n = alpha[0]

for i in range(matrix_size):
    alpha_n = 1.2*alpha_n
    alpha.append(alpha_n)
print(alpha_n)


working_directory = "./input/input.txt"
working_directory_ = "./input.txt"

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
HF_in.write("p=\n")
HF_in.write("1\n")
HF_in.write("cp1\n")
HF_in.write("1\n")
HF_in.write("0\n")
HF_in.write("thr_SCF =\n")
HF_in.write("1e-8\n")
HF_in.write("--SCF_End--\n")

HF_in.close()

print("My work here is done.")
