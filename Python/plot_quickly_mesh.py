from numpy import array, append
import matplotlib.pylab as plt

#################################################
#   Code to plot quickly from multiple file     #
#       by Timothee Jamin                       #
#   The program can be copied and edited by the #
#   user.                                       #
#################################################

#############
#   Menu    #
#############

Files           = ["../Fortran/output/output_mesh.txt"]
# Sentence to start the read
Start_read      = "    E       Coeff1      Coeff2"
# Number of lines to skip after the start read
Number_line_skip= 0
# Sentence to end the read
End_read        = "STOP"
# Number of columns to read
Number_column   = 3
# Column index to plot
Column_index    = [1,2,0]
        # The first index is the x axis.
        # The second index is the y axis.

# For now all the plots will be on the same graph
y_label         = "$Coefficient 1$"
x_label         = "$Coefficient 2$"
mesh_label      = "$Energy (Hartree)$"
title_graph     = "Value of the energy for each combination of coefficient"
labels          = []    # left voided to have no labels on the graph, else, we need as much label as we have curves.
#labels = ['RHF','UHF multiplicity 1','UHF multiplicity 3']
#colors = ["orange","blue","green"]
colors          = []    # left voided to have the default colors on the graph, else, we need as much colors as we have curves.
linestyle       = ""
marker          = "o"
save_format     = "png" # png,pdf,...
Name_save       = "Mesh_coefficient.png"



#######################################################################################################################


# Declaration

Start_read_compute  = Start_read.strip().split()
End_read_compute    = End_read.strip().split()

# We create the mess where all the datas from all the files will be taken.
Data_mess   = [[array([])]*Number_column]*len(Files)
count_data_tot = []

############ Program ############

#   Read files    #
print("------Starting to read the files, please, be patient------")
for i in range(len(Files)):
    print("I will open: "+str(Files[i]))
    Do_we_start = False
    jump_count  = 0
    count_data  = 0
    Do_we_end   = False
    f = open(Files[i])

    # Collect data #

    # Data_mess[File_index][Column_index][Data_index(to append())]
    
    for ligne in f:
        words           = ligne.replace(",",".").strip().split()
        # Verification to skip the read
        if (len(words) == 0):
            continue
        for j in range(len(words)):
            if (words[j] == "nan"):
                continue
        len_data_start= 0
        for j in range(min(len(Start_read_compute),len(words))):
            if (words[j] == Start_read_compute[j] and not Do_we_end and not Do_we_start):
                len_data_start += 1
                if (len_data_start == len(Start_read_compute)):
                    print("I am starting to read")
                    Do_we_start = True

        if Do_we_start:
            if jump_count > Number_line_skip:
                for j in range(len(End_read_compute)):
                    if (words[j] == End_read_compute[j]):
                        Do_we_end = True
            if jump_count <= Number_line_skip:
                jump_count += 1
                continue
        if not Do_we_start:
            continue

        if Do_we_end:
            print("End of the file")
            count_data_tot.append(count_data)
            break


        # Reading the data mess
        for j in range(Number_column):
            if ((float(words[1]) >= 1.0 and float(words[1]) <= 3.0) or (float(words[2]) >= 1.0 and float(words[2]) <= 3.0)):
                if j == 0:
                    count_data += 1
                Data_mess[i][j]               = append(Data_mess[i][j],float(words[j]))
            #print("words =",words,i,j)
            #print(Data_mess[i][j])
    f.close()

# Plot all datas

plt.figure(0)

#if (len(colors) == 0 or len(colors) > len(Files)):
#        colors = []
#        cmap = plt.get_cmap("hsv")
#        colors = [cmap(j) for j in linspace(0,1,len(Files))]
#        #colors = [cmap(j) for j in linspace(0,1,5)]
#print(colors)
for i in range(len(Files)):
    x = array([])
    y = array([])
    z = array([])
    if len(Files) > 1:
        all_data_count = i*count_data_tot[i]
    else:
        all_data_count = 0
    if (len(labels) == 0):
        label_used = ""
    else:
        label_used = labels[i]
        
    for j in range(count_data_tot[i]):
        x = append(x,Data_mess[i][Column_index[0]][j+all_data_count])
        y = append(y,Data_mess[i][Column_index[1]][j+all_data_count])
        z = append(z,Data_mess[i][Column_index[2]][j+all_data_count])

        plt.scatter(x=x,y=y,c=z,s=20,alpha=0.9,marker=marker,label=label_used,vmax=max(z),vmin=min(z),cmap="nipy_spectral")

#plt.xlim(0,2)
#plt.ylim(0,2)
plt.ylabel(y_label)
plt.xlabel(x_label)
plt.title(title_graph)
if (len(labels) != 0):
    plt.legend(loc=4,fontsize=8)

#h = plt.contourf(x,y,z)
#plt.axis("scaled")
plt.colorbar()

plt.savefig(Name_save,format=save_format)
print("I am done!")
#plt.show()
plt.close(0)