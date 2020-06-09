import numpy as np
import matplotlib.pyplot as plt
############################################################################################################
#This bit of code is from week1 it makes an array for each label in Galaxy1.txt
############################################################################################################
f = open('Galaxy1 (4).txt','r')
list_make = f.readlines()

All_data = []
Radias  = []
velocity = []
change_in_radias = []
change_in_velocity = []
Mass = []

for line in list_make:
    make_each_list = line.split('\t')
    All_data.append(make_each_list)

for i in All_data:
    Radias.append(i[0])
    velocity.append(i[1])
    change_in_radias.append(i[2])
    change_in_velocity.append(i[3])
    Mass.append(i[4])
    
Radias_Q = (Radias[1:])
velocity_Q = (velocity[1:])
change_in_radias_Q = (change_in_radias[1:])
change_in_velocity_Q = (change_in_velocity[1:])
Mass_data_n_Q = (Mass[1:])

Mass_data_Q = []

for i in Mass_data_n_Q:
    Just_data = i.replace('\n','')
    Mass_data_Q.append(Just_data)

Radias_data = []
velocity_data = []
change_in_radias_data = []
change_in_velocity_data = []
Mass_data = []

for i in Radias_Q:
    Radias_data.append(float(i))

for i in velocity_Q:
    velocity_data.append(float(i))

for i in change_in_radias_Q:
    change_in_radias_data.append(float(i))

for i in change_in_velocity_Q:
    change_in_velocity_data.append(float(i))

for i in Mass_data_Q:
    Mass_data.append(float(i))

Radias_data_array = np.array(Radias_data)
velocity_data_array = np.array(velocity_data)
change_in_radias_data_array = np.array(change_in_radias_data)
change_in_velocity_data_array = np.array(change_in_velocity_data)
Mass_data_array = np.array(Mass_data)
###############################################################################################################
###############################################################################################################


###############################################################################################################
#This is from week 2 it is plots for - v(r) for v observed, v visual, v sum. 
###############################################################################################################
plt.plot(Radias_data_array,velocity_data_array)
plt.xlabel ('Radias / kpc')
plt.ylabel ('velocity / km/s')
plt.title ('Observed velocity and calculated velocites against radius')

velocity_Visual_array = np.sqrt(((4.30*10**-6)*Mass_data_array)/Radias_data_array)

plt.plot(Radias_data_array,velocity_Visual_array)

Mass_DM_old = (4*np.pi*(100*10**6)*(1.87**2))*(Radias_data_array - (1.87*np.arctan(Radias_data_array / 1.87)))

Mass_sum = (Mass_data_array + Mass_DM_old)

velocity_sum_old = np.sqrt(((4.30*10**-6)*Mass_sum)/Radias_data_array) 

plt.plot(Radias_data_array,velocity_sum_old)
##################################################################################################################
##################################################################################################################

##################################################################################################################
#This is week 3. 
#1. Find X^2 for old p0 we were using - 100*10*6
#2. Iterate through range (guesses for optimum p0) and for each find X^2
#Then compare X^2 to old p0 X^2 and if X^2 is smaller means better that old p0 so that p0 and its X^2 added to better lists.
#Find smallest X^2 aka optimum p0.
#3. Using optimum p0 find optimum DM mass - find optimum sum mass - find optimum calculated v - use optimum v to plot optimum graph.
##################################################################################################################
#1.
Xsquared_values_old =  ((velocity_data_array - velocity_sum_old)**2)/(change_in_velocity_data_array **2)

Xsquared_old = np.sum(Xsquared_values_old)

better_guesses = []
better_Xsquared = []
##################################################################################################################
#2.

for density in np.arange((1*10**6), (200*10**6), 1000):

  Mass_density_DM_test = (4*np.pi*(density)*(1.87**2))*(Radias_data_array - (1.87*np.arctan(Radias_data_array / 1.87)))

  Cobmined_mass_test = (Mass_density_DM_test + Mass_data_array)

  v_model_test = np.sqrt(((4.30*10**-6)*Cobmined_mass_test)/Radias_data_array)

  Xsquared_values =  ((velocity_data_array - v_model_test)**2)/(change_in_velocity_data_array **2)

  X_squared = np.sum(Xsquared_values)

  if(X_squared < Xsquared_old):

    better_guesses.append(density)

    better_Xsquared.append(X_squared)
      
closest_0_X_squared = min(better_Xsquared, key=abs)

index_of_lowest_Xtwo = better_Xsquared.index(closest_0_X_squared)

Optimum_density =  (better_guesses[index_of_lowest_Xtwo])
#######################################################################################################################
#3.

Optimum_mass_DM = (4*np.pi*(Optimum_density)*(1.87**2))*(Radias_data_array - (1.87*np.arctan(Radias_data_array / 1.87)))

Optimum_Combined_mass = (Optimum_mass_DM + Mass_data_array)

v_optimum = np.sqrt(((4.30*10**-6)*Optimum_Combined_mass)/Radias_data_array)

plt.plot(Radias_data_array,v_optimum)

plt.show()
########################################################################################################################