"""
Created on apr27 2021

@author: anup
"""
from matplotlib import rc
from matplotlib.pyplot import *
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import savetxt
from matplotlib.animation import FuncAnimation
from celluloid import Camera

####################defnition########################


def running_mean(y, N):
    cumsum = np.cumsum(np.insert(y, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


N = 10

frametotime = 10

fig, axs = plt.subplots(1, 1, figsize=(6, 3.375))


"""

Set1 entire

"""
##################reading file########################

filepath = "/media/anup/Backup Plus/simulation/anup_data/u3.x_.15M/pep_micelle/u3.5_dimer/r7au3.5/antibeta_dimer/analysis/"

file = "ss_traj.dat"



filename = filepath + file

f = open(filename, "r")

### read and make list of line(\n) and remove space by replace()###########

data = f.read().replace(" ", "").split("\n")[:-1]

#######create array1 and put data in array ###########
array1 = [[] for i in range(len(data))]
for i in range(len(data)):
    for line in data[i][1:]:
        data1 = line
        array1[i].append(data1)


############Total resid  #######################

Tot_resid = len(array1[0])

############ array2 Total resid* frame  for sec_str timline##########


array2 = np.arange(len(data)*Tot_resid).reshape(Tot_resid, len(data))


array2 = array2.astype(str)

for i in range(len(data)):  # iteration of frames
    for j in range(Tot_resid):  # iteration of total resid
        line = data[i][j]
        (array2[:, i]) = line



array3 = np.zeros(len(array1)*8).reshape(len(array1), 8)


for j in range(len(array1)):
    array3[j, 1] = (array1[j].count("H")/Tot_resid)*100
    array3[j, 2] = (array1[j].count("G")/Tot_resid)*100
    array3[j, 3] = (array1[j].count("I")/Tot_resid)*100
    array3[j, 4] = (array1[j].count("E")/Tot_resid)*100
    array3[j, 5] = (array1[j].count("B")/Tot_resid)*100
    array3[j, 6] = (array1[j].count("T")/Tot_resid)*100
    array3[j, 7] = (array1[j].count("C")/Tot_resid)*100
    array3[j, 0] = indexs = j + 1

 ############   Averaging and Plotting  #########################

X = array3[:, 0]  # /10    ######frame to nanosecond

Y1 = array3[:, 1]
Y1 = running_mean(Y1, N)
Y2 = array3[:, 2]
Y2 = running_mean(Y2, N)
avg_file1 = Y1 + Y2


Y4 = array3[:, 4]
Y4 = running_mean(Y4, N)
Y5 = array3[:, 5]
Y5 = running_mean(Y5, N)
avg_file2 = Y4 + Y5


X = X[:len(avg_file1)]
#avg_file1  =  np.append(avg_file1,Y1[len(Y1)-N:-1])


###################Ranges#######################


axs.set_xlim([0, 500])
axs.set_ylim([0, 30])


########## Minor tick + Major ##############
axs.spines['top'].set_visible(False)
axs.spines['right'].set_visible(False)
axs.xaxis.set_minor_locator(MultipleLocator(50))
axs.xaxis.set_major_locator(MultipleLocator(100))
axs.xaxis.set_tick_params(which='both',direction= 'in',labelsize=12)
axs.yaxis.set_minor_locator(MultipleLocator(2.5))
axs.yaxis.set_major_locator(MultipleLocator(5))
axs.yaxis.set_tick_params(which='both', direction= 'in',labelsize=12)
    
############# subplot title#####################        

# axs[0].set_title('(a)', fontsize=12, fontweight='bold')
# axs[1].set_title('(b)', fontsize=12, fontweight='bold')


################### axis labels  ######################
axs.set_xlabel('Time (ns)', fontsize=12)  # X label
axs.set_ylabel('Secondary structure (%)', fontsize=12)



# line, = axs.plot([], [], lw = 3) 

# def init(): 
#     line.set_data([], [])
#     return line,
   
# def animate(i):
#     x = X *i
   
#     # plots a sine graph
#     y = (avg_file1)
#     line.set_data(x, y)
      
#     return line,

# anim = FuncAnimation(fig, animate, init_func = init,
#                      frames = 200, interval = 50, blit = False)
     
# anim.save('continuousSineWave.mp4', 
#           writer = 'ffmpeg', fps = 30)

"""
Animation
"""

x = X
y = avg_file1
y1 = avg_file2

# Initiate camera
camera = Camera(fig)

# Create individual frames
for j in range(1, 501):
    
    # Projectile's trajectory
        
    # Show Projectile's location
    axs.plot(x[:j][-1], y[:j][-1], marker='o', markersize=8, markeredgecolor='r', markerfacecolor='r')
    axs.plot(x[:j][-1], y1[:j][-1], marker='o', markersize=8, markeredgecolor='r', markerfacecolor='r')    
    
    # Show Projectile's trajectory
    t = axs.plot(x[:j], y[:j], color='red', lw=3, label="\u03B1-helix")
    t1 = axs.plot(x[:j], y1[:j], color='lime', lw=3, label="\u03B2-sheet" )
    # Hide the right and top spines
    
    j = j + 5
    
    # Capture frame
    camera.snap()

labels = ["\u03B1-helix","\u03B2-sheet"]
colors = ['red', 'lime']
handles = []
for c, l in zip(colors, labels):
    handles.append(Line2D([0], [0], color = c, label = l))

plt.legend(handles = handles, loc = 'upper left', fontsize = 12)


# Create animation
print(1)
anim = camera.animate(interval = 40)
print(2)
anim.save('animation_u3.5r7a.mp4', 
          writer = 'ffmpeg', dpi = 200)