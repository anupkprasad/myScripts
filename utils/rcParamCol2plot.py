"""
Created on Sun May 30 04:38:48 2021

@author: anup
"""
import matplotlib as mpl 
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/media/anup/BackupPlus/scripts/python")
############ restore default  
plt.rcdefaults()


plt.rcParams["figure.dpi"] = 300
### for one coloumn 
rc_fonts = { "font.size": 12, 'axes.titlesize': 12, "axes.labelsize": 10, "xtick.labelsize": 8,
            "ytick.labelsize": 8, "lines.linewidth": 1}

###font.size work for legend

mpl.rcParams.update(rc_fonts)
mpl.rcParams["font.weight"] = "bold" 
mpl.rcParams["axes.labelweight"] = "bold"

# set tick width
mpl.rcParams['xtick.major.size'] = 1.5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 1
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 1.5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 1
mpl.rcParams['ytick.minor.width'] = 1


# # set tick width
# mpl.rcParams['xtick.major.size'] = 1.5*2
# mpl.rcParams['xtick.major.width'] = 1
# mpl.rcParams['xtick.minor.size'] = 1*2
# mpl.rcParams['xtick.minor.width'] = 1
# mpl.rcParams['ytick.major.size'] = 1.5*2
# mpl.rcParams['ytick.major.width'] = 1
# mpl.rcParams['ytick.minor.size'] = 1*2
# mpl.rcParams['ytick.minor.width'] = 1






"""
print font family and font

# set the font globally
plt.rcParams.update({'font.family':'sans-serif'}) ### update font family
plt.rcParams.update({'font.sans-serif':'Arial'})  ### Update the font style
print(plt.rcParams['font.family'])   ### shows cuurent font family
print(plt.rcParams['font.sans-serif'])  # shows list of fonts in sans-serif
"""

# font_family = plt.rcParams['font.family']
# print("font_family is " + font_family[0])
# cur_font = 'font.' + font_family[0]
# plt.rcParams.update({'font.sans-serif':'Arial'})
# print("font style is " + plt.rcParams[cur_font][0])


# from matplotlib.font_manager import findfont, FontProperties
# font = findfont(FontProperties(family=['Arial']))
# mpl.font_manager.findfont(mpl.font_manager.FontProperties(family=font))

# print("Current font:" + font)


from matplotlib import font_manager

font_dirs = ['/media/anup/BackupPlus/scripts/python/Arial_Font/']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)

for font_file in font_files:
    font_manager.fontManager.addfont(font_file)

### set font
plt.rcParams['font.family'] = 'Arial'

plt.rcParams['mathtext.fontset'] = 'cm'

#plt.rcParams["mathtext.default"]

# import matplotlib.mathtext as mathtext
# mathtext.custom.sup2 = 0.5  # default 0.5



plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.rm'] = 'Arial'

# import matplotlib.mathtext as mathtext
# mathtext.sup1 = 0.3  # default 0.5


# import matplotlib
# matplotlib.mathtext.SHRINK_FACTOR = 0.9
# matplotlib.mathtext.GROW_FACTOR = 1 / 0.9
################# mathtext.default:regular leads to mathfont same as outside $$ fonts

# import matplotlib.pyplot as plt
# from matplotlib import mathtext

# plt.rcParams.update({
#     'text.usetex':False,
#     'font.family':'sans-serif',
#     'font.sans-serif':['Arial'],
#     'mathtext.default':'regular',
#     })


fig, ax = plt.subplots(figsize = (3.2,3.2))

x  =  np.arange(0,10,1)
y  =  2*x

ax.plot(x, y)
ax.set_xlabel(xlabel='$G-axis$')
ax.set_ylabel("Y-label")
ax.legend(["G_not work", "orange"])
ax.set_title("blueG")

# ax.set_xlabel(xlabel='$G-axis$', fontsize = 80)
# ax.legend(["G_not work", "orange"], fontsize = 30)
# ax.set_title("blueG", fontsize = 80)

plt.tight_layout()
plt.savefig('rc_param.png', dpi=300)