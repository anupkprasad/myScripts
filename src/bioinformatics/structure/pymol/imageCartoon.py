import pymol
from pymol import cmd

# cmd.select("high_bfactor", "resi 50-100 and b > 30")
for idx in ["AvgProtStrBeta"]:cmd.load(str(idx) + ".pdb","allframe")

cmd.select("domain", "resi 356-613")  
cmd.color("gray70", "domain")
cmd.show("cartoon", "domain")
cmd.set("cartoon_transparency", 0.3, "domain")




# Selection tail and apply transparency
cmd.select("tail", "resi 614-648")  # Modify selection as needed
cmd.color("red", "tail")
# Enable putty representation
cmd.show("cartoon", "tail")
cmd.set("cartoon_transparency", 0, "tail")

# Selection tail and apply transparency
cmd.select("phos", "resi 632 or resi 638")  # Modify selection as needed
cmd.show("sticks", "phos")
cmd.util.cbag("phos")



# Selection R and apply transparency
cmd.select("Rint", "resi 554 or resi 558")  # Modify selection as needed
cmd.show("sticks", "Rint")
cmd.util.cbag("Rint")



# Find and show hydrogen bonds
cmd.distance("hbonds", "phos and name O1P+O2P+OT", "Rint", cutoff=3, mode=2)
# Hide labels on H-bonds
cmd.hide("labels", "hbonds")

# Customize bond appearance
cmd.set("dash_color", "magenta", "hbonds")  # Color bonds blue
cmd.set("dash_width", 2.5)               # Increase bond width
cmd.set("dash_gap", 0.5)  # Default is 0.5

# Show the structure
cmd.show("sticks", "phos or Rint")


cmd.hide("everything", "phos and elem H")
cmd.hide("everything", "Rint and elem H")
