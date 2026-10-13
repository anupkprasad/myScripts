import pymol
from pymol import cmd

# Load your structure first
# cmd.load("your_structure.pdb")

# Select your region of interest (optional)
cmd.select("domain", "all")

# Show cartoon in putty mode
cmd.show("cartoon", "domain")
cmd.cartoon("putty", "domain")

# Color by B-factor (gradient from blue to red)
cmd.spectrum("b", "blue_white_red", "domain")

# Adjust putty parameters (thickness variation)
cmd.set("cartoon_putty_scale_min", 0.5)   # minimum thickness
cmd.set("cartoon_putty_scale_max", 3.0)   # maximum thickness
cmd.set("cartoon_putty_scale_power", 2)   # scaling curve
cmd.set("cartoon_putty_transform", 7)     # smooth interpolation

# Rendering options
cmd.set("cartoon_transparency", 0.1, "domain")
cmd.set("ray_opaque_background", 0)  # transparent background
cmd.set("antialias", 2)
cmd.set("depth_cue", 0)
cmd.set("specular", 1)
cmd.set("shininess", 50)
cmd.set("ambient", 0.5)

# Render high resolution
cmd.ray(1000, 1000)
cmd.png("avgStrPutty.png", dpi=600)

print("Putty image saved as 'avgStrPutty.png' with B-factor coloring")

