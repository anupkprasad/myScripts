import pymol
from pymol import cmd

# cmd.select("high_bfactor", "resi 50-100 and b > 30")


# Selection tail and apply transparency
cmd.select("tail", "resi 614-648")  # Modify selection as needed
cmd.color("red", "tail")
# Enable putty representation
cmd.show("cartoon", "tail")
cmd.cartoon("putty", "tail")
cmd.set("cartoon_putty_transform", 7)
#cmd.set("cartoon_transparency", 0.3, "tail")


# Selection domain and apply transparency
cmd.select("domain", "resi 356-613")  # Modify selection as needed
cmd.color("gray50", "domain")
# Enable putty representation
cmd.show("cartoon", "domain")
cmd.set("cartoon_transparency", 0.3, "domain")
cmd.cartoon("putty", "domain")
cmd.set("cartoon_putty_transform", 7)
# Set high-quality rendering options
cmd.set("ray_opaque_background", 0)  # Transparent background
cmd.set("antialias", 2)  # Smoother image
cmd.set("depth_cue", 0)  # Removes depth shading
cmd.set("specular", 1)  # Enhances highlights
cmd.set("shininess", 50)  # Adjust shininess
cmd.set("ambient", 0.5)  # Adjust lighting
cmd.ray(2000, 2000)  # High resolution
cmd.png("avgStrPutty.png", dpi=600)  # Save as PNG with high DPI

print("High-quality image saved as 'avgStrPutty.png'")