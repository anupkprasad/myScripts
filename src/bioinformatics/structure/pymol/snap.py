import pymol
from pymol import cmd


# Set high-quality rendering options
cmd.set("ray_opaque_background", 0)  # Transparent background
cmd.set("antialias", 2)  # Smoother image
cmd.set("depth_cue", 0)  # Removes depth shading
cmd.set("specular", 1)  # Enhances highlights
cmd.set("shininess", 50)  # Adjust shininess
cmd.set("ambient", 0.5)  # Adjust lighting
cmd.ray(2000, 2000)  # High resolution
cmd.png("avgStrPutty.png", dpi=600)  # Save as PNG with high DPI

print("High-quality image saved as 'hbond.png'")
