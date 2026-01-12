# ----- Load structure -----
load kinase_with_tail.pdb, kinase

# ----- Define selections -----
# Adjust residue numbers for your protein
select kinase_domain, kinase and not resi 250-330     # example tail residues
select tail, kinase and resi 250-330                  # example tail region
select highlight_residues, kinase and resi 3+264+276+277+283

# ----- Create objects (better than selections for rendering) -----
create obj_domain, kinase_domain
create obj_tail, tail
create obj_highlight, highlight_residues

# ----- Visual styles -----

# Cartoon outline for paper-quality look
set cartoon_outline, on
set cartoon_outline_color, black
set cartoon_outline_width, 2

# Color & transparency
color white, obj_domain
set transparency, 0.5, obj_domain       # 50% transparent domain

color red, obj_tail
set transparency, 0.0, obj_tail         # opaque tail

# Highlight residues as licorice
show sticks, obj_highlight
color yellow, obj_highlight
set stick_radius, 0.25

# Improve cartoon quality
set cartoon_fancy_helices, on
set cartoon_smooth_loops, on
set cartoon_highlight_color, grey

# ----- Background transparent -----
bg_color white
set ray_opaque_background, off

# ----- Lighting settings -----
set ambient, 0.5
set direct, 0.8
set specular, 0.3
set shininess, 20

# ----- Zoom to protein -----
zoom kinase

# ----- High-resolution ray tracing -----
ray 3000, 3000
png kinase_render.png, dpi=600, ray=1

