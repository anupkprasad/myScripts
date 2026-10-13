import pymol
from pymol import cmd

cmd.set_view(cmd.get_view())  # Lock current view
cmd.save("view_matrix.pse")
