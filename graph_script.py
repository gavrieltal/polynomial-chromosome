from sympy import symbols
from sympy.plotting import plot
import sys
import re

# argv: [./graph_script.py *savefile* *hash* *frame* *bound_x* *bound_y*]
x = symbols('x')

if len(sys.argv) != 6:
    print("Please only call me with exactly five parameters")
    sys.exit()

hash_name = sys.argv[2]
frame_num = sys.argv[3]
bound_x   = sys.argv[4]
bound_y   = sys.argv[5]

with open(sys.argv[1]) as f:
    exprs = []
    for line in f:
        exprs.append(line[:-1])
    exec("plot(" + ", ".join(exprs) + ", show = False, xlim = (-" +
             bound_x + ", +" + bound_x + "), ylim = (-" + bound_y + ", +" +
             bound_y + ")).save('" + hash_name + "_" + frame_num + "')")
