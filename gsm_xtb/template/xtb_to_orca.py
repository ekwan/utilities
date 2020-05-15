# converts an xtb output file into a fake ORCA output file for use by gsm
# Eugene Kwan, May 2020
#
# usage:
# python xtb_to_orca.py abc.xtbout
#
# generates abc.engrad and abc.out

import sys, os

# define elements
element_symbols=["Bq","H","He","Li","Be","B","C","N","O","F","Ne",
                 "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
                 "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge",
                 "As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc",
                 "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                 "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb",
                 "Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os",
                 "Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr",
                 "Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf",
                 "Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
                 "Ds","Rg","Cn","Ut","Fl","Up","Lv","Us","Uo"]

element_symbols_lc = [ symbol.lower() for symbol in element_symbols ]

# get xtb filename
argcount = len(sys.argv)
assert argcount == 2, f"Expected 2 command line arguments, but got {argcount}!"
basename = str(sys.argv[1])
xtb_filename = f"{basename}.xtbout"

# check necessary files are present
assert os.path.isfile(xtb_filename), f"Could not find xtb output file ({xtb_filename})!"
assert os.path.isfile("energy"), f"Could not find xtb energy file!"
assert os.path.isfile("gradient"), f"COuld not find xtb gradient file!"

# read energy
with open("energy", "r") as f:
    lines = f.readlines()
    energy = lines[1].split()[1]
    energy = float(energy)

# read geometry and gradient
with open("gradient", "r") as f:
    gradlines = f.readlines()
    n_atoms = ( len(gradlines)-3 ) / 2
    n_atoms = int(n_atoms)
    geometry_lines = gradlines[2:2+n_atoms]
    gradient_lines = gradlines[2+n_atoms:2+2*n_atoms]

# write fake ORCA engrad file
output_filename = f"{basename}.engrad".replace("structure","orcain")
with open(output_filename, "w") as engrad:
    engrad.write("#\n# Number of atoms\n#\n")
    engrad.write(f"{n_atoms:3d}\n")
    engrad.write(f"#\n# The current total energy in Eh\n#\n{energy:20.12f}\n")
    engrad.write("#\n# The current gradient in Eh/bohr\n#\n")
    for gl in gradient_lines:
        glval = gl.strip().split()
        for gv in glval:
            gv = gv.replace("D","E")
            gv = float(gv)
            engrad.write(f"{gv:21.12f}\n")
    engrad.write(u"#\n# The atomic numbers and current coordinates in Bohr\n#\n")
    for cl in geometry_lines:
        cl = cl.replace("D","E")
        cls = cl.strip().split()
        coords = [float(cls[0]),float(cls[1]),float(cls[2])]
        symbol = cls[3]
        assert symbol.lower() in element_symbols_lc, f"Unrecognized element symbol: {symbol}"
        atn = element_symbols_lc.index(symbol.lower())
        engrad.write(f"{atn:4d} {coords[0]:13.7f} {coords[1]:13.7f} {coords[2]:13.7f}\n")

# write fake ORCA output file
output_filename = f"{basename}.out".replace("structure","orcain")
with open(output_filename, "w") as outfile:
    outfile.write("ORCA-Dummy output for GSM TS Optimizer. Not a real ORCA-Output\n")
    outfile.write(f"Total Energy       :  {energy:20.8f}\n")
    outfile.write("------------------\nCARTESIAN GRADIENT\n------------------\n\n")

    for i in range(n_atoms):
        gl = gradient_lines[i].replace("D","E")
        cl = geometry_lines[i].replace("D","E")
        cls = cl.strip().split()
        gls = gl.strip().split()
        grads = [float(gls[0]),float(gls[1]),float(gls[2])]
        outfile.write(f"{i:4d}{cls[3]:>4}   :{grads[0]:15.9f}{grads[1]:15.9f}{grads[2]:15.9f}\n")

# rename original files
os.rename("energy", f"{basename}.energy")
os.rename("gradient", f"{basename}.gradient")
