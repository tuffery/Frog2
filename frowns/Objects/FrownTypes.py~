"""Default element properties:

The information in an element is:
  number == atomic number
  symbol == the short symbol (eg, He for helium)
  name == the full (American) English name
  mass == the atomic mass, in amu, natural propensities
  negativity == Pauling negativity value
  valences == list of possible valence occupancies

This was pretty much lifted from element.py
Thanks Andrew!
"""
# XXX FIX ME
# Let's make dictionaries instead.  That will
#  be cleaner
class AtomType:
    """AtomType
    This class essentiall reads a bunch of named arguments in the form
    a=2, b=3, c=[]
    and makes a class O with attributes O.a=2, O.b=3 O.c=[]
    Python is cool."""
    def __init__(self, symbol, **args):
        for kw, val in args.items():
            self.__dict__[kw] = val
        self.symbol = symbol
    def __str__(self):
        return "%s"%self.symbol
    
defaultAtomTypes = {}

# these represent valences in various rows of the periodic table.
valence_unknown = None
valence_1 =(1,)
valence_B =(3,)
valence_C =(4, 5)  # MODIF T. BOHME, P. TUFFERY, 22 May 2006 for the special case of doubel bond to oxygen in aromatic cycle, old: valence_C =(4,)
# Nitrogen should not be pentavalent
valence_N =(3, 5)
valence_O =(2, 3) # MODIF T. BOHME, P. TUFFERY, 22 May 2006 for the special case of oxygen in aromatic cycle, old: valence_O =(2,)
valence_P =(3, 5)
valence_S =(2, 4, 6)

num = 0
# the big list of properties
for element in (
  ("*" , 0, "unknown", 0.0),
  ("H" , 1.00794, "hydrogen", 2.20),          # 1
  ("He", 4.003, "helium", 0.0),
  ("Li", 6.941, "lithium", 0.98),
  ("Be", 9.0122, "beryllium", 1.57),
  ("B" , 10.81, "boron", 2.04),
  ("C" , 12.011, "carbon", 2.55),
  ("N" , 14.007, "nitrogen", 3.04),
  ("O" , 15.999, "oxygen", 3.44),
  ("F" , 18.998, "fluorine", 3.98),
  ("Ne", 20.179, "neon", 0.0),
  ("Na", 22.990, "sodium", 0.93),            # 11
  ("Mg", 24.305, "magnesium", 1.31),
  ("Al", 26.98, "aluminum", 1.61),
  ("Si", 28.086, "silicon", 1.90),
  ("P" , 30.974, "phosphorus", 2.19),
  ("S" , 32.066, "sulfer", 2.58),
  ("Cl", 35.453, "chlorine", 3.16),
  ("Ar", 39.948, "argon", 0.0),
  ("K" , 39.098, "potassium", 0.82),
  ("Ca", 40.08, "calcium", 1.00),
  ("Sc", 44.956, "scandium", 1.36),           # 21
  ("Ti", 47.88, "titanium", 1.54),
  ("V" , 50.94, "vanadium", 1.63),
  ("Cr", 51.996, "chromium", 1.66),
  ("Mn", 54.938, "manganese", 1.55),
  ("Fe", 55.847, "iron", 1.83),
  ("Co", 58.9332, "cobalt", 1.88),
  ("Ni", 58.69, "nickel", 1.91),
  ("Cu", 63.546, "copper", 1.90),
  ("Zn", 65.39, "zinc", 1.65),
  ("Ga", 69.72, "gallium", 1.81),             # 31
  ("Ge", 72.59, "germanium", 2.01),
  ("As", 74.922, "arsenic", 2.18),
  ("Se", 78.96, "selenium", 2.55),
  ("Br", 79.904, "bromine", 2.96),
  ("Kr", 83.80, "krypton", 0.0),
  ("Rb", 85.468, "rubidium", 0.82),
  ("Sr", 87.62, "strontium", 0.95),
  ("Y" , 88.9059, "yttrium", 1.22),
  ("Zr", 91.224, "zirconium", 1.33),
  ("Nb", 92.91, "niobium", 1.6),              # 41
  ("Mo", 95.94, "molybdenum", 2.16),
  ("Tc", 98., "technetium", 1.9),
  ("Ru", 101.07, "ruthenium", 2.2),
  ("Rh", 102.906, "rhodium", 2.28),
  ("Pd", 106.42, "palladium", 2.20),
  ("Ag", 107.868, "silver", 1.93),
  ("Cd", 112.41, "cadmium", 1.69),
  ("In", 114.82, "indium", 1.78),
  ("Sn", 118.71, "tin", 1.96),
  ("Sb", 121.76, "antimony", 2.05),           # 51
  ("Te", 127.60, "tellurium", 2.1),
  ("I" , 126.9045, "iodine", 2.66),
  ("Xe", 131.29, "xenon", 0.0),
  ("Cs", 132.91, "cesium", 0.79),
  ("Ba", 137.33, "barium", 0.89),
  ("La", 138.906, "lanthanum", 1.10),
  ("Ce", 140.12, "cerium", 1.12),
  ("Pr", 140.908, "praseodymium", 1.13),
  ("Nd", 144.24, "neodymium", 1.14),
  ("Pm", 145., "promethium", 0.0),            # 61
  ("Sm", 150.36, "samarium", 1.17),
  ("Eu", 151.96, "europium", 0.0),
  ("Gd", 157.25, "gadolinium", 1.20),
  ("Tb", 158.925, "terbium", 0.0),
  ("Dy", 162.50, "dysprosium", 1.22),
  ("Ho", 164.93, "holmium", 1.23),
  ("Er", 167.26, "erbium", 1.24),
  ("Tm", 168.934, "thulium", 1.25),
  ("Yb", 173.04, "ytterbium", 0.0),
  ("Lu", 174.967, "lutetium", 1.27),          # 71
  ("Hf", 178.49, "hafnium", 1.3),
  ("Ta", 180.95, "tantalum", 1.5),
  ("W" , 183.84, "tungsten", 2.36),
  ("Re", 186.207, "rhenium", 1.9),
  ("Os", 190.2, "osmium", 2.2),
  ("Ir", 192.22, "iridium", 2.20),
  ("Pt", 195.08, "platinum", 2.28),
  ("Au", 196.967, "gold", 2.54),
  ("Hg", 200.59, "mercury", 2.00),
  ("Tl", 204.383, "thallium", 1.62),            # 81
  ("Pb", 207.2, "lead", 1.8),
  ("Bi", 208.98, "bismuth", 2.02),
  ("Po", 209., "polonium", 2.0),
  ("At", 210., "astatine", 2.2),
  ("Rn", 222., "radon", 0.0),
  ("Fr", 223., "francium", 0.7),
  ("Ra", 226.025, "radium", 0.9),
  ("Ac", 227.028, "actinium", 1.1),
  ("Th", 232.038, "thorium", 1.3),
  ("Pa", 231.036, "protactinium", 1.5),          # 91
  ("U" , 238.029, "uranium", 1.38),
  ("Np", 237.048, "neptunium", 1.36),
  ("Pu", 244., "plutonium", 1.28),
  ("Am", 243., "americium", 1.3),
  ("Cm", 247., "curium", 1.3),
  ("Bk", 247., "berkelium", 1.3),
  ("Cf", 251., "califorium", 1.3),
  ("Es", 252., "einsteinium", 1.3),
  ("Fm", 257., "fermium", 1.3),
  ("Md", 258., "mendelevium", 1.3),           # 101
  ("No", 259., "nobelium", 1.3),
  ("Lr", 260., "lawrencium", 0.0), 
  ("Rf", 261., "rutherfordium", 0.0),
  ("Ha", 262., "hahnium", 0.0),     # also called "dubnium"
  ("Sg", 263., "seagorbium", 0.0),  # once 'unnilhexium'
  ("Ns", 269., "bohrium", 0.0),     # or "nielsbohrium"
  ("Hs", 268., "hassium", 0.0),     #  so what names do you want?
  ("Mt", 266., "meitnerium", 0.0),
  ("Uun", 269., "ununnilium", 0.0),
  ("Uuu", 272., "unununium", 0.0),             #  111
  ("Uub", 277., "ununbium", 0.0),
  # ("Uut", 0.0, "ununtrium", 0.0),   # enter when they are
  # ("Uuq", 0.0, "ununquadium", 0.0), # discovered
  # ("Uup", 0.0, "", 0.0),
  # ("Uuh", 0.0, "", 0.0),
  # ("Uus", 0.0, "", 0.0),
  # ("Uuo", 0.0, "", 0.0),
  ("R", 0, "R Group", 0.0),
  ):
    if num in (1, 9, 17, 35, 53):
	valences = valence_1
    elif num == 5:
	valences = valence_B
    elif num == 6:
	valences = valence_C
    elif num == 7:
	valences = valence_N
    elif num == 8:
	valences = valence_O
    elif num == 15:
	valences = valence_P
    elif num == 16:
	valences = valence_S
    else:
	valences = valence_unknown

    # ok, we have the info, now create the class and
    #  add it to the table
    ele = (element[0],
           num,          # number
           element[2],   # fullname
           element[1],   # mass
           element[3],   # negativity
           valences,     # valences
           num,          # equiv_class
           )
    defaultAtomTypes[element[0]] = ele
    num = num + 1


del valence_unknown
del valence_1
del valence_B
del valence_C
del valence_N
del valence_O
del valence_P
del valence_S
del num
del ele
del valences
