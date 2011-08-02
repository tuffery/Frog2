import string

# This isn't used by the tokenizers, but may be needed for
# the future.  Some of the data comes from CEX, which was
# released into the public domain.  Other's come from the
# Daylight documentation, or from experience with the
# toolkit.

class Element:
    # the name "need_quote" comes from CEX
    # "hcount" is the same as CEX's "normh()"
    _fields = ("number", "symbol", "need_quote", "mass", "name",
               "electronegativity", "is_aromatic", "hcount")
    def __init__(self, **kwargs):
        for k in kwargs.keys():
            assert k in self._fields
        self.__dict__.update(kwargs)
    def __repr__(self):
        x = []
        for k in self._fields:
            x.append("%s = %s" % (k, repr(getattr(self, k))))
        return 'Element(%s)' % string.join(x, ", ")

# symbol, need_quote, mass, name, electonegativity
_element_data = (
  ("*" , 1, 0, "unknown", 0.0),
  ("H" , 1, 1.00794, "hydrogen", 2.20),          # 1
  ("He", 1, 4.003, "helium", 0.0),
  ("Li", 1, 6.941, "lithium", 0.98),
  ("Be", 1, 9.0122, "beryllium", 1.57),
  ("B" , 0, 10.81, "boron", 2.04),
  ("C" , 0, 12.001, "carbon", 2.55),
  ("N" , 0, 14.007, "nitrogen", 3.04),
  ("O" , 0, 15.999, "oxygen", 3.44),
  ("F" , 0, 18.998, "flourine", 3.98),
  ("Ne", 1, 20.179, "neon", 0.0),
  ("Na", 1, 22.990, "sodium", 0.93),            # 11
  ("Mg", 1, 24.305, "magnesium", 1.31),
  ("Al", 1, 26.98, "aluminum", 1.61),
  ("Si", 1, 28.086, "silicon", 1.90),
  ("P" , 0, 30.974, "phosphorus", 2.19),
  ("S" , 0, 32.06, "sulpher", 2.58),
  ("Cl", 0, 35.453, "chlorine", 3.16),
  ("Ar", 1, 39.948, "argon", 0.0),
  ("K" , 1, 39.098, "potassium", 0.82),
  ("Ca", 1, 40.08, "calcium", 1.00),
  ("Sc", 1, 44.956, "scandium", 1.36),           # 21
  ("Ti", 1, 47.88, "titanium", 1.54),
  ("V" , 1, 50.94, "vanadium", 1.63),
  ("Cr", 1, 51.996, "chromium", 1.66),
  ("Mn", 1, 54.938, "manganese", 1.55),
  ("Fe", 1, 55.847, "iron", 1.83),
  ("Co", 1, 58.9332, "cobolt", 1.88),
  ("Ni", 1, 58.69, "nickel", 1.91),
  ("Cu", 1, 63.546, "copper", 1.90),
  ("Zn", 1, 65.39, "zinc", 1.65),
  ("Ga", 1, 69.72, "gallium", 1.81),             # 31
  ("Ge", 1, 72.59, "germanium", 2.01),
  ("As", 1, 74.922, "arsinide", 2.18),
  ("Se", 1, 78.96, "selenium", 2.55),
  ("Br", 0, 79.904, "bromine", 2.96),
  ("Kr", 1, 83.80, "krypton", 0.0),
  ("Rb", 1, 85.468, "rubiduim", 0.82),
  ("Sr", 1, 87.62, "strontium", 0.95),
  ("Y" , 1, 88.9059, "yttrium", 1.22),
  ("Zr", 1, 91.224, "zirconium", 1.33),
  ("Nb", 1, 92.91, "niobium", 1.6),              # 41
  ("Mo", 1, 95.94, "molybdenum", 2.16),
  ("Tc", 1, 98., "technetium", 1.9),
  ("Ru", 1, 101.07, "ruthenium", 2.2),
  ("Rh", 1, 102.906, "rhodium", 2.28),
  ("Pd", 1, 106.42, "palladium", 2.20),
  ("Ag", 1, 107.868, "silver", 1.93),
  ("Cd", 1, 112.41, "cadmium", 1.69),
  ("In", 1, 114.82, "indium", 1.78),
  ("Sn", 1, 118.71, "tin", 1.96),
  ("Sb", 1, 121.75, "antimony", 2.05),           # 51
  ("Te", 1, 127.60, "tellurium", 2.1),
  ("I" , 0, 126.905, "iodine", 2.66),
  ("Xe", 1, 131.29, "xenon", 0.0),
  ("Cs", 1, 132.91, "cesium", 0.79),
  ("Ba", 1, 137.33, "barium", 0.89),
  ("La", 1, 138.906, "lanthanium", 1.10),
  ("Ce", 1, 140.12, "cerium", 1.12),
  ("Pr", 1, 140.908, "praseodymium", 1.13),
  ("Nd", 1, 144.24, "neodymium", 1.14),
  ("Pm", 1, 145., "promethium", 0.0),           # 61
  ("Sm", 1, 150.36, "samarium", 1.17),
  ("Eu", 1, 151.96, "europium", 0.0),
  ("Gd", 1, 157.25, "gadolinium", 1.20),
  ("Tb", 1, 158.925, "terbium", 0.0),
  ("Dy", 1, 162.50, "dysprosium", 1.22),
  ("Ho", 1, 164.93, "holmium", 1.23),
  ("Er", 1, 167.26, "erbium", 1.24),
  ("Tm", 1, 168.934, "thulium", 1.25),
  ("Yb", 1, 173.04, "ytterbium", 0.0),
  ("Lu", 1, 174.967, "lutetium", 1.27),          # 71
  ("Hf", 1, 178.49, "hafnium", 1.3),
  ("Ta", 1, 180.95, "tantalum", 1.5),
  ("W" , 1, 183.85, "tungsten", 2.36),
  ("Re", 1, 186.207, "rhenium", 1.9),
  ("Os", 1, 190.2, "osmium", 2.2),
  ("Ir", 1, 192.22, "iridium", 2.20),
  ("Pt", 1, 195.08, "platinum", 2.28),
  ("Au", 1, 196.967, "gold", 2.54),
  ("Hg", 1, 200.59, "mercury", 2.00),
  ("Tl", 1, 204.383, "thallium", 1.62),            # 81
  ("Pb", 1, 207.2, "lead", 1.8),
  ("Bi", 1, 208.98, "bismuth", 2.02),
  ("Po", 1, 209., "polonium", 2.0),
  ("At", 1, 210., "astatine", 2.2),
  ("Rn", 1, 222., "radon", 0.0),
  ("Fr", 1, 223., "francium", 0.7),
  ("Ra", 1, 226.025, "radium", 0.9),
  ("Ac", 1, 227.028, "actinium", 1.1),
  ("Th", 1, 232.038, "technetium", 1.3),
  ("Pa", 1, 231.036, "palladium", 1.5),          # 91
  ("U" , 1, 238.029, "uranium", 1.38),
  ("Np", 1, 237.048, "neptunium", 1.36),
  ("Pu", 1, 244., "plutonium", 1.28),
  ("Am", 1, 243., "americium", 1.3),
  ("Cm", 1, 247., "curium", 1.3),
  ("Bk", 1, 247., "berkelium", 1.3),
  ("Cf", 1, 251., "californium", 1.3),
  ("Es", 1, 252., "einsteinium", 1.3),
  ("Fm", 1, 257., "fermium", 1.3),
  ("Md", 1, 258., "mendelevium", 1.3),           # 101
  ("No", 1, 259., "nobelium", 1.3),
  ("Lr", 1, 260., "lawrencium", 0.0),
  ("Rf", 1, 261., "rutherfordium", 0.0),
  ("Ha", 1, 262., "hahnium", 0.0),     # also called "dubnium"
  ("Sg", 1, 263., "seagorbium", 0.0),  # once 'unnilhexium'
  ("Ns", 1, 269., "bohrium", 0.0),     # or "nielsbohrium"
  ("Hs", 1, 268., "hassium", 0.0),     #  so what names do you want?
  ("Mt", 1, 266., "meitnerium", 0.0),
  ("Uun", 1, 269., "ununnilium", 0.0),
  ("Uuu", 1, 272., "unununium", 0.0),             #  111
  ("Uub", 1, 277., "ununbium", 0.0),
  # ("Uut", 1, 0.0, "ununtrium", 0.0),   # enter when they are
  # ("Uuq", 1, 0.0, "ununquadium", 0.0), # discovered
  # ("Uup", 1, 0.0, "", 0.0),
  # ("Uuh", 1, 0.0, "", 0.0),
  # ("Uus", 1, 0.0, "", 0.0),
  # ("Uuo", 1, 0.0, "", 0.0),
)
# The information for the "hcount" field comes from CEX as
# H(1), B(3), C(4), N(3,5), O(2), F(1), P(3,5), S(2,4,6), Cl(1), Br(1), I(1)

_hcount_data = {
    "H": (1,),
    "B": (3,),
    "C": (4,),
    "N": (3, 5),
    "O": (2,),
    "F": (1,),
    "P": (3, 5),
    "S": (2, 4, 6),
    "Cl": (1,),
    "Br": (1,),
    "I": (1,),
}

# "P" can be aromatic, says Daylight toolkit
_aromatic_data = ("C", "N", "O", "S", "P")

elements = []
elements_by_symbol = {}

def init():
    i = 0
    for element in _element_data:
        symbol, need_quote, mass, name, electonegativity = element
        is_aromatic = symbol in _aromatic_data
        hc = _hcount_data.get(symbol, ())
        kwargs = {}
        for k, v in map(None, Element._fields,
                        (i,) + element + (is_aromatic, hc)):
            kwargs[k] = v
        e = apply(Element, (), kwargs)
        elements.append(e)
        elements_by_symbol[symbol] = e
        i = i + 1

init()

del _element_data, _hcount_data, init
