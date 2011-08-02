from frowns import MDL
from frowns.mdl_parsers import sdfile


file = open("data\\bad_compounds.sdf")

for mol, error, record in MDL.sdin(file):
    print type(mol), `error`, len(record)

