from smiles_parsers import Smiles, BuildMol
from perception import RingDetection, BasicAromaticity, sssr

def smilin(smiles, transforms=[sssr.sssr,
                               BasicAromaticity.aromatize],
           enable_vfgraph=0):

    """(smiles)->molecule
    Convert a smiles string into a molecule representation"""
    builder = BuildMol.BuildMol(enable_vfgraph)
    Smiles.tokenize(smiles, builder)
    mol = builder.mol

        
    for transform in transforms:
        mol = transform(mol)

## implicit hcount doesn't make any sense anymore...
    for atom in mol.atoms:
        if not atom.has_explicit_hcount:
            atom.imp_hcount = atom.hcount - atom.explicit_hcount

    return mol
