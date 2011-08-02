"""figueras ring detection implementation.

Usage
molecule = sssr(molecule)

adds rings and cycles to a molecule.
"""

try:
    # try to import the C++ version
    from pysssr import sssr
except ImportError:
    # otherwise import the straight python version
    from figueras import sssr
