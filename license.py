#!/usr/bin/python

"""
   This software is part of Frog, a chemo informatics class able to build 
   3D coordinates for small compounds
    Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier

    Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon

    Using this software, please cite:
        Frog2: Efficient 3D conformation ensemble generator for small compounds.
        Miteva MA, Guyon F, Tuffery P.
        Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

"""
 Small tool to manage licenses
 for i in *.py; do echo $i; ../license.py $i > $i.lic && mv $i.lic $i; done 
 && cat toto > $i
"""

import sys, os

license = [
'"""\n',
'   This software is part of Frog, a chemo informatics class able to build \n'
'   3D coordinates for small compounds\n',
'    Copyright (C) 2006-2007 P. Tuffery, B.O. Villoutreix, Th. Bohme Leite, D. Gomes, M. Miteva, J. Chomilier\n',
'\n',
'    Frog2 (C) 2009-2010 by P. Tuffery, M. Miteva, F. Guyon\n',
'\n',
'    Using this software, please cite:\n',
'        Frog2: Efficient 3D conformation ensemble generator for small compounds.\n',
'        Miteva MA, Guyon F, Tuffery P.\n',
'        Nucleic Acids Res. 2010 Jul;38(Web Server issue):W622-7. Epub 2010 May 5.\n',
'\n',
'    This program is free software: you can redistribute it and/or modify\n',
'    it under the terms of the GNU General Public License as published by\n',
'    the Free Software Foundation, either version 3 of the License, or\n',
'    (at your option) any later version.\n',
'\n',
'    This program is distributed in the hope that it will be useful,\n',
'    but WITHOUT ANY WARRANTY; without even the implied warranty of\n',
'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n',
'    GNU General Public License for more details.\n',
'\n',
'    You should have received a copy of the GNU General Public License\n',
'    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n',
'"""\n'
]

def detag(lines):
    isOpen = 0
    for i,l in enumerate(lines):
        if l.count("This softawre is part of Frog") or l.count("This software is part of Frog"):
            start = i-1
            isOpen = 1
        if isOpen and l.count('"""'):
            stop = i
            break
    if not isOpen: # no previous license
        if lines[0].count("python"):
            start = 2
            stop = 2
        else:
            start = 0
            stop = -1
    return lines[:start], lines[stop+1:]
            
        
def main(argv):
    f = open(argv[1])
    lines = f.readlines()
    f.close()

    head, tail = detag(lines)
    print "".join(head + license + tail)

if __name__ == "__main__":
    main(sys.argv)
