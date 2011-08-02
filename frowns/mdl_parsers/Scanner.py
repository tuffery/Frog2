"""A selection of utilities to allow for fast scanning of a mol
file"""
from __future__ import generators
import re

def molSplitter(file):
    """(file) -> for m in molSplitter(file) ... do something with mol"""
    lines = []
    for line in file:
        if line.find("$$$$") == 0:
            yield "\n".join(lines)
            lines = []
        else:
            lines.append(line)


FIELDPAT = re.compile(">\s+<([^>]+)>\s+\(*([^)]*)")
                   
def getFields(file):
    c = 0
    fields = {}
    iterator = iter(file)
    for line in iterator:
        pat = FIELDPAT.match(line)
        if pat:
            groups = pat.groups()
            field = groups[0]
            list = fields.get(field, [])
            fields[field] = list
            try:
                line = iterator.next()
                list.append(line.strip())
            except:
                # end of data
                break

    fields = fields.items()
    fields.sort()
    return fields
