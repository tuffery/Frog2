from frowns import MDL
from frowns.mdl_parsers import sdfile


file = open("data\single2.mol")
##text = file.read()
##text.replace("\r", "\n")
##file.close()
##file = open("data\single2.mol", 'w')
##file.write(text)
##file.close()

reader = sdfile.MolReader(file)
print reader.read_one()

