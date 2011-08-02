file = open("Graph.py")
lines = file.readlines()

for line in lines:
    pos = line.find(">>>")
    if pos != -1:
        print line.rstrip()[pos+4:]
