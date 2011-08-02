from frowns import MDL

for mol, text, error in MDL.sdin(open("../../test/data/bad.sdf")):
    if not mol:
        print "cannot parse sd file"
        print error
        print text
    else:
        print mol.cansmiles()
        print "data"
        for key, value in mol.fields.items():
            print "\t%s: %s"%(key, value)

