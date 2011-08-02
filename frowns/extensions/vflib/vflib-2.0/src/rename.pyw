import os

files = os.listdir(os.curdir)

for f in files:
    if f[-2:] == "cc":
        newf = f[:-2] + "cpp"
        print f, newf
        os.system("copy %s %s"%(f, newf))
        #os.remove(f)
