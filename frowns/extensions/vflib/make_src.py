import os, shutil
for f in os.listdir("vflib-2.0/src"):
    head, ext = os.path.splitext(f)
    if ext == ".cc":
        src = os.path.join('vflib-2.0/src', f)
        dst = os.path.join('vflib-2.0/src', "%s.cxx"%head)
        shutil.copy(src, dst)
