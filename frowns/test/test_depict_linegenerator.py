from Numeric import *
from Tkinter import *
from frowns.Depict.LineGenerator import *
import time, random
r45 = R(0.785398163397)
print matrixmultiply(array([0,1,1]), r45)
tx = T(-1,-1)
print matrixmultiply(array([0,1,1]), tx)


test = [
    ((50,50),(100,100)),
    ((20,15), (5,5)),
    ((7,105), (10,10)),
    ((10,50), (100,50))
    ]

##   import random
for i in range(10000):
    test.append((
        (random.random()*100, random.random()*100),
        (random.random()*100, random.random()*100)))
start = time.time()
for t1, t2 in test:
    a,b,c,d = make_double_line(t1,t2)
end = time.time()
print (end-start)/len(test)

class Test:
    def __init__(self):
        self.tk = Tk()
        self.canvas = Canvas(self.tk, height=200, width=200)
        self.canvas.pack()
        import time
        tot = 0.0
        for t1, t2 in test:
            x1, y1 = t1
            x2, y2 = t2
            tt1 = time.time()
            a,b,c,d = make_double_line(t1, t2)
            tt2 = time.time()
            tot += tt2-tt1
            self.canvas.create_line(x1, y1, x2, y2)
            self.canvas.create_line(x1, y1, a[0], a[1])
            self.canvas.create_line(x1, y1, b[0], b[1])
            self.canvas.create_line(x2, y2, c[0], c[1])
            self.canvas.create_line(x2, y2, d[0], d[1])
        print tot/len(test)


Test()
mainloop()
