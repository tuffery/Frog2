"""Routines to make lines
"""
from Numeric import *
from math import cos, sin, pi, atan2


def R(theta):
    array = zeros([3,3], Float)
    array[0,0] = array[1,1] = cos(theta)
    array[0,1] = sin(theta)
    array[1,0] = -sin(theta)
    array[2,2] = 1
    return array

def T(tx, ty):
    array = zeros([3,3], Float)
    array[0,0] = array[1,1] = array[2,2] = 1.0
    array[2,0] = tx
    array[2,1] = ty
    return array

offset = 0.2
d1 = array((0,offset,1))
d2 = array((0,-offset,1))


def make_double_line_slow(p1, p2):
    x1,y1 = p1
    x2,y2 = p2
    dx = x2-x1
    dy = y2-y1
    theta = atan2(dy,dx)
    d = math.sqrt(dx*dx + dy*dy)
    d3 = array((d,offset,1))
    d4 = array((d,-offset,1))
    M = matrixmultiply(R(theta), T(x1,y1))
    a1 = matrixmultiply(d1, M)[0:2]
    a2 = matrixmultiply(d2, M)[0:2]
    a3 = matrixmultiply(d3, M)[0:2]
    a4 = matrixmultiply(d4, M)[0:2]
    return a1, a2, a3, a4

R90 = R(pi/2.0)

def make_double_line(p1, p2, split=3):
    x1,y1 = p1
    x2,y2 = p2
    dx = (x2-x1)
    dy = (y2-y1)
    d = math.sqrt(dx*dx + dy*dy)
    dx *= 0.05
    dy *= 0.05

    ox, oy = matrixmultiply([dx,dy,1], R90)[0:2]
    return ( (x1+ox, y1+oy),
             (x1-ox, y1-oy),
             (x2+ox, y2+oy),
             (x2-ox, y2-oy))

