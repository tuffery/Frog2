
def lineIntersection(A1, A2, B1, B2):
    x1, y1 = A1
    x2, y2 = A2

    x3, y3 = B1
    x4, y4 = B2

    denom = float((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1))
    if denom == 0:
        return 1
        
    ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3))/denom
    ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3))/denom

    x = x1 + ua * (x2-x1)
    y = y1 + ub * (y2-y1)

    if ua >0 and ua <1 and ub>0 and ub<1:
        return x,y, ua, ub
    else:
        return None

if __name__ == "__main__":
    A1 = 0,-1
    A2 = 0,1
    B1 = -1,0
    B2 = 1,0
    print lineIntersection(A1,A2,B1,B2)

    B1 = A1
    B2 = A2[0]+10, A2[1]
    print A1, A2
    print B1, B2

    print lineIntersection(A1,A2,B1,B2)
