
import numpy as np
import time

EPS = 0.00001


def CompareEq(d1,d2, EPS = EPS):
    return (abs(d1-d2) < EPS)

def CompareGt(d1,d2):
    return (d1 >= d2+EPS)

def CompareGeq(d1,d2):
    return (d1 > d2 - EPS)




def Determinant(vec1, vec2, vec3):

    assert len(vec1) == 3, "wektor musi być 3-elementowy"
    assert len(vec2) == 3, "wektor musi być 3-elementowy"
    assert len(vec3) == 3, "wektor musi być 3-elementowy"

    return vec1[0]*vec2[1]*vec3[2]+vec2[0]*vec3[1]*vec1[2]+vec3[0]*vec1[1]*vec2[2]-(
                vec1[2]*vec2[1]*vec3[0]+vec2[2]*vec3[1]*vec1[0]+vec3[2]*vec1[1]*vec2[0])




def ChainRead(plik):

    f = open(plik, 'r')

    g = f.readlines()

    chain = []

    for i in range(len(g)):

        r = g[i].replace('\n', '')
        s = r.split()

        pozycje = {'x':float(s[1]), 'y':float(s[2]), 'z': float(s[3])}
        atom = {'A':pozycje, 'id': int(s[0])}
        chain.append(atom)

    return chain

def intersection(trian, l):
    '''
    mat = np.array([[l[0]['x'] - l[1]['x'], trian[1]['x'] - trian[0]['x'], trian[2]['x'] - trian[0]['x']],
                    [l[0]['y'] - l[1]['y'], trian[1]['y'] - trian[0]['y'], trian[2]['y'] - trian[0]['y']],
                    [l[0]['z']-l[1]['z'], trian[1]['z']-trian[0]['z'], trian[2]['z']-trian[0]['z']]])
    '''
    mat = [[l[0]['x'] - l[1]['x'], trian[1]['x'] - trian[0]['x'], trian[2]['x'] - trian[0]['x']],
    [l[0]['y'] - l[1]['y'], trian[1]['y'] - trian[0]['y'], trian[2]['y'] - trian[0]['y']],
    [l[0]['z'] - l[1]['z'], trian[1]['z'] - trian[0]['z'], trian[2]['z'] - trian[0]['z']]]

    detmat = Determinant(mat[0],mat[1],mat[2])
    #detmat = np.linalg.det(mat)
    #if CompareEq(detmat, 0):
    #    print((("normalny: ", detmat1), ("numpy: ", detmat)))
    #print(mat)

    #if CompareEq(detmat, 0, 0.000000000001):
    #    return 1
    if detmat == 0:
        return 1

    else:
        #matinv = np.linalg.inv(mat)
        matinv = [[None,None,None],[None,None,None],[None,None,None]]
        matinv[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / detmat
        matinv[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) / detmat
        matinv[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / detmat


        matinv[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) / detmat
        matinv[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) / detmat
        matinv[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) / detmat

        matinv[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / detmat
        matinv[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) / detmat
        matinv[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / detmat
        #print(matinv)

        t = matinv[0][0] * (l[0]['x'] - trian[0]['x']) + matinv[0][1] * (l[0]['y'] - trian[0]['y']) +matinv[0][2] * (
                    l[0]['z'] - trian[0]['z'])

        u = matinv[1][0] * (l[0]['x'] - trian[0]['x']) + matinv[1][1] * (l[0]['y'] - trian[0]['y']) + matinv[1][2] * (
                    l[0]['z'] - trian[0]['z'])

        v = matinv[2][0] * (l[0]['x'] - trian[0]['x']) + matinv[2][1] * (l[0]['y'] - trian[0]['y']) + matinv[2][2] * (
                    l[0]['z'] - trian[0]['z'])
        '''
        if t > 0 and t < 1 and u > 0 and u < 1 and v > 0 and v < 1 and (u + v < 1):
            return 1
        elif (u == 0) or (v == 0) or (u + v == 1):
            return 2
        else:
            return 0
        '''
        if CompareGt(t, 0) and CompareGt(1, t) and CompareGt(u, 0) and CompareGt(1, u) and CompareGt(v,0) and CompareGt(1, v) and CompareGt(1, u + v):
            return 1
        elif (CompareEq(u,0) and CompareGeq(v,0) and CompareGeq(1,v) and CompareGeq(t,0) and CompareGeq(1,t)) or \
                (CompareEq(v, 0) and CompareGeq(u, 0) and CompareGeq(1, u) and CompareGeq(t, 0) and CompareGeq(1, t)) or \
                (CompareEq(u + v, 1) and CompareGeq(u, 0) and CompareGeq(1, u) and CompareGeq(v, 0) and CompareGeq(1,v) and CompareGeq(t, 0) and CompareGeq(1, t)) or \
                (CompareEq(t, 0) and CompareGeq(u, 0) and CompareGeq(v, 0) and CompareGeq(1, u + v)) or \
                (CompareEq(t, 1) and CompareGeq(u, 0) and CompareGeq(v, 0) and CompareGeq(1, u + v)):
            return 2
        else:
            return 0



def ChainReduce(chain, closed = True):


    atomWasRemoved = 1
    interType = 1
    max = len(chain)
    triangle = [None,None,None]

    edge = [None,None]

    if not closed:
        max = max-2

    suma = 0
    while atomWasRemoved == 1 and len(chain) > 3:

        for i in range(max):

            for j in range(3):
                triangle[j] = chain[(i+j)%(len(chain))]['A']



            interType = 0
            s = 0
            if (i+2)%len(chain) < 2:
                s = (i+3)%len(chain)
            if i >= 2:
                for j in range(s, i-1):

                    edge[0]=chain[j]['A']
                    edge[1]=chain[j+1]['A']
                    interType = intersection(triangle,edge)

                    if interType > 0:
                        break

            if i <= len(chain)-5 and interType == 0:

                for j in range(i+3, len(chain)-1):

                    edge[0]=chain[j]['A']
                    edge[1]=chain[j+1]['A']
                    interType = intersection(triangle,edge)

                    if interType > 0:
                        break

            if closed and i < len(chain)-3 and i > 0 and interType == 0:


                edge[0] = chain[len(chain)-1]['A']
                edge[1] = chain[0]['A']
                interType = intersection(triangle, edge)

            if interType == 0:
                break

        atomWasRemoved = 0

        if interType == 0:
            #print(i+1)
            del chain[(i+1)%len(chain)]
            atomWasRemoved = 1
            max -= 1


    #print('Ile razy przeszedł while: ' + str(x) + ' W sumie ile wykonanych operacji: ' + str(suma) + ' Średnio ' + str(suma/x) + ' operacji na pętle')
    return((chain, len(chain)))