import time
import random
import numpy as np
import os
EPS = 0.00001

def CompareEq(d1,d2, EPS = EPS):         #Funkcje porównujące epsilonowo
    return (abs(d1-d2) < EPS)

def CompareGt(d1,d2):
    return (d1 >= d2+EPS)

def CompareGeq(d1,d2):
    return (d1 > d2 - EPS)

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
    #chain = np.array(chain)
    return chain


def strip(chain):
    '''
    pozbądź się zbędnych wartości, tak żeby zostały same liczby
    daje listę postaci [[x,y,z], [x,y,z] ... ,[x,y,z]] gdzie każda wewnętrzna lista to wsp. kolejnego atomu
    '''
    stripped = [list(i['A'].values()) for i in chain]
    return np.array(stripped)



def pierwszy_wiersz(punkty):      #Tworzy macierze dla pojedynczego trójkąta
    '''
    Tworzy tensor Nx3x3 dla danego trójkąta oraz macierz Nx3, gdzie elementami są listy liczb do liczenia tuv
    Trójkąt[i] ma macierz w wiersz[i] oraz liczby do tuv w do_liczenia_tuv[i] a odpowiada to wszystko punktowi[i+1]
    '''
    N = len(punkty)

    wiersz = np.empty((N,3,3),dtype=np.float64)

    do_liczenia_tuv = np.empty((N,3),dtype=np.float64)


    tr = punkty[:3,:]

    kr_1 = punkty
    kr_2 = np.roll(punkty, -3)

    do_liczenia_tuv = kr_1 - tr[0]

    wiersz[:,:,1:] = np.array([tr[1]-tr[0],tr[2]-tr[0]]).T
    wiersz[:,:,0] = kr_1-kr_2

    return wiersz, do_liczenia_tuv




def odwroc_wiersz(wiersz):
    '''
    Dla każdej macierzy w wierszu liczy wyznacznik
    Porównuje te wyznaczniki z 0 omijając wiersz[0] i wiersz[1], jeśli gdzieś jest 0 to nie odwraca
    Odwraca resztę macierzy i je zwraca
    '''

    detmat = np.linalg.det(wiersz)
    comp = CompareEq(detmat[2:], 0)

    if np.any(comp):
        return None
    else:
        inv_wiersz = np.linalg.inv(wiersz[2:])
        return inv_wiersz


def intersection(inv_wiersz, do_liczenia_tuv):
    '''
    Bierze odwrócone macierze i liczby do liczenia tuv
    Jeśli gdzieś poza pierwszymi dwoma macierzami był wyznacznik == 0 to zwraca 1
    W przeciwnym razie tworzy listy wartości: t u v
    i liczy warunek na przecięcie dla każdej krawędzi
    Jeśli jakaś krawędź ma wszystkie True to zwracamy 1 = przecięcie, jeśli nie to 0 = można usunąć
    '''
    if inv_wiersz is None:
        return 1
    else:
        t = (inv_wiersz[:,0]*do_liczenia_tuv[2:]).sum(axis=1)
        u = (inv_wiersz[:,1]*do_liczenia_tuv[2:]).sum(axis=1)
        v = (inv_wiersz[:,2]*do_liczenia_tuv[2:]).sum(axis=1)

        x = np.array([CompareGt(t, 0) ,CompareGt(1, t) , CompareGt(u, 0) , CompareGt(1, u) ,CompareGt(v,0) , CompareGt(1, v) , CompareGt(1, u + v)])
        y = np.all(x, axis=0)

        if np.any(y):
            return 1
        else:
            return 0


def redukcja(plik, test=0):

    if test == 0:


        chain = ChainRead(plik)

        stripped = strip(chain)
    else:
        stripped = plik

    max = len(stripped)

    atomWasRemoved = 1
    while len(stripped) >= 3 and atomWasRemoved:

        for i in range(max):

            wiersz, do_liczenia_tuv = pierwszy_wiersz(stripped)

            inv_wiersz = odwroc_wiersz(wiersz)

            inter_type = intersection(inv_wiersz, do_liczenia_tuv)

            if inter_type > 0:
                stripped = np.roll(stripped,-3)

            else:
                break

        atomWasRemoved = 0

        if inter_type == 0:

            stripped = np.delete(stripped,1,0)
            atomWasRemoved = 1
            max -= 1

    return stripped, len(stripped)
