import numpy as np

def Sx(s = 1/2.):

    n = int(2 * s + 1) 
    S = np.zeros((n,n), dtype = complex)

    for row_ix, a in enumerate(range(1,n+1)):
        for column_ix, b in enumerate(range(1, n+1)):
            S[row_ix, column_ix] = (int(a == (b+1)) + int((a+1) == b)) * np.lib.scimath.sqrt((s + 1) * (a + b - 1) - (a*b))

    return S

def Sy(s = 1/2.):

    n = int(2 * s + 1) 
    S = np.zeros((n,n), dtype = complex)

    for row_ix, a in enumerate(range(1,n+1)):
        for column_ix, b in enumerate(range(1, n+1)):
            S[row_ix, column_ix] = (int(a == (b+1)) - int((a+1) == b)) * np.lib.scimath.sqrt((s + 1) * (a + b - 1) - (a*b))

    return S

def Sz(s = 1/2.):

    n = int(2 * s + 1) 
    S = np.zeros((n,n), dtype = complex)

    for row_ix, a in enumerate(range(1,n+1)):
        for column_ix, b in enumerate(range(1, n+1)):
            S[row_ix, column_ix] = (s+1-a) * int(a == b)

    return S


if __name__ == '__main__':
    s = 7/2.
    S = Sx(s)
    print('Sx:')
    print(S)
    S = Sy(s)
    print('Sy:')
    print(S)
    S = Sz(s)
    print('Sz:')
    print(S)

