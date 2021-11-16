import math as m


def calc_discharge1(b, h, k_st, m_bank, S):
    A = h * (b + (h * m_bank))
    P = b + (2 * h) * m.sqrt(m_bank ** 2 + 1)
    R = A / P
    Q = k_st * m.sqrt(S) * R ** (2 / 3) * A
    return Q

#Flexibility

def calc_discharge2(b, h, m_bank, S, **kwargs):
    for char in kwargs.items():
        if 'k_st' in char[0]:
            k_st = float(char[1])
        if 'n_m' in kwargs:
            k_st = 1 / float(char[1])
        if 'D_90' in char[0]:
            k_st = 26 / (float(char[1])) ** (1/6)
    A = h * (b + (h * m_bank))
    P = b + (2 * h) * m.sqrt(m_bank ** 2 + 1)
    R = A / P
    Q = k_st * m.sqrt(S) * R ** (2 / 3) * A
    return Q


if __name__ == '__main__':
    # input parameters
    Q = 15.5        # discharge in (m3/s)
    b = 5.1         # bottom channel width (m)
    m_bank = 2.5    # bank slope
    k_st = 20       # Strickler value
    n_m = 1 / k_st  # Manning's n
    S_0 = 0.005 # channel slope
    h = 2
    # call the solver with user-defined channel geometry and discharge
print(calc_discharge1(b,h,k_st,m_bank,S_0))
print(calc_discharge2(b,h,m_bank,S_0,k_st= 20))



def interpolate_h(Q,b,m,S,**kwargs):
    for char in kwargs.items():
        if 'k_st' in char[0]:
            k_st = float(char[1])
        if 'n_m' in kwargs:
            k_st = 1 / float(char[1])
        if 'D_90' in char[0]:
            k_st = 26 / (float(char[1])) ** (1 / 6)
    eps = 1
    h = 1
    while eps > .001 :
        A = h * (b + (h * m))
        P = b + (2 * h) * (m ** 2 + 1) ** .5
        Q_k = (A ** (5/3) * S ** (.5) * k_st)/(P ** (2/3))
        eps = abs(Q - Q_k)/ Q
        F = (Q * P ** (2/3) * (1/k_st)) - (A ** (5/3) * S ** .5)
        F_d = ((Q * (2 / 3) * (1/k_st) * P ** (-1/3) * (2 * (m ** 2 + 1) ** .5)) - ((5/3) * A ** (2/3) * S ** (.5) * (b + 2 * m * h)))
        h = abs(h - F/ F_d)

    return h



print(interpolate_h(Q,b,m_bank,S_0,k_st = 20))
