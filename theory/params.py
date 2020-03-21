#!/usr/bin/env python
import numpy as np


MZ   = 91.1876
MW   = 80.385
GF   = 1.1663787e-5
MN   = 0.938272046

MZ2  = MZ * MZ
MW2  = MW * MW
MN2  = MN * MN
GF2  = GF * GF


def get_stot(mA,EA,mB,EB):
    """
    mA : mass of incoming beam A 
    mB : mass of incoming beam B
    EA : energy of incoming beam A 
    EB : energy of incoming beam B
    """
    pA = np.sqrt(EA**2-mA**2)
    pB = np.sqrt(EB**2-mB**2)
    return mA**2+mB**2 + 2*(EA*EB + pA*pB)


if __name__=="__main__":


    me,Ee=0,18.0
    mN,EN=MN, 275.
    print(get_stot(me,Ee,mN,EN))




