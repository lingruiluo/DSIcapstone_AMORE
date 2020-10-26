# isoprene_rates
# isoprene_rates_v5.rtf
import math
def EXP(x):
    return(math.exp(x))

def LOG10(x):
    return(math.log10(x))
'''
KPP_REAL FUNCTION TUN( A0,B0,C0 )
      REAL A0,B0,C0
      TUN =  DBLE(A0) * EXP(-DBLE(B0)/TEMP) * EXP(DBLE(C0)/TEMP**3)
   END FUNCTION TUN
'''
def TUN(A0, B0, C0):
    return(A0 * EXP(B0/TEMP) * EXP(C0/TEMP**3))

'''
KPP_REAL FUNCTION ALK ( A0,B0,C0,n,X0,Y0)
      REAL A0,B0,C0,n,X0,Y0
      REAL(kind=dp) K0, K1, K2, K3, K4
      K0 = 2.0E-22_dp * EXP(DBLE(n))
      K1 = 4.3E-1_dp*(TEMP/298.0_dp)**(-8)
      K0 = K0*CFACTOR
      K1 = K0/K1
      K2 = (K0/(1.0_dp+K1))*   &
           4.1E-1_dp**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
      K3 = DBLE(C0)/(K2+DBLE(C0))
      K4 = DBLE(A0)*(DBLE(X0)-TEMP*DBLE(Y0))
      ALK = K4 * EXP(DBLE(B0)/TEMP) * K3
   END FUNCTION ALK
'''
def ALK(A0, B0, C0, n, X0, Y0):
    K0 = 2.0E-22 * EXP(n)
    K1 = 4.3E-1 * (TEMP/298.0) ** (-8)
    K0 = K0 * CFACTOR
    K1 = K0/K1
    K2 = (K0/(1.0 + K1)) * 4.1E-1 ** (1.0/(1.0 + (LOG10(K1)) ** 2))
    K3 = C0/(K2 + C0)
    K4 = A0 * (X0 - TEMP * Y0)
    ALK = K4 * EXP(B0/TEMP) * K3
    return(ALK)

'''
KPP_REAL FUNCTION NIT ( A0,B0,C0,n,X0,Y0)
      REAL A0,B0,C0,n,X0,Y0
      REAL(kind=dp) K0, K1, K2, K3, K4
      K0 = 2.0E-22_dp * EXP(DBLE(n))
      K1 = 4.3E-1_dp*(TEMP/298.0_dp)**(-8)
      K0 = K0*CFACTOR
      K1 = K0/K1
      K2 = (K0/(1.0_dp+K1))*   &
           4.1E-1_dp**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
      K3 = K2/(K2+DBLE(C0))
      K4 = DBLE(A0)*(DBLE(X0)-TEMP*DBLE(Y0))
      NIT = K4 * EXP(DBLE(B0)/TEMP) * K3
   END FUNCTION NIT
'''
def NIT(A0, B0, C0, n, X0, Y0):
    K0 = 2.0E-22 * EXP(n)
    K1 = 4.3E-1 *(TEMP/298.0) ** (-8)
    K0 = K0 * CFACTOR
    K1 = K0/K1
    K2 = (K0/(1.0 + K1)) * 4.1E-1 ** (1.0 /(1.0 +(LOG10(K1)) ** 2))
    K3 = K2/(K2 + C0)
    K4 = A0*(X0 - TEMP * Y0)
    NIT = K4 * EXP(B0/TEMP) * K3
    return(NIT)
