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


'''
KPP_REAL FUNCTION ISO1( A0,B0,C0,D0,E0,F0,G0 )
      REAL A0,B0,C0,D0,E0,F0,G0
      REAL(kind=dp) K0, K1, K2
      K0 = DBLE(D0)*EXP(DBLE(E0)/TEMP)*EXP(1.E8/TEMP**3)
      K1 = DBLE(F0)*EXP(DBLE(G0)/TEMP)
      K2 = DBLE(C0)*K0/(K0+K1)
      ISO1 =  DBLE(A0) * EXP(DBLE(B0)/TEMP) * (1.-K2)
   END FUNCTION ISO1 
'''
def ISO1(A0, B0, C0, D0, E0, F0, G0)
    K0 = D0 * EXP(E0/TEMP) * EXP(1.E8/TEMP**3)
    K1 = F0 * EXP(G0/TEMP)
    K2 = C0 * K0/(K0+K1)
    ISO1 = A0 * EXP(B0/TEMP) * (-(K2-1))
    return ISO1


'''
KPP_REAL FUNCTION ISO2( A0,B0,C0,D0,E0,F0,G0 )
  REAL A0,B0,C0,D0,E0,F0,G0
  REAL(kind=dp) K0, K1, K2
  K0 = DBLE(D0)*EXP(DBLE(E0)/TEMP)*EXP(1.E8/TEMP**3)
  K1 = DBLE(F0)*EXP(DBLE(G0)/TEMP)
  K2 = DBLE(C0)*K0/(K0+K1)
  ISO2 =  DBLE(A0) * EXP(DBLE(B0)/TEMP) * K2
END FUNCTION ISO2   
'''
def ISO2(A0, B0, C0, D0, E0, F0, G0)
    K0 = D0 * EXP(E0/TEMP) * EXP(1.E8/TEMP**3)
    K1 = F0 * EXP(G0/TEMP)
    K2 = C0 * K0/(K0+K1)
    ISO1 = A0 * EXP(B0/TEMP) * K2
    return ISO2


'''
KPP_REAL FUNCTION EPO(A1,E1,M1) 
  REAL A1, E1, M1
  REAL(kind=dp) K1      
  K1 = 1.0_dp/(DBLE(M1) * CFACTOR + 1.0_dp)
  EPO = DBLE(A1) * EXP(DBLE(E1)/TEMP) *  K1
END FUNCTION EPO
'''
def EPO(A1, E1, M1):
    K1 = 1 / (M1 * CFACTOR + 1)
    EPO = A1 * EXP(E1/TEMP) * K1
    return EPO


'''
KPP_REAL FUNCTION KCO(A1,M1) 
  REAL A1, M1
  KCO = DBLE(A1) * (1.0_dp + (CFACTOR / DBLE(M1)))
END FUNCTION KCO
'''
def KCO(A1, M1):
    KCO = A1 * (1 + (CFACTOR/M1))
    return KCO


'''
KPP_REAL FUNCTION FALL ( A0,B0,C0,A1,B1,C1,CF)
  REAL A0,B0,C0,A1,B1,C1,CF
  REAL(kind=dp) K0, K1     
  K0 = DBLE(A0) * EXP(DBLE(B0)/TEMP) * (TEMP/300.0_dp)**DBLE(C0)
  K1 = DBLE(A1) * EXP(DBLE(B1)/TEMP) * (TEMP/300.0_dp)**DBLE(C1)
  K0 = K0*CFACTOR
  K1 = K0/K1
  FALL = (K0/(1.0_dp+K1))*   &
       DBLE(CF)**(1.0_dp/(1.0_dp+(LOG10(K1))**2))
END FUNCTION FALL
'''


def FALL(A0,B0,C0,A1,B1,C1,CF):
    K0 = A0 * EXP(B0/TEMP) * (TEMP/300)**C0
    K1 = A1 * EXP(B1/TEMP) * (TEMP/300)**(C1)
    K0 = K0*CFACTOR
    K1 = K0/K1
    FALL = (K0 / (1.00+K1) * CF**(1 / (1 + LOG10(K1)) **2))
    return FALL


'''
KPP_REAL FUNCTION TROE ( A0,B0,C0,A1,B1,C1,CF)
  REAL A0,B0,C0,A1,B1,C1,CF
  REAL(kind=dp) K0, K1, KR, NC, F     
  K0 = DBLE(A0) * EXP(DBLE(B0)/TEMP) * (TEMP/300.0_dp)**DBLE(C0)
  K1 = DBLE(A1) * EXP(DBLE(B1)/TEMP) * (TEMP/300.0_dp)**DBLE(C1)
  K0 = K0*CFACTOR
  KR = K0/K1
  NC = 0.75_dp-1.27_dp*(LOG10(DBLE(CF)))
  F  = 10.0_dp**(LOG10(DBLE(CF))/(1+(LOG10(KR)/NC)**2))
  TROE = K0*K1*F/(K0+K1)
END FUNCTION TROE
'''
def TROE(A0, B0, C0, A1, B1, C1, CF):
    K0 = A0 * EXP(B0/TEMP) * (TEMP/300) * C0
    K1 = A1 * EXP(B1/TEMP) * (TEMP/300) * C1
    K0 = K0 * CFACTOR
    KR = K0/K1
    NC = 0.75 - 1.27 * LOG10(CF)
    F = 10 ** ((LOG10(CF)) / (1+(LOG10(KR)/NC)**2))
    TROE = K0*K1*F / (K0+K1)
    return TROE
