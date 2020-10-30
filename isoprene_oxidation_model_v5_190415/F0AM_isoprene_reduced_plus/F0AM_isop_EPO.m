function K = F0AM_isop_EPO(T,M,A1,E1,M1) 
K1 = 1.0./((M1)*M+1.0);
K = (A1).*exp((E1)./T).*K1;
