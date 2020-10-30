function K = F0AM_isop_TROE(A0,B0,C0,A1,B1,C1,CF,T,M)
K0 = (A0)*exp((B0)./T).*(T./300.0).^(C0);
K1 = (A1)*exp((B1)./T).*(T./300.0).^(C1);
K0 = K0.*M;
KR = K0./K1;
NC = 0.75-1.27*(log10((CF)));
F  = 10.0.^(log10((CF))./(1+(log10(KR)./NC).^2));
K = K0.*K1.*F./(K0+K1);
