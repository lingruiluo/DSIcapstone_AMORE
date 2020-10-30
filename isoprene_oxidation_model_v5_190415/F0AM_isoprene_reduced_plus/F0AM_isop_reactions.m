% Isoprene_reduced_plus
% with MCM inorganics
% 20190201
% # of species = 155
% # of reactions = 429

SpeciesToAdd = {...
'HCHO'; 'CH3NO3'; 'CH3OH'; 'O1D'; 'O3'; 'HO2NO2'; 'NO3'; 'N2O5'; 'H2O2'; 'NO'; ...
'NA'; 'HO2'; 'NO2'; 'HSO3'; 'CO'; 'O'; 'HNO3'; 'SO3'; 'SO2'; 'CH3O'; ...
'OH'; 'H2'; 'HONO'; 'CH3OONO2'; 'CH3OOH'; 'SA'; 'CH3OO'; 'HCOCO'; 'HCOOH'; 'HAC';...
'MGLY'; 'HPETHNL'; 'CH3CO3'; 'GLYC'; 'GLYX'; 'HPAC'; 'PYRAC'; 'SCI'; 'HMHP'; 'PROPNN';...
'ETHLN'; 'NCH2CO3'; 'HOCH2CO3'; 'HOCH2CO2H'; 'HOCH2CO3H'; 'PHAN'; 'NCH2CO2H'; 'NCH2CO3H'; 'PNAN'; 'CH3CO2H';...
'CH3CO3H'; 'PAN'; 'HMML'; 'ISOP'; 'IHOO1'; 'IHOO4'; 'ISOP1CO4OOH'; 'ISOP1OOH4CO'; 'ISOP1CO4OH'; 'ISOP1OH4CO';...
'ISOP3CO4OH'; 'ISOP1CO2OOH'; 'ISOP3OOH4CO'; 'ICHOO'; 'ISOP1OH2N'; 'ISOP1OH4N'; 'ISOP1N4OH'; 'ISOP3N4OH'; 'ISOP1OH2OOH'; 'ISOP3OOH4OH';...
'ISOP1OH4OOH'; 'ISOP1OOH4OH'; 'IHPOO1'; 'IHPOO2'; 'IHPOO3'; 'IEPOXt'; 'IEPOXc'; 'IEPOXD'; 'ICHE'; 'IEPOXAOO';... 
'IEPOXBOO'; 'MPAN'; 'ISOP1CO4CO'; 'C4HVP1'; 'C4HVP2'; 'HPALD1OO'; 'HPALD2OO'; 'ISOPNOO1'; 'ISOPNOO2'; 'INO2B'; ...
'INO2D'; 'INPB'; 'INPD'; 'ISOP1N4CO'; 'ISOP1CO4N'; 'ISOP3CO4N'; 'IHNB'; 'INO'; 'IDN'; 'IDHNBOO';...
'IDHNDOO1'; 'IDHNDOO2'; 'IHPNBOO'; 'IHPNDOO'; 'IHNE'; 'ICNE'; 'IHNEOO'; 'ICN1OO'; 'ICN2OO'; 'ICN3OO';...
'ICN4OO'; 'ICN5OO'; 'ICHNP'; 'IHNDP'; 'IDHDN'; 'IDHPN'; 'IDHCN'; 'ICPDH'; 'IDHDP'; 'IDHPE';...
'IHPDN'; 'ICHDN'; 'IHNPE'; 'IDCHP'; 'MVK'; 'MACR'; 'MVK3CO4N'; 'MVK3OH4OH'; 'MACR2OH3OH'; 'MVK3OOH4OH';...
'MACR2OOH3OH'; 'MVK3CO4OH'; 'MVK3OH4CO'; 'MACR2OH3CO'; 'MVKOHOO'; 'MVK3OH4OOH'; 'MVK3OH4N'; 'MVK3N4OH'; 'MCROHOO'; 'MACR1OO';...
'MACR1OOH'; 'MACR1OH'; 'MACR2N3OH'; 'MVK3OOH4CO'; 'MVKENOL'; 'MCRENOL'; 'MVK3OOH4N'; 'MACR2OOH3N'; 'MACR2OH3N'; 'MACR2OH3OOH';...
'IDNOO'; 'MACRNO2'; 'MACRNOH'; 'MACRNOOH'; 'MPANHN'};

% RO2ToAdd = {'CH3OO'; 'CH3CO3'};

AddSpecies

% MCM INORGANIC & C1 SCHEME

i=i+1;
Rnames{i} = 'O = O3';
k(:,i) = 5.6e-34.*.78.*M.*(T./300).^-2.6.*.21.*M;
Gstr{i,1} = 'O'; 
fO(i)=fO(i)-1; fO3(i)=fO3(i)+1; 

i=i+1;
Rnames{i} = 'O = O3';
k(:,i) = 6.0e-34.*.21.*M.*(T./300).^-2.6.*.21.*M;
Gstr{i,1} = 'O'; 
fO(i)=fO(i)-1; fO3(i)=fO3(i)+1; 

i=i+1;
Rnames{i} = 'O + O3 =';
k(:,i) = 8.0e-12.*exp(-2060./T);
Gstr{i,1} = 'O'; Gstr{i,2} = 'O3'; 
fO(i)=fO(i)-1; fO3(i)=fO3(i)-1; 

i=i+1;
Rnames{i} = 'O + NO = NO2';
k(:,i) = KMT01;
Gstr{i,1} = 'O'; Gstr{i,2} = 'NO'; 
fO(i)=fO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'O + NO2 = NO';
k(:,i) = 5.5e-12.*exp(188./T);
Gstr{i,1} = 'O'; Gstr{i,2} = 'NO2'; 
fO(i)=fO(i)-1; fNO2(i)=fNO2(i)-1; fNO(i)=fNO(i)+1; 

i=i+1;
Rnames{i} = 'O + NO2 = NO3';
k(:,i) = KMT02;
Gstr{i,1} = 'O'; Gstr{i,2} = 'NO2'; 
fO(i)=fO(i)-1; fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)+1; 

i=i+1;
Rnames{i} = 'O1D = O';
k(:,i) = 3.2e-11.*exp(67./T).*.21.*M;
Gstr{i,1} = 'O1D'; 
fO1D(i)=fO1D(i)-1; fO(i)=fO(i)+1; 

i=i+1;
Rnames{i} = 'O1D = O';
k(:,i) = 2.0e-11.*exp(130./T).*.78.*M;
Gstr{i,1} = 'O1D'; 
fO1D(i)=fO1D(i)-1; fO(i)=fO(i)+1; 

i=i+1;
Rnames{i} = 'NO + O3 = NO2';
k(:,i) = 1.4e-12.*exp(-1310./T);
Gstr{i,1} = 'NO'; Gstr{i,2} = 'O3'; 
fNO(i)=fNO(i)-1; fO3(i)=fO3(i)-1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NO2 + O3 = NO3';
k(:,i) = 1.4e-13.*exp(-2470./T);
Gstr{i,1} = 'NO2'; Gstr{i,2} = 'O3'; 
fNO2(i)=fNO2(i)-1; fO3(i)=fO3(i)-1; fNO3(i)=fNO3(i)+1; 

i=i+1;
Rnames{i} = 'NO + NO = NO2 + NO2';
k(:,i) = 3.3e-39.*exp(530./T).*.21.*M;
Gstr{i,1} = 'NO'; Gstr{i,2} = 'NO'; 
fNO(i)=fNO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NO + NO3 = NO2 + NO2';
k(:,i) = 1.8e-11.*exp(110./T);
Gstr{i,1} = 'NO'; Gstr{i,2} = 'NO3'; 
fNO(i)=fNO(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NO2 + NO3 = NO + NO2';
k(:,i) = 4.50e-14.*exp(-1260./T);
Gstr{i,1} = 'NO2'; Gstr{i,2} = 'NO3'; 
fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)-1; fNO(i)=fNO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NO2 + NO3 = N2O5';
k(:,i) = KMT03;
Gstr{i,1} = 'NO2'; Gstr{i,2} = 'NO3'; 
fNO2(i)=fNO2(i)-1; fNO3(i)=fNO3(i)-1; fN2O5(i)=fN2O5(i)+1; 

i=i+1;
Rnames{i} = 'O1D = OH + OH';
k(:,i) = 2.14e-10.*H2O;
Gstr{i,1} = 'O1D'; 
fO1D(i)=fO1D(i)-1; fOH(i)=fOH(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'OH + O3 = HO2';
k(:,i) = 1.70e-12.*exp(-940./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'O3'; 
fOH(i)=fOH(i)-1; fO3(i)=fO3(i)-1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + H2 = HO2';
k(:,i) = 7.7e-12.*exp(-2100./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'H2'; 
fOH(i)=fOH(i)-1; fH2(i)=fH2(i)-1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + CO = HO2';
k(:,i) = KMT05;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'CO'; 
fOH(i)=fOH(i)-1; fCO(i)=fCO(i)-1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + H2O2 = HO2';
k(:,i) = 2.9e-12.*exp(-160./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'H2O2'; 
fOH(i)=fOH(i)-1; fH2O2(i)=fH2O2(i)-1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HO2 + O3 = OH';
k(:,i) = 2.03e-16.*(T./300).^4.57.*exp(693./T);
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'O3'; 
fHO2(i)=fHO2(i)-1; fO3(i)=fO3(i)-1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'OH + HO2 =';
k(:,i) = 4.8e-11.*exp(250./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'HO2'; 
fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)-1; 

i=i+1;
Rnames{i} = 'HO2 + HO2 = H2O2';
k(:,i) = 2.20e-13.*KMT06.*exp(600./T);
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'HO2'; 
fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)-1; fH2O2(i)=fH2O2(i)+1; 

i=i+1;
Rnames{i} = 'HO2 + HO2 = H2O2';
k(:,i) = 1.90e-33.*M.*KMT06.*exp(980./T);
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'HO2'; 
fHO2(i)=fHO2(i)-1; fHO2(i)=fHO2(i)-1; fH2O2(i)=fH2O2(i)+1; 

i=i+1;
Rnames{i} = 'OH + NO = HONO';
k(:,i) = KMT07;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'NO'; 
fOH(i)=fOH(i)-1; fNO(i)=fNO(i)-1; fHONO(i)=fHONO(i)+1; 

i=i+1;
Rnames{i} = 'OH + NO2 = HNO3';
k(:,i) = KMT08;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'NO2'; 
fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)-1; fHNO3(i)=fHNO3(i)+1; 

i=i+1;
Rnames{i} = 'OH + NO3 = HO2 + NO2';
k(:,i) = 2.0e-11;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'NO3'; 
fOH(i)=fOH(i)-1; fNO3(i)=fNO3(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'HO2 + NO = OH + NO2';
k(:,i) = 3.45e-12.*exp(270./T);
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'NO'; 
fHO2(i)=fHO2(i)-1; fNO(i)=fNO(i)-1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'HO2 + NO2 = HO2NO2';
k(:,i) = KMT09;
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'NO2'; 
fHO2(i)=fHO2(i)-1; fNO2(i)=fNO2(i)-1; fHO2NO2(i)=fHO2NO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + HO2NO2 = NO2';
k(:,i) = 3.2e-13.*exp(690./T).*1.0;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'HO2NO2'; 
fOH(i)=fOH(i)-1; fHO2NO2(i)=fHO2NO2(i)-1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'HO2 + NO3 = OH + NO2';
k(:,i) = 4.0e-12;
Gstr{i,1} = 'HO2'; Gstr{i,2} = 'NO3'; 
fHO2(i)=fHO2(i)-1; fNO3(i)=fNO3(i)-1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + HONO = NO2';
k(:,i) = 2.5e-12.*exp(260./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'HONO'; 
fOH(i)=fOH(i)-1; fHONO(i)=fHONO(i)-1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + HNO3 = NO3';
k(:,i) = KMT11;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'HNO3'; 
fOH(i)=fOH(i)-1; fHNO3(i)=fHNO3(i)-1; fNO3(i)=fNO3(i)+1; 

i=i+1;
Rnames{i} = 'O + SO2 = SO3';
k(:,i) = 4.0e-32.*exp(-1000./T).*M;
Gstr{i,1} = 'O'; Gstr{i,2} = 'SO2'; 
fO(i)=fO(i)-1; fSO2(i)=fSO2(i)-1; fSO3(i)=fSO3(i)+1; 

i=i+1;
Rnames{i} = 'OH + SO2 = HSO3';
k(:,i) = KMT12;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'SO2'; 
fOH(i)=fOH(i)-1; fSO2(i)=fSO2(i)-1; fHSO3(i)=fHSO3(i)+1; 

i=i+1;
Rnames{i} = 'HSO3 = HO2 + SO3';
k(:,i) = 1.3e-12.*exp(-330./T).*.21.*M;
Gstr{i,1} = 'HSO3'; 
fHSO3(i)=fHSO3(i)-1; fHO2(i)=fHO2(i)+1; fSO3(i)=fSO3(i)+1; 

i=i+1;
Rnames{i} = 'HNO3 = NA';
k(:,i) = 6.00e-06;
Gstr{i,1} = 'HNO3'; 
fHNO3(i)=fHNO3(i)-1; fNA(i)=fNA(i)+1; 

i=i+1;
Rnames{i} = 'N2O5 = NA + NA';
k(:,i) = 4.00e-04;
Gstr{i,1} = 'N2O5'; 
fN2O5(i)=fN2O5(i)-1; fNA(i)=fNA(i)+1; fNA(i)=fNA(i)+1; 

i=i+1;
Rnames{i} = 'SO3 = SA';
k(:,i) = 1.20e-15.*H2O;
Gstr{i,1} = 'SO3'; 
fSO3(i)=fSO3(i)-1; fSA(i)=fSA(i)+1; 

i=i+1;
Rnames{i} = 'O3 + hv = O1D';
k(:,i) = J1;
Gstr{i,1} = 'O3'; 
fO3(i)=fO3(i)-1; fO1D(i)=fO1D(i)+1; 

i=i+1;
Rnames{i} = 'O3 + hv = O';
k(:,i) = J2;
Gstr{i,1} = 'O3'; 
fO3(i)=fO3(i)-1; fO(i)=fO(i)+1; 

i=i+1;
Rnames{i} = 'H2O2 + hv = OH + OH';
k(:,i) = J3;
Gstr{i,1} = 'H2O2'; 
fH2O2(i)=fH2O2(i)-1; fOH(i)=fOH(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'NO2 + hv = NO + O';
k(:,i) = J4;
Gstr{i,1} = 'NO2'; 
fNO2(i)=fNO2(i)-1; fNO(i)=fNO(i)+1; fO(i)=fO(i)+1; 

i=i+1;
Rnames{i} = 'NO3 + hv = NO';
k(:,i) = J5;
Gstr{i,1} = 'NO3'; 
fNO3(i)=fNO3(i)-1; fNO(i)=fNO(i)+1; 

i=i+1;
Rnames{i} = 'NO3 + hv = NO2 + O';
k(:,i) = J6;
Gstr{i,1} = 'NO3'; 
fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fO(i)=fO(i)+1; 

i=i+1;
Rnames{i} = 'HONO + hv = OH + NO';
k(:,i) = J7;
Gstr{i,1} = 'HONO'; 
fHONO(i)=fHONO(i)-1; fOH(i)=fOH(i)+1; fNO(i)=fNO(i)+1; 

i=i+1;
Rnames{i} = 'HNO3 + hv = OH + NO2';
k(:,i) = J8;
Gstr{i,1} = 'HNO3'; 
fHNO3(i)=fHNO3(i)-1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'N2O5 = NO2 + NO3';
k(:,i) = KMT04;
Gstr{i,1} = 'N2O5'; 
fN2O5(i)=fN2O5(i)-1; fNO2(i)=fNO2(i)+1; fNO3(i)=fNO3(i)+1; 

i=i+1;
Rnames{i} = 'HO2NO2 = HO2 + NO2';
k(:,i) = KMT10;
Gstr{i,1} = 'HO2NO2'; 
fHO2NO2(i)=fHO2NO2(i)-1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + HO2 = CH3OOH';
k(:,i) = 3.8e-13.*exp(780./T).*(1-1./(1+498.*exp(-1160./T)));
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'HO2'; 
fCH3OO(i)=fCH3OO(i)-1; fHO2(i)=fHO2(i)-1; fCH3OOH(i)=fCH3OOH(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + HO2 = HCHO';
k(:,i) = 3.8e-13.*exp(780./T).*(1./(1+498.*exp(-1160./T)));
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'HO2'; 
fCH3OO(i)=fCH3OO(i)-1; fHO2(i)=fHO2(i)-1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + NO = CH3NO3';
k(:,i) = 2.3e-12.*exp(360./T).*0.001;
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'NO'; 
fCH3OO(i)=fCH3OO(i)-1; fNO(i)=fNO(i)-1; fCH3NO3(i)=fCH3NO3(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + NO = CH3O + NO2';
k(:,i) = 2.3e-12.*exp(360./T).*0.999;
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'NO'; 
fCH3OO(i)=fCH3OO(i)-1; fNO(i)=fNO(i)-1; fCH3O(i)=fCH3O(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + NO2 = CH3OONO2';
k(:,i) = KMT13;
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'NO2'; 
fCH3OO(i)=fCH3OO(i)-1; fNO2(i)=fNO2(i)-1; fCH3OONO2(i)=fCH3OONO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + NO3 = CH3O + NO2';
k(:,i) = 1.2e-12;
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'NO3'; 
fCH3OO(i)=fCH3OO(i)-1; fNO3(i)=fNO3(i)-1; fCH3O(i)=fCH3O(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OO + CH3OO = CH3O + CH3O';
k(:,i) = 2.*1.03e-13.*exp(365./T).*7.18.*exp(-885./T);
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'CH3OO';
fCH3OO(i)=fCH3OO(i)-2; fCH3O(i)=fCH3O(i)+2; 

i=i+1;
Rnames{i} = 'CH3OO + CH3OO = CH3OH + HCHO';
k(:,i) = 2.*1.03e-13.*exp(365./T).*(1-7.18.*exp(-885./T));
Gstr{i,1} = 'CH3OO'; Gstr{i,2} = 'CH3OO';
fCH3OO(i)=fCH3OO(i)-2; fCH3OH(i)=fCH3OH(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'CH3OOH + hv = CH3O + OH';
k(:,i) = J41;
Gstr{i,1} = 'CH3OOH'; 
fCH3OOH(i)=fCH3OOH(i)-1; fCH3O(i)=fCH3O(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'OH + CH3OOH = CH3OO';
k(:,i) = 5.3e-12.*exp(190./T).*0.6;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'CH3OOH'; 
fOH(i)=fOH(i)-1; fCH3OOH(i)=fCH3OOH(i)-1; fCH3OO(i)=fCH3OO(i)+1; 

i=i+1;
Rnames{i} = 'OH + CH3OOH = HCHO + OH';
k(:,i) = 5.3e-12.*exp(190./T).*0.4;
Gstr{i,1} = 'OH'; Gstr{i,2} = 'CH3OOH'; 
fOH(i)=fOH(i)-1; fCH3OOH(i)=fCH3OOH(i)-1; fHCHO(i)=fHCHO(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'HCHO + hv = CO + HO2 + HO2';
k(:,i) = J11;
Gstr{i,1} = 'HCHO'; 
fHCHO(i)=fHCHO(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HCHO + hv = H2 + CO';
k(:,i) = J12;
Gstr{i,1} = 'HCHO'; 
fHCHO(i)=fHCHO(i)-1; fH2(i)=fH2(i)+1; fCO(i)=fCO(i)+1; 

i=i+1;
Rnames{i} = 'NO3 + HCHO = HNO3 + CO + HO2';
k(:,i) = 5.5e-16;
Gstr{i,1} = 'NO3'; Gstr{i,2} = 'HCHO'; 
fNO3(i)=fNO3(i)-1; fHCHO(i)=fHCHO(i)-1; fHNO3(i)=fHNO3(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + HCHO = HO2 + CO';
k(:,i) = 5.4e-12.*exp(135./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'HCHO'; 
fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)-1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+1; 

i=i+1;
Rnames{i} = 'CH3NO3 + hv = CH3O + NO2';
k(:,i) = J51;
Gstr{i,1} = 'CH3NO3'; 
fCH3NO3(i)=fCH3NO3(i)-1; fCH3O(i)=fCH3O(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'OH + CH3NO3 = HCHO + NO2';
k(:,i) = 4.0e-13.*exp(-845./T);
Gstr{i,1} = 'OH'; Gstr{i,2} = 'CH3NO3'; 
fOH(i)=fOH(i)-1; fCH3NO3(i)=fCH3NO3(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3O = HCHO + HO2';
k(:,i) = 7.2e-14.*exp(-1080./T).*.21.*M;
Gstr{i,1} = 'CH3O'; 
fCH3O(i)=fCH3O(i)-1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OONO2 = CH3OO + NO2';
k(:,i) = KMT14;
Gstr{i,1} = 'CH3OONO2'; 
fCH3OONO2(i)=fCH3OONO2(i)-1; fCH3OO(i)=fCH3OO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3OH + OH = HO2 + HCHO';
k(:,i) = 2.85e-12.*exp(-345./T);
Gstr{i,1} = 'CH3OH'; Gstr{i,2} = 'OH'; 
fCH3OH(i)=fCH3OH(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

% CALTECH ISOPRENE REDUCED-PLUS SCHEME

i=i+1;
Rnames{i} = 'ISOP + OH = IHOO1';
k(:,i) = KIHOO1;
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'OH'; 
fISOP(i)=fISOP(i)-1; fOH(i)=fOH(i)-1; fIHOO1(i)=fIHOO1(i)+1;

i=i+1;
Rnames{i} = 'ISOP + OH = IHOO4';
k(:,i) = KIHOO4;
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'OH'; 
fISOP(i)=fISOP(i)-1; fOH(i)=fOH(i)-1; fIHOO4(i)=fIHOO4(i)+1;

i=i+1;
Rnames{i} = 'ISOP + OH = 0.15ISOP1CO2OOH + 0.25ISOP1CO4OOH + 0.4HO2 + 0.6CO + 1.5OH + 0.3HCHO + 0.3MGLY + 0.3HPETHNL + 0.3CH3CO3';
k(:,i) = KISO1;
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'OH'; 
fISOP(i)=fISOP(i)-1; fOH(i)=fOH(i)+0.5; fISOP1CO2OOH(i)=fISOP1CO2OOH(i)+0.15; fISOP1CO4OOH(i)=fISOP1CO4OOH(i)+0.25; fHO2(i)=fHO2(i)+0.4; fCO(i)=fCO(i)+0.6; fHCHO(i)=fHCHO(i)+0.3; fMGLY(i)=fMGLY(i)+0.3; fHPETHNL(i)=fHPETHNL(i)+0.3; fCH3CO3(i)=fCH3CO3(i)+0.3; 

i=i+1;
Rnames{i} = 'ISOP + OH = 0.15ISOP3OOH4CO + 0.25ISOP1OOH4CO + 1.5OH + 0.3HCHO + 0.9CO + 0.7HO2 + 0.3MGLY + 0.3HPAC';
k(:,i) = KISO4;
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'OH'; 
fISOP(i)=fISOP(i)-1; fOH(i)=fOH(i)+0.5; fISOP3OOH4CO(i)=fISOP3OOH4CO(i)+0.15; fISOP1OOH4CO(i)=fISOP1OOH4CO(i)+0.25; fHO2(i)=fHO2(i)+0.7; fCO(i)=fCO(i)+0.9; fHCHO(i)=fHCHO(i)+0.3; fMGLY(i)=fMGLY(i)+0.3; fHPAC(i)=fHPAC(i)+0.3;

i=i+1;
Rnames{i} = 'IHOO1 + IHOO1 = 2MVK + 2HO2 + 2HCHO';
k(:,i) = (1.1644-T.*7.0485E-4)*6.92E-14;
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'IHOO1'; 
fIHOO1(i)=fIHOO1(i)-2; fMVK(i)=fMVK(i)+2; fHCHO(i)=fHCHO(i)+2; fHO2(i)=fHO2(i)+2;

i=i+1;
Rnames{i} = 'IHOO4 + IHOO4 = 2MACR + 2HO2 + 2HCHO';
k(:,i) = (1.2038-T.*9.0435E-4)*5.74E-12; 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-2; fMACR(i)=fMACR(i)+2; fHCHO(i)=fHCHO(i)+2; fHO2(i)=fHO2(i)+2;

i=i+1;
Rnames{i} = 'IHOO1 + IHOO4 = MACR + MVK + 2HO2 + 2HCHO';
k(:,i) = (2.3682-T.*1.6092E-3)/2*3.08E-12; 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-1; fIHOO1(i)=fIHOO1(i)-1; fMVK(i)=fMVK(i)+1; fMACR(i)=fMACR(i)+1; fHCHO(i)=fHCHO(i)+2; fHO2(i)=fHO2(i)+2;

i=i+1;
Rnames{i} = 'IHOO1 + IHOO1 = 0.5HO2 + 0.5ISOP1CO4OH + 0.5CO + 0.5OH + 0.5MVK3OOH4OH';
k(:,i) = (T.*7.0485E-4-0.1644)*2.49E-12; 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'IHOO1'; 
fIHOO1(i)=fIHOO1(i)-2; fHO2(i)=fHO2(i)+0.5; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.5; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.5; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5; 

i=i+1;
Rnames{i} = 'IHOO4 + IHOO4 = 0.5HO2 + 0.5ISOP1OH4CO + 0.5CO + 0.5OH + 0.5MACR2OOH3OH';
k(:,i) = (T.*9.0435E-4-0.2038)*3.94E-12; 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-2; fHO2(i)=fHO2(i)+0.5; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.5; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.5; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5; 

i=i+1;
Rnames{i} = 'IHOO1 + IHOO4 = 0.5HO2 + 0.25ISOP1CO4OH + 0.25ISOP1OH4CO + 0.5CO + 0.5OH + 0.25MVK3OOH4OH + 0.25MACR2OOH3OH';
k(:,i) = (T.*1.6092E-3-0.3682)/2*3.08E-12; 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-1; fIHOO1(i)=fIHOO1(i)-1; fHO2(i)=fHO2(i)+0.5; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.25; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.25; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.25; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.25;

i=i+1;
Rnames{i} = 'IHOO1 + NO = ISOP1OH2N';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.19,6,1.1644,7.05E-4);
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'NO'; 
fIHOO1(i)=fIHOO1(i)-1; fNO(i)=fNO(i)-1; fISOP1OH2N(i)=fISOP1OH2N(i)+1; 

i=i+1;
Rnames{i} = 'IHOO1 + NO = NO2 + MVK + HO2 + HCHO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.19,6,1.1644,7.05E-4);
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'NO'; 
fIHOO1(i)=fIHOO1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fMVK(i)=fMVK(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'IHOO1 + HO2 = 0.063MVK + 0.063OH + 0.063HO2 + 0.063HCHO + 0.937ISOP1OH2OOH';
k(:,i) = (1.1644-T.*7.0485E-4).*2.12E-13.*exp(1300./T);
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'HO2'; 
fIHOO1(i)=fIHOO1(i)-1; fHO2(i)=fHO2(i)-0.937; fISOP1OH2OOH(i)=fISOP1OH2OOH(i)+0.937; fHCHO(i)=fHCHO(i)+0.063; fMVK(i)=fMVK(i)+0.063; fOH(i)=fOH(i)+0.063;

i=i+1;
Rnames{i} = 'IHOO1 + NO = ISOP1OH4N';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.421,6,-0.1644,-7.05E-4);
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'NO'; 
fIHOO1(i)=fIHOO1(i)-1; fNO(i)=fNO(i)-1; fISOP1OH4N(i)=fISOP1OH4N(i)+1; 

i=i+1;
Rnames{i} = 'IHOO1 + NO = NO2 + 0.45ISOP1CO4OH + 0.45HO2 + 0.55MVK3OOH4OH + 0.55CO + 0.55OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.421,6,-0.1644,-7.05E-4);
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'NO'; 
fIHOO1(i)=fIHOO1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.45; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.45; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.55; fOH(i)=fOH(i)+0.55; fCO(i)=fCO(i)+0.55; 

i=i+1;
Rnames{i} = 'IHOO1 + HO2 = ISOP1OH4OOH';
k(:,i) = (T.*7.0485E-4-0.1644).*2.12E-13.*exp(1300./T); 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'HO2'; 
fIHOO1(i)=fIHOO1(i)-1; fHO2(i)=fHO2(i)-0.937; fISOP1OH4OOH(i)=fISOP1OH4OOH(i)+1; 

i=i+1;
Rnames{i} = 'IHOO1 + CH3OO = MVK + 2HO2 + 2HCHO';
k(:,i) = (1.1644-T.*7.0485E-4)*2.00E-12 ; 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'CH3OO'; 
fIHOO1(i)=fIHOO1(i)-1; fCH3OO(i)=fCH3OO(i)-1; fMVK(i)=fMVK(i)+1; fHO2(i)=fHO2(i)+2; fHCHO(i)=fHCHO(i)+2; 

i=i+1;
Rnames{i} = 'IHOO1 + CH3OO = HCHO + 0.5ISOP1CO4OH + 1.5HO2 + 0.5MVK3OOH4OH + 0.5CO + 0.5OH';
k(:,i) = (T.*7.0485E-4-0.1644)*2.00E-12 ; 
Gstr{i,1} = 'IHOO1'; Gstr{i,2} = 'CH3OO'; 
fIHOO1(i)=fIHOO1(i)-1; fCH3OO(i)=fCH3OO(i)-1; fHCHO(i)=fHCHO(i)+1; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.5; fHO2(i)=fHO2(i)+1.5; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.5; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5; 

i=i+1;
Rnames{i} = 'IHOO1 = HCHO + OH + MVK';
k(:,i) = (1.1644-T.*7.0485E-4).*1.04E11.*exp(-9746./T); 
Gstr{i,1} = 'IHOO1'; 
fIHOO1(i)=fIHOO1(i)-1; fHCHO(i)=fHCHO(i)+1; fOH(i)=fOH(i)+1; fMVK(i)=fMVK(i)+1; 

i=i+1;
Rnames{i} = 'IHOO1 = 0.15ISOP1CO2OOH + 0.25ISOP1CO4OOH + 0.4HO2 + 0.6CO + 1.5OH + 0.3HCHO + 0.3MGLY + 0.3HPETHNL + 0.3CH3CO3';
k(:,i) = (T.*5.1242E-5-0.0128).*5.05E15.*exp(-12200./T).*exp(1E8./T.^3); 
Gstr{i,1} = 'IHOO1'; 
fIHOO1(i)=fIHOO1(i)-1; fISOP1CO2OOH(i)=fISOP1CO2OOH(i)+0.15; fISOP1CO4OOH(i)=fISOP1CO4OOH(i)+0.25; fHO2(i)=fHO2(i)+0.4; fCO(i)=fCO(i)+0.6; fOH(i)=fOH(i)+1.5; fHCHO(i)=fHCHO(i)+0.3; fMGLY(i)=fMGLY(i)+0.3; fHPETHNL(i)=fHPETHNL(i)+0.3; fCH3CO3(i)=fCH3CO3(i)+0.3; 

i=i+1;
Rnames{i} = 'IHOO4 + NO = ISOP3N4OH';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.297,6,1.2038,9.04E-4);
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'NO'; 
fIHOO4(i)=fIHOO4(i)-1; fNO(i)=fNO(i)-1; fISOP3N4OH(i)=fISOP3N4OH(i)+1; 

i=i+1;
Rnames{i} = 'IHOO4 + NO = NO2 + MACR + HO2 + HCHO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.297,6,1.2038,9.04E-4);
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'NO'; 
fIHOO4(i)=fIHOO4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMACR(i)=fMACR(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'IHOO4 + HO2 = 0.063MACR + 0.063OH + 0.063HO2 + 0.063HCHO + 0.937ISOP3OOH4OH';
k(:,i) = (1.2038-T.*9.0435E-4).*2.12E-13.*exp(1300./T); 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'HO2'; 
fIHOO4(i)=fIHOO4(i)-1; fHO2(i)=fHO2(i)-0.937; fMACR(i)=fMACR(i)+0.063; fOH(i)=fOH(i)+0.063; fHCHO(i)=fHCHO(i)+0.063; fISOP3OOH4OH(i)=fISOP3OOH4OH(i)+0.937; 

i=i+1;
Rnames{i} = 'IHOO4 + NO = ISOP1N4OH';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.421,6,-0.2038,-9.04E-4);
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'NO'; 
fIHOO4(i)=fIHOO4(i)-1; fNO(i)=fNO(i)-1; fISOP1N4OH(i)=fISOP1N4OH(i)+1;

i=i+1;
Rnames{i} = 'IHOO4 + NO = NO2 + 0.45HO2 + 0.45ISOP1OH4CO + 0.55MACR2OOH3OH + 0.55CO + 0.55OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.421,6,-0.2038,-9.04E-4);
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'NO'; 
fIHOO4(i)=fIHOO4(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.45; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.45; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.55; fCO(i)=fCO(i)+0.55; fOH(i)=fOH(i)+0.55; 

i=i+1;
Rnames{i} = 'IHOO4 + HO2 = ISOP1OOH4OH';
k(:,i) = (T.*9.0435E-4-0.2038).*2.12E-13.*exp(1300./T); 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'HO2'; 
fIHOO4(i)=fIHOO4(i)-1; fHO2(i)=fHO2(i)-1; fISOP1OOH4OH(i)=fISOP1OOH4OH(i)+1; 

i=i+1;
Rnames{i} = 'IHOO4 + CH3OO = MACR + 2HO2 + 2HCHO';
k(:,i) = (1.2038-T.*9.0435E-4)*2.00E-12 ; 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'CH3OO'; 
fIHOO4(i)=fIHOO4(i)-1; fCH3OO(i)=fCH3OO(i)-1; fMACR(i)=fMACR(i)+1; fHO2(i)=fHO2(i)+2; fHCHO(i)=fHCHO(i)+2; 

i=i+1;
Rnames{i} = 'IHOO4 + CH3OO = HCHO + 0.5ISOP1OH4CO + 1.5HO2 + 0.5MACR2OOH3OH + 0.5CO + 0.5OH';
k(:,i) = (T.*9.0435E-4-0.2038)*2.00E-12 ; 
Gstr{i,1} = 'IHOO4'; Gstr{i,2} = 'CH3OO'; 
fIHOO4(i)=fIHOO4(i)-1; fCH3OO(i)=fCH3OO(i)-1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1.5; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.5; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+1; fCO(i)=fCO(i)+0.5; fOH(i)=fOH(i)+0.5; 

i=i+1;
Rnames{i} = 'IHOO4 = MACR + OH + HCHO';
k(:,i) = (1.2038-T.*9.0435E-4).*1.88E11.*exp(-9752./T); 
Gstr{i,1} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-1; fMACR(i)=fMACR(i)+1; fOH(i)=fOH(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'IHOO4 = 0.15ISOP3OOH4CO + 0.25ISOP1OOH4CO + 1.5OH + 0.3HCHO + 0.9CO + 0.7HO2 + 0.3MGLY + 0.3HPAC';
k(:,i) = (T.*1.1346E-4-0.0306).*2.22E9.*exp(-7160./T).*exp(1E8./T.^3); 
Gstr{i,1} = 'IHOO4'; 
fIHOO4(i)=fIHOO4(i)-1; fISOP3OOH4CO(i)=fISOP3OOH4CO(i)+0.15; fISOP1OOH4CO(i)=fISOP1OOH4CO(i)+0.25; fOH(i)=fOH(i)+1.5; fHCHO(i)=fHCHO(i)+0.3; fCO(i)=fCO(i)+0.9; fHO2(i)=fHO2(i)+0.7; fMGLY(i)=fMGLY(i)+0.3; fHPAC(i)=fHPAC(i)+0.3; 

i=i+1;
Rnames{i} = 'ISOP1CO2OOH + OH = OH + 0.230MVK + 0.420CO + 0.190MVK3OH4OOH + 0.580ICHE';
k(:,i) = 2.2E-11*exp(390./T) ;
Gstr{i,1} = 'ISOP1CO2OOH'; Gstr{i,2} = 'OH'; 
fISOP1CO2OOH(i)=fISOP1CO2OOH(i)-1; fMVK(i)=fMVK(i)+0.23; fCO(i)=fCO(i)+0.42; fMVK3OH4OOH(i)=fMVK3OH4OOH(i)+0.19; fICHE(i)=fICHE(i)+0.58; 

i=i+1;
Rnames{i} = 'ISOP1CO2OOH = CO + OH + HO2 + MVK';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'ISOP1CO2OOH'; 
fISOP1CO2OOH(i)=fISOP1CO2OOH(i)-1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fMVK(i)=fMVK(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3OOH4CO + OH = OH + 0.770ICHE + 0.230CO + 0.090MACR2OH3OOH + 0.140MACR';
k(:,i) = 3.5E-11*exp(390./T) ;
Gstr{i,1} = 'ISOP3OOH4CO'; Gstr{i,2} = 'OH'; 
fISOP3OOH4CO(i)=fISOP3OOH4CO(i)-1; fICHE(i)=fICHE(i)+0.77; fCO(i)=fCO(i)+0.23; fMACR2OH3OOH(i)=fMACR2OH3OOH(i)+0.09; fMACR(i)=fMACR(i)+0.14;

i=i+1;
Rnames{i} = 'ISOP3OOH4CO = CO + OH + HO2 + MACR';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'ISOP3OOH4CO'; 
fISOP3OOH4CO(i)=fISOP3OOH4CO(i)-1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fMACR(i)=fMACR(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1CO4OH + OH = 1.065OH + 0.355CO2 + 0.355MGLY + 0.548CO + 0.193HO2 + 0.193MVK3OOH4OH + 0.452IEPOXAOO';
k(:,i) = 4.64E-12*exp(650./T);
Gstr{i,1} = 'ISOP1CO4OH'; Gstr{i,2} = 'OH'; 
fISOP1CO4OH(i)=fISOP1CO4OH(i)-1; fOH(i)=fOH(i)+0.065; fMGLY(i)=fMGLY(i)+0.355; fCO(i)=fCO(i)+0.548; fHO2(i)=fHO2(i)+0.193; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.193; fIEPOXAOO(i)=fIEPOXAOO(i)+0.452; 

i=i+1;
Rnames{i} = 'ISOP1OH4CO + OH = 1.065OH + 0.355CO2 + 0.355MGLY + 0.452HO2 + 0.807CO + 0.452MACR2OOH3OH + 0.193IEPOXBOO';
k(:,i) = 4.64E-12*exp(650./T);
Gstr{i,1} = 'ISOP1OH4CO'; Gstr{i,2} = 'OH'; 
fISOP1OH4CO(i)=fISOP1OH4CO(i)-1; fOH(i)=fOH(i)+0.065; fMGLY(i)=fMGLY(i)+0.355; fCO(i)=fCO(i)+0.807; fHO2(i)=fHO2(i)+0.452; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.452; fIEPOXBOO(i)=fIEPOXBOO(i)+0.193; 

i=i+1;
Rnames{i} = 'ISOP3CO4OH + OH = ICHOO';
k(:,i) = 2.7E-11*exp(390./T); 
Gstr{i,1} = 'ISOP3CO4OH'; Gstr{i,2} = 'OH'; 
fISOP3CO4OH(i)=fISOP3CO4OH(i)-1; fOH(i)=fOH(i)-1; fICHOO(i)=fICHOO(i)+1;

i=i+1;
Rnames{i} = 'ICHOO + HO2 = 0.35ICPDH + 0.65OH + 0.13MVK3CO4OH + 0.65HCHO + 0.65HO2 + 0.52HAC + 0.52CO';
k(:,i) = 2.38E-13*exp(1300./T); 
Gstr{i,1} = 'ICHOO'; Gstr{i,2} = 'HO2'; 
fICHOO(i)=fICHOO(i)-1; fHO2(i)=fHO2(i)-0.35; fICPDH(i)=fICPDH(i)+0.35; fOH(i)=fOH(i)+0.65; fMVK3CO4OH(i)=fMVK3CO4OH(i)+0.13; fHCHO(i)=fHCHO(i)+0.65; fHAC(i)=fHAC(i)+0.52; fCO(i)=fCO(i)+0.52; 

i=i+1;
Rnames{i} = 'ICHOO + NO = IDHCN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,13.098,8,1,0);
Gstr{i,1} = 'ICHOO'; Gstr{i,2} = 'NO'; 
fICHOO(i)=fICHOO(i)-1; fNO(i)=fNO(i)-1; fIDHCN(i)=fIDHCN(i)+1; 

i=i+1;
Rnames{i} = 'ICHOO + NO = NO2 + 0.8HAC + 0.8CO + HCHO + HO2 + 0.2MVK3CO4OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,13.098,8,1,0.);
Gstr{i,1} = 'ICHOO'; Gstr{i,2} = 'NO'; 
fICHOO(i)=fICHOO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHAC(i)=fHAC(i)+0.8; fCO(i)=fCO(i)+0.8; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; fMVK3CO4OH(i)=fMVK3CO4OH(i)+0.2; 

i=i+1;
Rnames{i} = 'ICHOO = HO2 + CO + CO + HAC + OH';
k(:,i) = 1.875E13*exp(-10000./T);
Gstr{i,1} = 'ICHOO'; 
fICHOO(i)=fICHOO(i)-1; fCO(i)=fCO(i)+2; fHO2(i)=fHO2(i)+1; fHAC(i)=fHAC(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH2OOH + OH = 0.655IHPOO3 + 0.345IHPOO1';
k(:,i) = 2.47E-12*exp(390./T); 
Gstr{i,1} = 'ISOP1OH2OOH'; Gstr{i,2} = 'OH'; 
fISOP1OH2OOH(i)=fISOP1OH2OOH(i)-1; fOH(i)=fOH(i)-1; fIHPOO3(i)=fIHPOO3(i)+0.655; fIHPOO1(i)=fIHPOO1(i)+0.345; 

i=i+1;
Rnames{i} = 'ISOP1OH2OOH + OH = OH + 0.67IEPOXt + 0.33IEPOXc';
k(:,i) = F0AM_isop_EPO(T,M,1.62E-11,390,4.77E-21);
Gstr{i,1} = 'ISOP1OH2OOH'; Gstr{i,2} = 'OH'; 
fISOP1OH2OOH(i)=fISOP1OH2OOH(i)-1; fIEPOXt(i)=fIEPOXt(i)+0.67; fIEPOXc(i)=fIEPOXc(i)+0.33; 

i=i+1;
Rnames{i} = 'ISOP3OOH4OH + OH = 0.655IHPOO3 + 0.345IHPOO2';
k(:,i) = 4.35E-12*exp(390./T); 
Gstr{i,1} = 'ISOP3OOH4OH'; Gstr{i,2} = 'OH'; 
fISOP3OOH4OH(i)=fISOP3OOH4OH(i)-1; fOH(i)=fOH(i)-1; fIHPOO3(i)=fIHPOO3(i)+0.655; fIHPOO2(i)=fIHPOO2(i)+0.345; 

i=i+1;
Rnames{i} = 'ISOP3OOH4OH + OH = OH + 0.68IEPOXt + 0.32IEPOXc';
k(:,i) = F0AM_isop_EPO(T,M,2.85E-11,390,4.77E-21);
Gstr{i,1} = 'ISOP3OOH4OH'; Gstr{i,2} = 'OH'; 
fISOP3OOH4OH(i)=fISOP3OOH4OH(i)-1; fIEPOXt(i)=fIEPOXt(i)+0.68; fIEPOXc(i)=fIEPOXc(i)+0.32; 

i=i+1;
Rnames{i} = 'ISOP1OH2OOH + OH = 0.75IHOO1 + 0.125MVK + 0.125MVK3OOH4OH + 0.25HO2 + 0.25CO';
k(:,i) = 6.1E-12*exp(200./T); 
Gstr{i,1} = 'ISOP1OH2OOH'; Gstr{i,2} = 'OH'; 
fISOP1OH2OOH(i)=fISOP1OH2OOH(i)-1; fOH(i)=fOH(i)-1; fIHOO1(i)=fIHOO1(i)+0.75; fMVK(i)=fMVK(i)+0.125; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.125; fHO2(i)=fHO2(i)+0.25; fCO(i)=fCO(i)+0.25; 

i=i+1;
Rnames{i} = 'ISOP3OOH4OH + OH = 0.51IHOO4 + 0.16OH + 0.16ISOP3CO4OH + 0.33CO + 0.33HO2 + 0.165MACR + 0.165MACR2OOH3OH';
k(:,i) = 4.1E-12*exp(200./T); 
Gstr{i,1} = 'ISOP3OOH4OH'; Gstr{i,2} = 'OH'; 
fISOP3OOH4OH(i)=fISOP3OOH4OH(i)-1; fOH(i)=fOH(i)-0.84; fIHOO4(i)=fIHOO4(i)+0.51; fISOP3CO4OH(i)=fISOP3CO4OH(i)+0.16; fCO(i)=fCO(i)+0.33; fHO2(i)=fHO2(i)+0.33; fMACR(i)=fMACR(i)+0.165; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.165; 

i=i+1;
Rnames{i} = 'ISOP1OH2OOH = MVK + HCHO + HO2 + OH';
k(:,i) = J41 ; % SUN*6.5E-6; 
Gstr{i,1} = 'ISOP1OH2OOH'; 
fISOP1OH2OOH(i)=fISOP1OH2OOH(i)-1; fMVK(i)=fMVK(i)+1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3OOH4OH = MACR + HCHO + HO2 + OH';
k(:,i) = J41 ; % SUN*6.5E-6; 
Gstr{i,1} = 'ISOP3OOH4OH'; 
fISOP3OOH4OH(i)=fISOP3OOH4OH(i)-1; fMACR(i)=fMACR(i)+1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH4OOH + OH = 0.595IHPOO1 + 0.03IHOO1 + 0.06ISOP1OH4CO + 0.024HO2 + 0.009ISOP1CO2OOH + 0.015ISOP1CO4OOH + 0.405OH + 0.036CO + 0.018HCHO + 0.018MGLY + 0.018HPETHNL + 0.018CH3CO3 + 0.255IEPOXD';
k(:,i) = 3.53E-11*exp(390./T);
Gstr{i,1} = 'ISOP1OH4OOH'; Gstr{i,2} = 'OH'; 
fISOP1OH4OOH(i)=fISOP1OH4OOH(i)-1; fOH(i)=fOH(i)-0.595; fIHPOO1(i)=fIHPOO1(i)+0.595; fIHOO1(i)=fIHOO1(i)+0.03; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.06; fHO2(i)=fHO2(i)+0.024; fISOP1CO2OOH(i)=fISOP1CO2OOH(i)+0.009; fISOP1CO4OOH(i)=fISOP1CO4OOH(i)+0.015; fCO(i)=fCO(i)+0.036; fHCHO(i)=fHCHO(i)+0.018; fMGLY(i)=fMGLY(i)+0.018; fHPETHNL(i)=fHPETHNL(i)+0.018; fCH3CO3(i)=fCH3CO3(i)+0.018; fIEPOXD(i)=fIEPOXD(i)+0.255; 

i=i+1;
Rnames{i} = 'ISOP1OOH4OH + OH = 0.255IHPOO2 + 0.03IHOO4 + 0.745OH + 0.06ISOP1CO4OH + 0.009ISOP3OOH4CO + 0.015ISOP1OOH4CO + 0.042HO2 + 0.018HCHO + 0.054CO + 0.018MGLY + 0.018HPAC + 0.595IEPOXD';
k(:,i) = 3.53E-11*exp(390./T);
Gstr{i,1} = 'ISOP1OOH4OH'; Gstr{i,2} = 'OH'; 
fISOP1OOH4OH(i)=fISOP1OOH4OH(i)-1; fOH(i)=fOH(i)-0.255; fIHPOO2(i)=fIHPOO2(i)+0.255; fIHOO4(i)=fIHOO4(i)+0.03; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.06; fISOP3OOH4CO(i)=fISOP3OOH4CO(i)+0.009; fHO2(i)=fHO2(i)+0.042; fHCHO(i)=fHCHO(i)+0.018; fCO(i)=fCO(i)+0.054; fMGLY(i)=fMGLY(i)+0.018; fHPAC(i)=fHPAC(i)+0.018; fISOP1OOH4CO(i)=fISOP1OOH4CO(i)+0.015; fIEPOXD(i)=fIEPOXD(i)+0.595; 

i=i+1;
Rnames{i} = 'IHPOO1 = IDHPE + OH';
k(:,i) = F0AM_isop_TUN(T,M,6.80E12,1.12E4,8.46E7); 
Gstr{i,1} = 'IHPOO1'; 
fIHPOO1(i)=fIHPOO1(i)-1; fIDHPE(i)=fIDHPE(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO1 + NO = 0.17MACR2OOH3OH + 0.17HCHO + 0.83HPETHNL + 0.83HAC + NO2 + HO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,2.1,9,1,0);
Gstr{i,1} = 'IHPOO1'; Gstr{i,2} = 'NO'; 
fIHPOO1(i)=fIHPOO1(i)-1; fNO(i)=fNO(i)-1; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.17; fHCHO(i)=fHCHO(i)+0.17; fHPETHNL(i)=fHPETHNL(i)+0.83; fHAC(i)=fHAC(i)+0.83; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO1 + NO = IDHPN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,2.1,9,1,0.);
Gstr{i,1} = 'IHPOO1'; Gstr{i,2} = 'NO'; 
fIHPOO1(i)=fIHPOO1(i)-1; fNO(i)=fNO(i)-1; fIDHPN(i)=fIDHPN(i)+1;

i=i+1;
Rnames{i} = 'IHPOO1 + HO2 = 0.59IDHDP + 0.032MACR2OOH3OH + 0.032HCHO + 0.378HPETHNL + 0.378HAC + 0.41OH + 0.41HO2';
k(:,i) = 2.47E-13*exp(1300./T); 
Gstr{i,1} = 'IHPOO1'; Gstr{i,2} = 'HO2'; 
fIHPOO1(i)=fIHPOO1(i)-1; fHO2(i)=fHO2(i)-0.59; fIDHDP(i)=fIDHDP(i)+0.59; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.032; fOH(i)=fOH(i)+0.41; fHCHO(i)=fHCHO(i)+0.032; fHPETHNL(i)=fHPETHNL(i)+0.378; fHAC(i)=fHAC(i)+0.378; 

i=i+1;
Rnames{i} = 'IHPOO2 = 0.32ICPDH + 0.68IDHPE + OH';
k(:,i) = F0AM_isop_TUN(T,M,3.73E12,1.88E4,1.82E8); 
Gstr{i,1} = 'IHPOO2'; 
fIHPOO2(i)=fIHPOO2(i)-1; fICPDH(i)=fICPDH(i)+0.32; fIDHPE(i)=fIDHPE(i)+0.68; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO2 + NO = 0.847MVK3OOH4OH + 0.847HCHO + 0.153GLYC + 0.153HPAC + NO2 + HO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,2.315,9,1,0);
Gstr{i,1} = 'IHPOO2'; Gstr{i,2} = 'NO'; 
fIHPOO2(i)=fIHPOO2(i)-1; fNO(i)=fNO(i)-1; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.847; fHCHO(i)=fHCHO(i)+0.847; fGLYC(i)=fGLYC(i)+0.153; fHPAC(i)=fHPAC(i)+0.153; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO2 + NO = IDHPN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,2.315,9,1,0);
Gstr{i,1} = 'IHPOO2'; Gstr{i,2} = 'NO'; 
fIHPOO2(i)=fIHPOO2(i)-1; fNO(i)=fNO(i)-1; fIDHPN(i)=fIDHPN(i)+1;

i=i+1;
Rnames{i} = 'IHPOO2 + HO2 = 0.76IDHDP + 0.17MVK3OOH4OH + 0.17HCHO + 0.07GLYC + 0.07HPAC + 0.24OH + 0.24HO2';
k(:,i) = 2.47E-13*exp(1300./T); 
Gstr{i,1} = 'IHPOO2'; Gstr{i,2} = 'HO2'; 
fIHPOO2(i)=fIHPOO2(i)-1; fHO2(i)=fHO2(i)-0.76; fIDHDP(i)=fIDHDP(i)+0.76; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.17; fHCHO(i)=fHCHO(i)+0.17; fGLYC(i)=fGLYC(i)+0.07; fHPAC(i)=fHPAC(i)+0.07; fOH(i)=fOH(i)+0.24; 

i=i+1;
Rnames{i} = 'IHPOO3 = IDHPE + OH';
k(:,i) = F0AM_isop_TUN(T,M,1.87E12,9.63E3,8.02E7); 
Gstr{i,1} = 'IHPOO3'; 
fIHPOO3(i)=fIHPOO3(i)-1; fIDHPE(i)=fIDHPE(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO3 + NO = GLYC + HAC + NO2 + OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,3.079,9,1,0);
Gstr{i,1} = 'IHPOO3'; Gstr{i,2} = 'NO'; 
fIHPOO3(i)=fIHPOO3(i)-1; fNO(i)=fNO(i)-1; fHAC(i)=fHAC(i)+1; fGLYC(i)=fGLYC(i)+1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'IHPOO3 + NO = IDHPN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,3.079,9,1,0);
Gstr{i,1} = 'IHPOO3'; Gstr{i,2} = 'NO'; 
fIHPOO3(i)=fIHPOO3(i)-1; fNO(i)=fNO(i)-1; fIDHPN(i)=fIDHPN(i)+1;

i=i+1;
Rnames{i} = 'IHPOO3 + HO2 = 0.35IDHDP + 0.65GLYC + 0.65HAC + 1.3OH';
k(:,i) = 2.47E-13*exp(1300./T); 
Gstr{i,1} = 'IHPOO3'; Gstr{i,2} = 'HO2'; 
fIHPOO3(i)=fIHPOO3(i)-1; fHO2(i)=fHO2(i)-1; fIDHDP(i)=fIDHDP(i)+0.35; fGLYC(i)=fGLYC(i)+0.65; fHAC(i)=fHAC(i)+0.65; fOH(i)=fOH(i)+1.3; 

i=i+1;
Rnames{i} = 'IEPOXD + OH = 0.75ICHE + 0.75HO2 + 0.25ICHOO';
k(:,i) = 3.22E-11*exp(-400./T);
Gstr{i,1} = 'IEPOXD'; Gstr{i,2} = 'OH'; 
fIEPOXD(i)=fIEPOXD(i)-1; fOH(i)=fOH(i)-1; fICHE(i)=fICHE(i)+0.75; fHO2(i)=fHO2(i)+0.75; fICHOO(i)=fICHOO(i)+0.25; 

i=i+1;
Rnames{i} = 'IEPOXt + OH = ICHE + HO2';
k(:,i) = 1.05E-11*exp(-400./T); 
Gstr{i,1} = 'IEPOXt'; Gstr{i,2} = 'OH'; 
fIEPOXt(i)=fIEPOXt(i)-1; fOH(i)=fOH(i)-1; fICHE(i)=fICHE(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXt + OH = 0.67IEPOXAOO + 0.33IEPOXBOO';
k(:,i) = F0AM_isop_EPO(T,M,5.82E-11,-400,1.14E-20);
Gstr{i,1} = 'IEPOXt'; Gstr{i,2} = 'OH'; 
fIEPOXt(i)=fIEPOXt(i)-1; fOH(i)=fOH(i)-1; fIEPOXAOO(i)=fIEPOXAOO(i)+0.67; fIEPOXBOO(i)=fIEPOXBOO(i)+0.33; 

i=i+1;
Rnames{i} = 'IEPOXc + OH = ICHE + HO2';
k(:,i) = 8.25E-12*exp(-400./T); 
Gstr{i,1} = 'IEPOXc'; Gstr{i,2} = 'OH'; 
fIEPOXc(i)=fIEPOXc(i)-1; fOH(i)=fOH(i)-1; fICHE(i)=fICHE(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXc + OH = 0.81IEPOXAOO + 0.19IEPOXBOO';
k(:,i) = F0AM_isop_EPO(T,M,3.75E-11,-400,8.91E-21);
Gstr{i,1} = 'IEPOXc'; Gstr{i,2} = 'OH'; 
fIEPOXc(i)=fIEPOXc(i)-1; fOH(i)=fOH(i)-1; fIEPOXAOO(i)=fIEPOXAOO(i)+0.81; fIEPOXBOO(i)=fIEPOXBOO(i)+0.19; 

i=i+1;
Rnames{i} = 'IEPOXBOO = IDCHP + HO2';
k(:,i) = 1.875E13*exp(-10000./T);
Gstr{i,1} = 'IEPOXBOO'; 
fIEPOXBOO(i)=fIEPOXBOO(i)-1; fIDCHP(i)=fIDCHP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXBOO = CO + OH + MACR2OH3OH';
k(:,i) = 1.0E7*exp(-5000./T);
Gstr{i,1} = 'IEPOXBOO';
fIEPOXBOO(i)=fIEPOXBOO(i)-1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+1; fMACR2OH3OH(i)=fMACR2OH3OH(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXBOO + NO = NO2 + HO2 + 0.8GLYX + 0.8HAC + 0.2CO + 0.2MACR2OH3OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,16.463,8,1,0);
Gstr{i,1} = 'IEPOXBOO'; Gstr{i,2} = 'NO'; 
fIEPOXBOO(i)=fIEPOXBOO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fGLYX(i)=fGLYX(i)+0.8; fHAC(i)=fHAC(i)+0.8; fCO(i)=fCO(i)+0.2; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.2; 

i=i+1;
Rnames{i} = 'IEPOXBOO + NO = IDHCN';
k(:,i) =  F0AM_isop_NIT(T,M,2.7E-12,350,16.463,8,1,0);
Gstr{i,1} = 'IEPOXBOO'; Gstr{i,2} = 'NO'; 
fIEPOXBOO(i)=fIEPOXBOO(i)-1; fNO(i)=fNO(i)-1; fIDHCN(i)=fIDHCN(i)+1;

i=i+1;
Rnames{i} = 'IEPOXBOO + HO2 = 0.13CO + 0.65OH + 0.65HO2 + 0.13MACR2OH3OH + 0.52HAC + 0.52GLYC + 0.35ICPDH';
k(:,i) = 2.38E-13*exp(1300./T);
Gstr{i,1} = 'IEPOXBOO'; Gstr{i,2} = 'HO2'; 
fIEPOXBOO(i)=fIEPOXBOO(i)-1; fHO2(i)=fHO2(i)-0.35; fCO(i)=fCO(i)+0.13; fOH(i)=fOH(i)+0.65; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.13; fHAC(i)=fHAC(i)+0.52; fGLYC(i)=fGLYC(i)+0.52; fICPDH(i)=fICPDH(i)+0.35; 

i=i+1;
Rnames{i} = 'IEPOXAOO = IDCHP + HO2';
k(:,i) = 1.875E13*exp(-10000./T);
Gstr{i,1} = 'IEPOXAOO';
fIEPOXAOO(i)=fIEPOXAOO(i)-1; fIDCHP(i)=fIDCHP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXAOO = OH + CO + MVK3OH4OH';
k(:,i) = 1.0E7*exp(-5000./T);
Gstr{i,1} = 'IEPOXAOO'; 
fIEPOXAOO(i)=fIEPOXAOO(i)-1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+1; fMVK3OH4OH(i)=fMVK3OH4OH(i)+1; 

i=i+1;
Rnames{i} = 'IEPOXAOO + HO2 = 0.13CO + 0.65OH + 0.65HO2 + 0.13MVK3OH4OH + 0.52GLYC + 0.52MGLY + 0.35ICPDH';
k(:,i) = 2.38E-13*exp(1300./T);
Gstr{i,1} = 'IEPOXAOO'; Gstr{i,2} = 'HO2'; 
fIEPOXAOO(i)=fIEPOXAOO(i)-1; fHO2(i)=fHO2(i)-0.35; fCO(i)=fCO(i)+0.13; fOH(i)=fOH(i)+0.65; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.13; fGLYC(i)=fGLYC(i)+0.52; fMGLY(i)=fMGLY(i)+0.52; fICPDH(i)=fICPDH(i)+0.35; 

i=i+1;
Rnames{i} = 'IEPOXAOO + NO = 0.2MVK3OH4OH + HO2 + NO2 + 0.2CO + 0.8GLYC + 0.8MGLY';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,13.098,8,1,0);
Gstr{i,1} = 'IEPOXAOO'; Gstr{i,2} = 'NO'; 
fIEPOXAOO(i)=fIEPOXAOO(i)-1; fNO(i)=fNO(i)-1; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.2; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.2; fGLYC(i)=fGLYC(i)+0.8; fMGLY(i)=fMGLY(i)+0.8; 

i=i+1;
Rnames{i} = 'IEPOXAOO + NO = IDHCN';
k(:,i) =  F0AM_isop_NIT(T,M,2.7E-12,350,13.098,8,1,0);
Gstr{i,1} = 'IEPOXAOO'; Gstr{i,2} = 'NO'; 
fIEPOXAOO(i)=fIEPOXAOO(i)-1; fNO(i)=fNO(i)-1; fIDHCN(i)=fIDHCN(i)+1;

i=i+1;
Rnames{i} = 'MVK3OH4OH + OH = 0.4MVK3OH4CO + 0.6MVK3CO4OH + HO2';
k(:,i) = 8.7E-12*exp(70./T); 
Gstr{i,1} = 'MVK3OH4OH'; Gstr{i,2} = 'OH'; 
fMVK3OH4OH(i)=fMVK3OH4OH(i)-1; fOH(i)=fOH(i)-1; fMVK3OH4CO(i)=fMVK3OH4CO(i)+0.4; fMVK3CO4OH(i)=fMVK3CO4OH(i)+0.6; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OH4CO + OH = OH + MGLY + CO2';
k(:,i) = 5.0E-12*exp(470./T); 
Gstr{i,1} = 'MVK3OH4CO'; Gstr{i,2} = 'OH'; 
fMVK3OH4CO(i)=fMVK3OH4CO(i)-1; fMGLY(i)=fMGLY(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OH4CO = 0.5GLYX + 1.5HO2 + 0.5CH3CO3 + 0.5CO + 0.5MGLY';
k(:,i) = J35 ; 
Gstr{i,1} = 'MVK3OH4CO'; 
fMVK3OH4CO(i)=fMVK3OH4CO(i)-1; fGLYX(i)=fGLYX(i)+0.5; fHO2(i)=fHO2(i)+1.5; fCH3CO3(i)=fCH3CO3(i)+0.5; fCO(i)=fCO(i)+0.5; fMGLY(i)=fMGLY(i)+0.5; 

i=i+1;
Rnames{i} = 'MVK3CO4OH + OH = 2CO + HO2 + CH3CO3';
k(:,i) = 2.0E-12*exp(70./T); 
Gstr{i,1} = 'MVK3CO4OH'; Gstr{i,2} = 'OH'; 
fMVK3CO4OH(i)=fMVK3CO4OH(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+2; fHO2(i)=fHO2(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'MVK3CO4OH = CO + HO2 + HCHO + CH3CO3';
k(:,i) = J35 ; 
Gstr{i,1} = 'MVK3CO4OH'; 
fMVK3CO4OH(i)=fMVK3CO4OH(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OH3OH + OH = 0.16MACR2OH3CO + HO2 + 0.21HAC + 0.84CO + 0.63OH + 0.63CO2 + 0.63CH3CO3';
k(:,i) = 2.4E-11*exp(70./T); 
Gstr{i,1} = 'MACR2OH3OH'; Gstr{i,2} = 'OH'; 
fMACR2OH3OH(i)=fMACR2OH3OH(i)-1; fOH(i)=fOH(i)-0.37; fMACR2OH3CO(i)=fMACR2OH3CO(i)+0.16; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+0.84; fHAC(i)=fHAC(i)+0.21; fCH3CO3(i)=fCH3CO3(i)+0.63; 

i=i+1;
Rnames{i} = 'MACR2OH3CO + OH = CO2 + OH + CO + HO2 + CH3CO3';
k(:,i) = 5.0E-12*exp(470./T);
Gstr{i,1} = 'MACR2OH3CO'; Gstr{i,2} = 'OH'; 
fMACR2OH3CO(i)=fMACR2OH3CO(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'MVK + OH = MVKOHOO';
k(:,i) = 2.6E-12*exp(610./T); 
Gstr{i,1} = 'MVK'; Gstr{i,2} = 'OH'; 
fMVK(i)=fMVK(i)-1; fOH(i)=fOH(i)-1; fMVKOHOO(i)=fMVKOHOO(i)+1; 

i=i+1;
Rnames{i} = 'MVKOHOO + HO2 = 0.360CH3CO3 + 0.360GLYC + 0.665OH + 0.305HO2 + 0.255MVK3CO4OH + 0.135MVK3OOH4OH + 0.200MVK3OH4OOH + 0.050MGLY + 0.050HCHO';
k(:,i) = 2.12E-13*exp(1300./T); 
Gstr{i,1} = 'MVKOHOO'; Gstr{i,2} = 'HO2'; 
fMVKOHOO(i)=fMVKOHOO(i)-1; fHO2(i)=fHO2(i)-0.695; fCH3CO3(i)=fCH3CO3(i)+0.36; fGLYC(i)=fGLYC(i)+0.36; fOH(i)=fOH(i)+0.665; fMVK3CO4OH(i)=fMVK3CO4OH(i)+0.255; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.135; fMVK3OH4OOH(i)=fMVK3OH4OOH(i)+0.2; fMGLY(i)=fMGLY(i)+0.05; fHCHO(i)=fHCHO(i)+0.05; 

i=i+1;
Rnames{i} = 'MVKOHOO + NO = 0.758CH3CO3 + 0.758GLYC + 0.242MGLY + 0.242HCHO + 0.242HO2 + NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,4.573,6,1,0);
Gstr{i,1} = 'MVKOHOO'; Gstr{i,2} = 'NO'; 
fMVKOHOO(i)=fMVKOHOO(i)-1; fNO(i)=fNO(i)-1; fCH3CO3(i)=fCH3CO3(i)+0.758; fGLYC(i)=fGLYC(i)+0.758; fMGLY(i)=fMGLY(i)+0.242; fHCHO(i)=fHCHO(i)+0.242; fHO2(i)=fHO2(i)+0.242; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MVKOHOO + NO = 0.438MVK3OH4N + 0.562MVK3N4OH';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,4.573,6,1,0);
Gstr{i,1} = 'MVKOHOO'; Gstr{i,2} = 'NO'; 
fMVKOHOO(i)=fMVKOHOO(i)-1; fNO(i)=fNO(i)-1; fMVK3OH4N(i)=fMVK3OH4N(i)+1; fMVK3N4OH(i)=fMVK3N4OH(i)+1; 

i=i+1;
Rnames{i} = 'MACR + OH = 0.036HPAC + 0.036CO + 0.036HO2 + 0.964MCROHOO';
k(:,i) = 4.4E-12*exp(380./T); 
Gstr{i,1} = 'MACR'; Gstr{i,2} = 'OH'; 
fMACR(i)=fMACR(i)-1; fOH(i)=fOH(i)-1; fHPAC(i)=fHPAC(i)+0.036; fCO(i)=fCO(i)+0.036; fHO2(i)=fHO2(i)+0.036; fMCROHOO(i)=fMCROHOO(i)+0.964; 

i=i+1;
Rnames{i} = 'MACR + OH = MACR1OO';
k(:,i) = 2.7E-12*exp(470./T); 
Gstr{i,1} = 'MACR'; Gstr{i,2} = 'OH'; 
fMACR(i)=fMACR(i)-1; fOH(i)=fOH(i)-1; fMACR1OO(i)=fMACR1OO(i)+1; 

i=i+1;
Rnames{i} = 'MCROHOO + HO2 = 0.41MACR2OOH3OH + 0.507HAC + 0.507CO + 0.507HO2 + 0.59OH + 0.59O2 + 0.083MGLY + 0.083HCHO';
k(:,i) = 2.12E-13*exp(1300./T) ;
Gstr{i,1} = 'MCROHOO'; Gstr{i,2} = 'HO2'; 
fMCROHOO(i)=fMCROHOO(i)-1; fHO2(i)=fHO2(i)-0.493; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.41; fHAC(i)=fHAC(i)+0.507; fCO(i)=fCO(i)+0.507; fOH(i)=fOH(i)+0.59; fMGLY(i)=fMGLY(i)+0.083; fHCHO(i)=fHCHO(i)+0.083; 

i=i+1;
Rnames{i} = 'MACR1OO + HO2 = 0.37MACR1OOH + 0.5HCHO + 0.325CO + 0.325CH3OO + 0.175CH3CO3 + 0.5CO2 + 0.5OH + 0.13O3 + 0.13MACR1OH';
k(:,i) = 3.14E-12*exp(580./T) ; 
Gstr{i,1} = 'MACR1OO'; Gstr{i,2} = 'HO2'; 
fMACR1OO(i)=fMACR1OO(i)-1; fHO2(i)=fHO2(i)-1; fMACR1OOH(i)=fMACR1OOH(i)+0.37; fHCHO(i)=fHCHO(i)+0.5; fCO(i)=fCO(i)+0.325; fCH3OO(i)=fCH3OO(i)+0.325; fCH3CO3(i)=fCH3CO3(i)+0.175; fOH(i)=fOH(i)+0.5; fO3(i)=fO3(i)+0.13; fMACR1OH(i)=fMACR1OH(i)+0.13; 

i=i+1;
Rnames{i} = 'MCROHOO = HAC + CO + OH';
k(:,i) = 2.9E7*exp(-5297./T); 
Gstr{i,1} = 'MCROHOO'; 
fMCROHOO(i)=fMCROHOO(i)-1; fHAC(i)=fHAC(i)+1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'MCROHOO + NO = 0.86HAC + 0.86CO + 0.86HO2 + NO2 + 0.14MGLY + 0.14HCHO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,2.985,6,1,0);
Gstr{i,1} = 'MCROHOO'; Gstr{i,2} = 'NO'; 
fMCROHOO(i)=fMCROHOO(i)-1; fNO(i)=fNO(i)-1; fHAC(i)=fHAC(i)+0.86; fCO(i)=fCO(i)+0.86; fHO2(i)=fHO2(i)+0.86; fNO2(i)=fNO2(i)+1; fMGLY(i)=fMGLY(i)+0.14; fHCHO(i)=fHCHO(i)+0.14; 

i=i+1;
Rnames{i} = 'MCROHOO + NO = MACR2N3OH';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,2.985,6,1,0);
Gstr{i,1} = 'MCROHOO'; Gstr{i,2} = 'NO'; 
fMCROHOO(i)=fMCROHOO(i)-1; fNO(i)=fNO(i)-1; fMACR2N3OH(i)=fMACR2N3OH(i)+1;

i=i+1;
Rnames{i} = 'MACR1OO + NO = 0.35CH3CO3 + 0.65CH3OO + 0.65CO + HCHO + CO2 + NO2';
k(:,i) = 8.7E-12*exp(290./T); 
Gstr{i,1} = 'MACR1OO'; Gstr{i,2} = 'NO'; 
fMACR1OO(i)=fMACR1OO(i)-1; fNO(i)=fNO(i)-1; fCH3CO3(i)=fCH3CO3(i)+0.35; fCH3OO(i)=fCH3OO(i)+0.65; fCO(i)=fCO(i)+0.65; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MACR1OO + NO2 = MPAN';
k(:,i) = F0AM_isop_TROE(2.591E-28,0,-6.87,1.125E-11,0,-1.105,0.3,T,M);
Gstr{i,1} = 'MACR1OO'; Gstr{i,2} = 'NO2'; 
fMACR1OO(i)=fMACR1OO(i)-1; fNO2(i)=fNO2(i)-1; fMPAN(i)=fMPAN(i)+1; 

i=i+1;
Rnames{i} = 'MPAN = MACR1OO + NO2';
k(:,i) = 1.58E16*exp(-13500./T);
Gstr{i,1} = 'MPAN'; 
fMPAN(i)=fMPAN(i)-1; fMACR1OO(i)=fMACR1OO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MPAN + OH = 0.75HMML + NO3 + 0.25HAC + 0.25CO';
k(:,i) = 2.9E-11; 
Gstr{i,1} = 'MPAN'; Gstr{i,2} = 'OH'; 
fMPAN(i)=fMPAN(i)-1; fOH(i)=fOH(i)-1; fHMML(i)=fHMML(i)+0.75; fNO3(i)=fNO3(i)+1; fHAC(i)=fHAC(i)+0.25; fCO(i)=fCO(i)+0.25; 

i=i+1;
Rnames{i} = 'ISOP1CO4OOH = 0.888CO + 1.662OH + 0.112HO2 + 0.112ISOP1CO4CO + 0.112MVK3OOH4CO + 0.552MVKENOL + 0.224C4HVP1';
k(:,i) = J20 ; % SUN*2.32E-4;  
Gstr{i,1} = 'ISOP1CO4OOH'; 
fISOP1CO4OOH(i)=fISOP1CO4OOH(i)-1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+0.888; fHO2(i)=fHO2(i)+1.662; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.112; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.112; fMVKENOL(i)=fMVKENOL(i)+0.552; fC4HVP1(i)=fC4HVP1(i)+0.224; 

i=i+1;
Rnames{i} = 'ISOP1OOH4CO = 0.818CO + 1.637OH + 0.182HO2 + 0.182ISOP1CO4CO + 0.182MVK3OOH4CO + 0.455MCRENOL + 0.182C4HVP2';
k(:,i) = J20 ; % SUN*2.2E-4;  
Gstr{i,1} = 'ISOP1OOH4CO';
fISOP1OOH4CO(i)=fISOP1OOH4CO(i)-1; fCO(i)=fCO(i)+0.818; fOH(i)=fOH(i)+1.637; fHO2(i)=fHO2(i)+0.182; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.182; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.182; fMCRENOL(i)=fMCRENOL(i)+0.455; fC4HVP2(i)=fC4HVP2(i)+0.182; 

i=i+1;
Rnames{i} = 'ISOP1CO4OOH + OH = 0.035MVK + 0.315HPALD1OO + 0.15ISOP1CO4CO + 0.35MVK3OH4OOH + 0.075HO2 + 0.075HCHO + 0.075MGLY + 0.075ICHE + 1.075OH + 0.46CO';
k(:,i) = 1.17E-11*exp(450./T) ; 
Gstr{i,1} = 'ISOP1CO4OOH'; Gstr{i,2} = 'OH'; 
fISOP1CO4OOH(i)=fISOP1CO4OOH(i)-1; fOH(i)=fOH(i)+0.075; fMVK(i)=fMVK(i)+0.035; fHPALD1OO(i)=fHPALD1OO(i)+0.315; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.15; fMVK3OH4OOH(i)=fMVK3OH4OOH(i)+0.35; fHO2(i)=fHO2(i)+0.075; fHCHO(i)=fHCHO(i)+0.075; fMGLY(i)=fMGLY(i)+0.075; fICHE(i)=fICHE(i)+0.075; fCO(i)=fCO(i)+0.46; 

i=i+1;
Rnames{i} = 'ISOP1OOH4CO + OH = 0.035MACR + 0.315HPALD2OO + 0.15ISOP1CO4CO + 0.15MACR2OH3OOH + 0.175HO2 + 0.175HCHO + 0.175MGLY + 0.175ICHE + 1.175OH + 0.36CO';
k(:,i) = 1.17E-11*exp(450./T) ; 
Gstr{i,1} = 'ISOP1OOH4CO'; Gstr{i,2} = 'OH'; 
fISOP1OOH4CO(i)=fISOP1OOH4CO(i)-1; fOH(i)=fOH(i)+0.175; fMACR(i)=fMACR(i)+0.035; fHPALD2OO(i)=fHPALD2OO(i)+0.315; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.15; fMACR2OH3OOH(i)=fMACR2OH3OOH(i)+0.15; fHO2(i)=fHO2(i)+0.175; fHCHO(i)=fHCHO(i)+0.175; fMGLY(i)=fMGLY(i)+0.175; fICHE(i)=fICHE(i)+0.175; fCO(i)=fCO(i)+0.36; 

i=i+1;
Rnames{i} = 'C4HVP1 + NO = NO2 + MVKOHOO';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'C4HVP1'; Gstr{i,2} = 'NO'; 
fC4HVP1(i)=fC4HVP1(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMVKOHOO(i)=fMVKOHOO(i)+1; 

i=i+1;
Rnames{i} = 'C4HVP1 + HO2 = OH + MVKOHOO';
k(:,i) = 1.93E-13*exp(1300./T); 
Gstr{i,1} = 'C4HVP1'; Gstr{i,2} = 'HO2'; 
fC4HVP1(i)=fC4HVP1(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+1; fMVKOHOO(i)=fMVKOHOO(i)+1; 

i=i+1;
Rnames{i} = 'C4HVP1 + NO2 = MVK3N4OH';
k(:,i) = 9.0E-12; 
Gstr{i,1} = 'C4HVP1'; Gstr{i,2} = 'NO2'; 
fC4HVP1(i)=fC4HVP1(i)-1; fNO2(i)=fNO2(i)-1; fMVK3N4OH(i)=fMVK3N4OH(i)+1; 

i=i+1;
Rnames{i} = 'C4HVP2 + NO = NO2 + MCROHOO';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'C4HVP2'; Gstr{i,2} = 'NO'; 
fC4HVP2(i)=fC4HVP2(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fMCROHOO(i)=fMCROHOO(i)+1; 

i=i+1;
Rnames{i} = 'C4HVP2 + HO2 = OH + MCROHOO';
k(:,i) = 1.93E-13*exp(1300./T); 
Gstr{i,1} = 'C4HVP2'; Gstr{i,2} = 'HO2'; 
fC4HVP2(i)=fC4HVP2(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+1; fMCROHOO(i)=fMCROHOO(i)+1; 

i=i+1;
Rnames{i} = 'C4HVP2 + NO2 = MACR2N3OH';
k(:,i) = 9.0E-12; 
Gstr{i,1} = 'C4HVP2'; Gstr{i,2} = 'NO2'; 
fC4HVP2(i)=fC4HVP2(i)-1; fNO2(i)=fNO2(i)-1; fMACR2N3OH(i)=fMACR2N3OH(i)+1; 

i=i+1;
Rnames{i} = 'MCRENOL = CO + PYRAC + 2OH';
k(:,i) = J35 ; % 2.5E-4*SUN; 
Gstr{i,1} = 'MCRENOL'; 
fMCRENOL(i)=fMCRENOL(i)-1; fCO(i)=fCO(i)+1; fPYRAC(i)=fPYRAC(i)+1; fOH(i)=fOH(i)+2; 

i=i+1;
Rnames{i} = 'MVKENOL = 0.5MGLY + 0.5HO2 + 0.5CO + 0.5CH3CO3 + 0.5GLYX + OH';
k(:,i) = J35 ; % 2.5E-4*SUN; 
Gstr{i,1} = 'MVKENOL'; 
fMVKENOL(i)=fMVKENOL(i)-1; fMGLY(i)=fMGLY(i)+0.5; fHO2(i)=fHO2(i)+0.5; fCO(i)=fCO(i)+0.5; fCH3CO3(i)=fCH3CO3(i)+0.5; fGLYX(i)=fGLYX(i)+0.5; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'MCRENOL + OH = CO + 0.87PYRAC + 0.87HO2 + 0.13OH + 0.13CO2 + 0.13CH3CO3';
k(:,i) = 3.83E-12*exp(983./T); 
Gstr{i,1} = 'MCRENOL'; Gstr{i,2} = 'OH'; 
fMCRENOL(i)=fMCRENOL(i)-1; fOH(i)=fOH(i)-0.87; fCO(i)=fCO(i)+1; fPYRAC(i)=fPYRAC(i)+0.87; fHO2(i)=fHO2(i)+0.87; fCH3CO3(i)=fCH3CO3(i)+0.13; 

i=i+1;
Rnames{i} = 'MVKENOL + OH = 0.25HO2 + 0.25MVK3OH4CO + 0.75OH + 0.75MGLY + 0.75HCOOH';
k(:,i) = 3.35E-12*exp(983./T); 
Gstr{i,1} = 'MVKENOL'; Gstr{i,2} = 'OH'; 
fMVKENOL(i)=fMVKENOL(i)-1; fOH(i)=fOH(i)-0.25; fHO2(i)=fHO2(i)+0.25; fMVK3OH4CO(i)=fMVK3OH4CO(i)+0.25; fMGLY(i)=fMGLY(i)+0.75; fHCOOH(i)=fHCOOH(i)+0.75; 

i=i+1;
Rnames{i} = 'MVK3OOH4CO + OH = OH + CO + MGLY';
k(:,i) = 5.0E-12*exp(470./T); 
Gstr{i,1} = 'MVK3OOH4CO'; Gstr{i,2} = 'OH'; 
fMVK3OOH4CO(i)=fMVK3OOH4CO(i)-1; fCO(i)=fCO(i)+1; fMGLY(i)=fMGLY(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OOH4CO = OH + 0.571CO + 0.571MGLY + 0.571HO2 + 0.429GLYX + 0.429CH3CO3';
k(:,i) = J35 ; 
Gstr{i,1} = 'MVK3OOH4CO'; 
fMVK3OOH4CO(i)=fMVK3OOH4CO(i)-1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+0.571; fMGLY(i)=fMGLY(i)+0.571; fHO2(i)=fHO2(i)+0.571; fGLYX(i)=fGLYX(i)+0.429; fCH3CO3(i)=fCH3CO3(i)+0.429; 

i=i+1;
Rnames{i} = 'HPALD1OO + NO = NO2 + OH + CO2 + MVK';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'HPALD1OO'; Gstr{i,2} = 'NO'; 
fHPALD1OO(i)=fHPALD1OO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1; fMVK(i)=fMVK(i)+1; 

i=i+1;
Rnames{i} = 'HPALD1OO + HO2 = OH + OH + CO2 + MVK';
k(:,i) = 2.38E-12*exp(1300./T); 
Gstr{i,1} = 'HPALD1OO'; Gstr{i,2} = 'HO2'; 
fHPALD1OO(i)=fHPALD1OO(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+2; fMVK(i)=fMVK(i)+1; 

i=i+1;
Rnames{i} = 'HPALD2OO + NO = NO2 + OH + CO2 + MACR';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'HPALD2OO'; Gstr{i,2} = 'NO'; 
fHPALD2OO(i)=fHPALD2OO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1; fMACR(i)=fMACR(i)+1; 

i=i+1;
Rnames{i} = 'HPALD2OO + HO2 = OH + OH + CO2 + MACR';
k(:,i) = 2.38E-12*exp(1300./T); 
Gstr{i,1} = 'HPALD2OO'; Gstr{i,2} = 'HO2'; 
fHPALD2OO(i)=fHPALD2OO(i)-1; fHO2(i)=fHO2(i)-1; fOH(i)=fOH(i)+2; fMACR(i)=fMACR(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH2N + OH = ISOPNOO1';
k(:,i) = 7.14E-12*exp(390./T); 
Gstr{i,1} = 'ISOP1OH2N'; Gstr{i,2} = 'OH'; 
fISOP1OH2N(i)=fISOP1OH2N(i)-1; fOH(i)=fOH(i)-1; fISOPNOO1(i)=fISOPNOO1(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH2N + OH = 0.67IEPOXt + 0.33IEPOXc + NO2';
k(:,i) = F0AM_isop_EPO(T,M,6.30E-12,390,1.62E-19);
Gstr{i,1} = 'ISOP1OH2N'; Gstr{i,2} = 'OH'; 
fISOP1OH2N(i)=fISOP1OH2N(i)-1; fOH(i)=fOH(i)-1; fIEPOXt(i)=fIEPOXt(i)+0.67; fIEPOXc(i)=fIEPOXc(i)+0.33; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3N4OH + OH = ISOPNOO2';
k(:,i) = 1.02E-11*exp(390./T); 
Gstr{i,1} = 'ISOP3N4OH'; Gstr{i,2} = 'OH'; 
fISOP3N4OH(i)=fISOP3N4OH(i)-1; fOH(i)=fOH(i)-1; fISOPNOO2(i)=fISOPNOO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3N4OH + OH = 0.67IEPOXt + 0.33IEPOXc + NO2';
k(:,i) = F0AM_isop_EPO(T,M,1.05E-11,390,2.49E-19);
Gstr{i,1} = 'ISOP3N4OH'; Gstr{i,2} = 'OH'; 
fISOP3N4OH(i)=fISOP3N4OH(i)-1; fOH(i)=fOH(i)-1; fIEPOXt(i)=fIEPOXt(i)+0.67; fIEPOXc(i)=fIEPOXc(i)+0.33; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOPNOO1 = ICHNP + HO2';
k(:,i) = 1.875E13*exp(-10000./T); 
Gstr{i,1} = 'ISOPNOO1'; 
fISOPNOO1(i)=fISOPNOO1(i)-1; fICHNP(i)=fICHNP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOPNOO1 + HO2 = 0.482IDHPN + 0.059MACR2N3OH + 0.059HCHO + 0.459GLYC + 0.459HAC + 0.059HO2 + 0.459NO2 + 0.518OH';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'ISOPNOO1'; Gstr{i,2} = 'HO2'; 
fISOPNOO1(i)=fISOPNOO1(i)-1; fHO2(i)=fHO2(i)-0.941; fIDHPN(i)=fIDHPN(i)+0.482; fMACR2N3OH(i)=fMACR2N3OH(i)+0.059; fHCHO(i)=fHCHO(i)+0.059; fGLYC(i)=fGLYC(i)+0.459; fHAC(i)=fHAC(i)+0.459; fNO2(i)=fNO2(i)+0.459; fOH(i)=fOH(i)+0.518; 

i=i+1;
Rnames{i} = 'ISOPNOO1 + NO = 0.272MACR2N3OH + 0.272HCHO + 0.728GLYC + 0.728HAC + 0.272HO2 + 1.728NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,6.32,11,1,0);
Gstr{i,1} = 'ISOPNOO1'; Gstr{i,2} = 'NO'; 
fISOPNOO1(i)=fISOPNOO1(i)-1; fNO(i)=fNO(i)-1; fMACR2N3OH(i)=fMACR2N3OH(i)+0.272; fHCHO(i)=fHCHO(i)+0.272; fGLYC(i)=fGLYC(i)+0.782; fHAC(i)=fHAC(i)+0.782; fHO2(i)=fHO2(i)+0.272; fNO2(i)=fNO2(i)+1.782; 

i=i+1;
Rnames{i} = 'ISOPNOO1 + NO = IDHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,6.32,11,1,0.);
Gstr{i,1} = 'ISOPNOO1'; Gstr{i,2} = 'NO'; 
fISOPNOO1(i)=fISOPNOO1(i)-1; fNO(i)=fNO(i)-1; fIDHDN(i)=fIDHDN(i)+1;

i=i+1;
Rnames{i} = 'ISOPNOO2 = ICHNP + HO2';
k(:,i) = 1.875E13*exp(-10000./T); 
Gstr{i,1} = 'ISOPNOO2'; 
fISOPNOO2(i)=fISOPNOO2(i)-1; fICHNP(i)=fICHNP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOPNOO2 + HO2 = 0.401IDHPN + 0.599MVK3N4OH + 0.599HCHO + 0.599HO2 + 0.599OH';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'ISOPNOO2'; Gstr{i,2} = 'HO2'; 
fISOPNOO2(i)=fISOPNOO2(i)-1; fHO2(i)=fHO2(i)-0.401; fIDHPN(i)=fIDHPN(i)+0.401; fMVK3N4OH(i)=fMVK3N4OH(i)+0.599; fHCHO(i)=fHCHO(i)+0.599; fOH(i)=fOH(i)+0.599; 

i=i+1;
Rnames{i} = 'ISOPNOO2 + NO = MVK3N4OH + HCHO + HO2 + NO2';
k(:,i) =F0AM_isop_ALK(T,M,2.7E-12,350,7.941,11,1,0.);
Gstr{i,1} = 'ISOPNOO2'; Gstr{i,2} = 'NO'; 
fISOPNOO2(i)=fISOPNOO2(i)-1; fNO(i)=fNO(i)-1; fMVK3N4OH(i)=fMVK3N4OH(i)+1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOPNOO2 + NO = IDHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,7.941,11,1,0.);
Gstr{i,1} = 'ISOPNOO2'; Gstr{i,2} = 'NO'; 
fISOPNOO2(i)=fISOPNOO2(i)-1; fNO(i)=fNO(i)-1; fIDHDN(i)=fIDHDN(i)+1; 

i=i+1;
Rnames{i} = 'ISOP + O3 = 0.416MACR + 0.177MVK + 0.58SCI + 0.28OH + 0.16HO2 + 0.827HCHO + 0.407CO2 + 0.407CO + 0.407CH3OO + 0.013H2O2';
k(:,i) =  1.3E-17 ;
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'O3'; 
fISOP(i)=fISOP(i)-1; fO3(i)=fO3(i)-1; fMACR(i)=fMACR(i)+0.416; fMVK(i)=fMVK(i)+0.177; fSCI(i)=fSCI(i)+0.58; fOH(i)=fOH(i)+0.28; fHO2(i)=fHO2(i)+0.16; fHCHO(i)=fHCHO(i)+0.827; fCO(i)=fCO(i)+0.407; fCH3OO(i)=fCH3OO(i)+0.407; fH2O2(i)=fH2O2(i)+0.013; 

i=i+1;
Rnames{i} = 'SCI + H2O = 0.730HMHP + 0.210HCOOH + 0.060HCHO + 0.060H2O2';
k(:,i) = 1.7E-15*H2O ;
Gstr{i,1} = 'SCI'; 
fSCI(i)=fSCI(i)-1; fHMHP(i)=fHMHP(i)+0.73; fHCOOH(i)=fHCOOH(i)+0.21; fHCHO(i)=fHCHO(i)+0.06; fH2O2(i)=fH2O2(i)+0.06; 

i=i+1;
Rnames{i} = 'SCI + H2O + H2O = 0.400HMHP + 0.540HCOOH + 0.060HCHO + 0.060H2O2';
k(:,i) = 2.88E-35*exp(1391./T).*H2O.*H2O ;
Gstr{i,1} = 'SCI'; 
fSCI(i)=fSCI(i)-1; fHMHP(i)=fHMHP(i)+0.4; fHCOOH(i)=fHCOOH(i)+0.54; fHCHO(i)=fHCHO(i)+0.06; fH2O2(i)=fH2O2(i)+0.06; 

i=i+1;
Rnames{i} = 'SCI + O3 = HCHO';
k(:,i) = 2.0E-12*0.7; 
Gstr{i,1} = 'SCI'; Gstr{i,2} = 'O3'; 
fSCI(i)=fSCI(i)-1; fO3(i)=fO3(i)-1; fHCHO(i)=fHCHO(i)+1;

i=i+1;
Rnames{i} = 'HMHP + OH = 0.5HCHO + 0.5HO2 + 0.5HCOOH + 0.5OH';
k(:,i) = 1.3E-12*exp(500./T); 
Gstr{i,1} = 'HMHP'; Gstr{i,2} = 'OH'; 
fHMHP(i)=fHMHP(i)-1; fOH(i)=fOH(i)-0.5; fHCHO(i)=fHCHO(i)+0.5; fHO2(i)=fHO2(i)+0.5; fHCOOH(i)=fHCOOH(i)+0.5; 

i=i+1;
Rnames{i} = 'ISOP + NO3 = 0.465INO2B + 0.535INO2D';
k(:,i) = 2.95E-12*exp(-450./T) ; 
Gstr{i,1} = 'ISOP'; Gstr{i,2} = 'NO3'; 
fISOP(i)=fISOP(i)-1; fNO3(i)=fNO3(i)-1; fINO2B(i)=fINO2B(i)+0.465; fINO2D(i)=fINO2D(i)+0.535; 

i=i+1;
Rnames{i} = 'INO2B + HO2 = 0.473INPB + 0.048MACR + 0.479MVK + 0.527OH + 0.527HCHO + 0.527NO2';
k(:,i) = 2.47E-13*exp(1300./T); 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'HO2'; 
fINO2B(i)=fINO2B(i)-1; fHO2(i)=fHO2(i)-1; fINPB(i)=fINPB(i)+0.473; fMACR(i)=fMACR(i)+0.048; fMVK(i)=fMVK(i)+0.479; fOH(i)=fOH(i)+0.527; fHCHO(i)=fHCHO(i)+0.527; fNO2(i)=fNO2(i)+0.527; 

i=i+1;
Rnames{i} = 'INO2D + HO2 = INPD';
k(:,i) = 2.47E-13*exp(1300./T); 
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'HO2'; 
fINO2D(i)=fINO2D(i)-1; fHO2(i)=fHO2(i)-1; fINPD(i)=fINPD(i)+1;

i=i+1;
Rnames{i} = 'INO2B + INO2D = 0.336ISOP1N4CO + 0.399IHNB + 0.544MVK + 0.563HCHO + 0.563NO2 + 0.474INO + 0.019MACR + 0.152ISOP1CO4N + 0.032ISOP1N4OH + 0.006ISOP1OH4N + 0.038ISOP3CO4N + 0.089HO2';
k(:,i) = 2.56E-12; 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'INO2D'; 
fINO2B(i)=fINO2B(i)-1; fINO2D(i)=fINO2D(i)-1; fISOP1N4CO(i)=fISOP1N4CO(i)+0.336; fIHNB(i)=fIHNB(i)+0.399; fMVK(i)=fMVK(i)+0.544; fHCHO(i)=fHCHO(i)+0.563; fNO2(i)=fNO2(i)+0.563; fINO(i)=fINO(i)+0.474; fMACR(i)=fMACR(i)+0.019; fISOP1CO4N(i)=fISOP1CO4N(i)+0.152; fISOP1N4OH(i)=fISOP1N4OH(i)+0.032; fISOP1OH4N(i)=fISOP1OH4N(i)+0.006; fISOP3CO4N(i)=fISOP3CO4N(i)+0.038; fHO2(i)=fHO2(i)+0.089; 

i=i+1;
Rnames{i} = 'INO2D + INO2D  = 0.064HO2 + 0.191ISOP1CO4N + 0.340INO + 0.67ISOP1N4CO + 0.671ISOP1N4OH + 0.127ISOP1OH4N';
k(:,i) = 3.71E-12; 
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'INO2D'; 
fINO2D(i)=fINO2D(i)-2; fHO2(i)=fHO2(i)+0.064; fISOP1CO4N(i)=fISOP1CO4N(i)+0.191; fINO(i)=fINO(i)+0.34; fISOP1N4CO(i)=fISOP1N4CO(i)+0.671; fISOP1N4OH(i)=fISOP1N4OH(i)+0.671; fISOP1OH4N(i)=fISOP1OH4N(i)+0.127; 

i=i+1;
Rnames{i} = 'INO2B + INO2B = 1.737MVK + 0.123MACR + 1.860HCHO + 1.860NO2 + 0.070IHNB + 0.070ISOP3CO4N';
k(:,i) = 1.61E-12 ; 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'INO2B'; 
fINO2B(i)=fINO2B(i)-2; fMVK(i)=fMVK(i)+1.737; fMACR(i)=fMACR(i)+0.123; fHCHO(i)=fHCHO(i)+1.86; fNO2(i)=fNO2(i)+1.86; fIHNB(i)=fIHNB(i)+0.07; fISOP3CO4N(i)=fISOP3CO4N(i)+0.07; 

i=i+1;
Rnames{i} = 'INO2B + CH3OO = 0.355IHNB + 0.583MVK + 0.028MACR + 0.034ISOP3CO4N + 0.611HO2 + 1.577HCHO + 0.611NO2 + 0.034CH3OH';
k(:,i) = 2.80E-13 ; 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'CH3OO'; 
fINO2B(i)=fINO2B(i)-1; fCH3OO(i)=fCH3OO(i)-1; fIHNB(i)=fIHNB(i)+0.355; fMVK(i)=fMVK(i)+0.583; fMACR(i)=fMACR(i)+0.028; fISOP3CO4N(i)=fISOP3CO4N(i)+0.034; fHO2(i)=fHO2(i)+0.611; fHCHO(i)=fHCHO(i)+1.577; fNO2(i)=fNO2(i)+0.611; fCH3OH(i)=fCH3OH(i)+0.034; 

i=i+1;
Rnames{i} = 'INO2D + CH3OO = 0.102ISOP1CO4N + 0.298ISOP1N4OH + 0.057ISOP1OH4N + 0.244INO + 0.299ISOP1N4CO + 0.355CH3OH + 0.336HO2 + 0.645HCHO';
k(:,i) = 1.18E-12 ; 
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'CH3OO'; 
fINO2D(i)=fINO2D(i)-1; fCH3OO(i)=fCH3OO(i)-1; fISOP1CO4N(i)=fISOP1CO4N(i)+0.102; fISOP1N4OH(i)=fISOP1N4OH(i)+0.298; fISOP1OH4N(i)=fISOP1OH4N(i)+0.057; fINO(i)=fINO(i)+0.244; fISOP1N4CO(i)=fISOP1N4CO(i)+0.299; fCH3OH(i)=fCH3OH(i)+0.355; fHO2(i)=fHO2(i)+0.336; fHCHO(i)=fHCHO(i)+0.645; 

i=i+1;
Rnames{i} = 'INO2B + CH3CO3 = HCHO + NO2 + CH3OO + CO2 + 0.903MVK + 0.097MACR';
k(:,i) = 1.92E-12 ; 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'CH3CO3'; 
fINO2B(i)=fINO2B(i)-1; fCH3CO3(i)=fCH3CO3(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; fCH3OO(i)=fCH3OO(i)+1; fMVK(i)=fMVK(i)+0.903; fMACR(i)=fMACR(i)+0.097; 

i=i+1;
Rnames{i} = 'INO2D + CH3CO3 = CH3OO + CO2 + 0.841INO + 0.159HO2 + 0.159ISOP1CO4N';
k(:,i) = 7.71E-12 ; 
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'CH3CO3'; 
fINO2D(i)=fINO2D(i)-1; fCH3CO3(i)=fCH3CO3(i)-1; fCH3OO(i)=fCH3OO(i)+1; fINO(i)=fINO(i)+0.841; fHO2(i)=fHO2(i)+0.159; fISOP1CO4N(i)=fISOP1CO4N(i)+0.159; 

i=i+1;
Rnames{i} = 'INO2B + NO3 = HCHO + 2NO2 + 0.903MVK + 0.097MACR';
k(:,i) = 2.3E-12 ; 
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'NO3'; 
fINO2B(i)=fINO2B(i)-1; fNO3(i)=fNO3(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+2; fMVK(i)=fMVK(i)+0.903; fMACR(i)=fMACR(i)+0.097; 

i=i+1;
Rnames{i} = 'INO2D + NO3 = NO2 + 0.841INO + 0.159HO2 + 0.159ISOP1CO4N';
k(:,i) = 2.3E-12 ; 
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'NO3'; 
fINO2D(i)=fINO2D(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fINO(i)=fINO(i)+0.841; fHO2(i)=fHO2(i)+0.159; fISOP1CO4N(i)=fISOP1CO4N(i)+0.159; 

i=i+1;
Rnames{i} = 'INO2B + NO = 2NO2 + HCHO + 0.096MACR + 0.904MVK';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,12.915,9,1,0.);
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'NO'; 
fINO2B(i)=fINO2B(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+2; fHCHO(i)=fHCHO(i)+1; fMACR(i)=fMACR(i)+0.096; fMVK(i)=fMVK(i)+0.904; 

i=i+1;
Rnames{i} = 'INO2B + NO = IDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,12.915,9,1,0.);
Gstr{i,1} = 'INO2B'; Gstr{i,2} = 'NO'; 
fINO2B(i)=fINO2B(i)-1; fNO(i)=fNO(i)-1; fIDN(i)=fIDN(i)+1; 

i=i+1;
Rnames{i} = 'INO2D + NO = NO2 + 0.159HO2 + 0.159ISOP1CO4N + 0.841INO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.412,9,1,0.);
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'NO'; 
fINO2D(i)=fINO2D(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.159; fISOP1CO4N(i)=fISOP1CO4N(i)+0.159; fINO(i)=fINO(i)+0.841; 

i=i+1;
Rnames{i} = 'INO2D + NO = IDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.412,9,1,0);
Gstr{i,1} = 'INO2D'; Gstr{i,2} = 'NO'; 
fINO2D(i)=fINO2D(i)-1; fNO(i)=fNO(i)-1; fIDN(i)=fIDN(i)+1; 

i=i+1;
Rnames{i} = 'INO + O2 = ISOP1N4CO + HO2';
k(:,i) = 2.5E-14*exp(-300./T)*0.79.*M ; 
Gstr{i,1} = 'INO'; 
fINO(i)=fINO(i)-1; fISOP1N4CO(i)=fISOP1N4CO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'INO = IDHNBOO';
k(:,i) = 1.0E20*exp(-10000./T) ; 
Gstr{i,1} = 'INO';
fINO(i)=fINO(i)-1; fIDHNBOO(i)=fIDHNBOO(i)+1; 

i=i+1;
Rnames{i} = 'IHNB + OH = IDHNBOO';
k(:,i) = 8.72E-12*exp(390./T); 
Gstr{i,1} = 'IHNB'; Gstr{i,2} = 'OH'; 
fIHNB(i)=fIHNB(i)-1; fOH(i)=fOH(i)-1; fIDHNBOO(i)=fIDHNBOO(i)+1;

i=i+1;
Rnames{i} = 'ISOP1N4OH + OH = IEPOXD + NO2';
k(:,i) = F0AM_isop_EPO(T,M,1.55E-11,390,2.715E-19);
Gstr{i,1} = 'ISOP1N4OH'; Gstr{i,2} = 'OH'; 
fISOP1N4OH(i)=fISOP1N4OH(i)-1; fOH(i)=fOH(i)-1; fIEPOXD(i)=fIEPOXD(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1N4OH + OH = IDHNDOO1';
k(:,i) = 2.04E-11*exp(390./T); 
Gstr{i,1} = 'ISOP1N4OH'; Gstr{i,2} = 'OH'; 
fISOP1N4OH(i)=fISOP1N4OH(i)-1; fOH(i)=fOH(i)-1; fIDHNDOO1(i)=fIDHNDOO1(i)+1;

i=i+1;
Rnames{i} = 'ISOP1OH4N + OH = IEPOXD + NO2';
k(:,i) = F0AM_isop_EPO(T,M,9.52E-12,390,2.715E-19);
Gstr{i,1} = 'ISOP1OH4N'; Gstr{i,2} = 'OH'; 
fISOP1OH4N(i)=fISOP1OH4N(i)-1; fOH(i)=fOH(i)-1; fIEPOXD(i)=fIEPOXD(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH4N + OH = IDHNDOO2';
k(:,i) = 2.95E-11*exp(390./T); 
Gstr{i,1} = 'ISOP1OH4N'; Gstr{i,2} = 'OH'; 
fISOP1OH4N(i)=fISOP1OH4N(i)-1; fOH(i)=fOH(i)-1; fIDHNDOO2(i)=fIDHNDOO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1OH4N + OH = 0.6OH + 0.6CO + 0.6MVK3OOH4N + 0.4HO2 + 0.4ISOP1CO4N';
k(:,i) = 7.5E-12*exp(20./T);
Gstr{i,1} = 'ISOP1OH4N'; Gstr{i,2} = 'OH'; 
fISOP1OH4N(i)=fISOP1OH4N(i)-1; fOH(i)=fOH(i)-0.4; fCO(i)=fCO(i)+0.6; fHO2(i)=fHO2(i)+0.4; fMVK3OOH4N(i)=fMVK3OOH4N(i)+0.6; fISOP1CO4N(i)=fISOP1CO4N(i)+0.4;

i=i+1;
Rnames{i} = 'ISOP1N4OH + OH = 0.6OH + 0.6CO + 0.6MACR2OOH3N + 0.4HO2 + 0.4ISOP1N4CO';
k(:,i) = 7.5E-12*exp(20./T);
Gstr{i,1} = 'ISOP1N4OH'; Gstr{i,2} = 'OH'; 
fISOP1N4OH(i)=fISOP1N4OH(i)-1; fOH(i)=fOH(i)-0.4; fCO(i)=fCO(i)+0.6; fHO2(i)=fHO2(i)+0.4; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.6; fISOP1N4CO(i)=fISOP1N4CO(i)+0.4; 

i=i+1;
Rnames{i} = 'IDHNDOO1 + HO2 = 0.418IDHPN + 0.551PROPNN + 0.551GLYC + 0.031MACR2OH3N + 0.031HCHO + 0.582HO2 + 0.582OH';
k(:,i) = 2.6E-13*exp(1300./T);
Gstr{i,1} = 'IDHNDOO1'; Gstr{i,2} = 'HO2'; 
fIDHNDOO1(i)=fIDHNDOO1(i)-1; fHO2(i)=fHO2(i)-0.418; fIDHPN(i)=fIDHPN(i)+0.418; fPROPNN(i)=fPROPNN(i)+0.551; fGLYC(i)=fGLYC(i)+0.551; fMACR2OH3N(i)=fMACR2OH3N(i)+0.031; fHCHO(i)=fHCHO(i)+0.031; fOH(i)=fOH(i)+0.582; 

i=i+1;
Rnames{i} = 'IDHNDOO1 + NO = 0.935PROPNN + 0.935GLYC + 0.065MACR2OH3N + 0.065HCHO + HO2 + NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,4.712,11,1,0);
Gstr{i,1} = 'IDHNDOO1'; Gstr{i,2} = 'NO'; 
fIDHNDOO1(i)=fIDHNDOO1(i)-1; fNO(i)=fNO(i)-1; fPROPNN(i)=fPROPNN(i)+0.935; fGLYC(i)=fGLYC(i)+0.935; fMACR2OH3N(i)=fMACR2OH3N(i)+0.065; fHCHO(i)=fHCHO(i)+0.065; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'IDHNDOO1 + NO = IDHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,4.712,11,1,0);
Gstr{i,1} = 'IDHNDOO1'; Gstr{i,2} = 'NO'; 
fIDHNDOO1(i)=fIDHNDOO1(i)-1; fNO(i)=fNO(i)-1; fIDHDN(i)=fIDHDN(i)+1; 

i=i+1;
Rnames{i} = 'IDHNDOO1 = ICHNP + HO2';
k(:,i) = 1.256E13*exp(-10000./T);
Gstr{i,1} = 'IDHNDOO1';
fIDHNDOO1(i)=fIDHNDOO1(i)-1; fICHNP(i)=fICHNP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IDHNDOO2 + HO2 = 0.494IDHPN + 0.441HAC + 0.441ETHLN + 0.065MVK3OH4N + 0.065HCHO + 0.506OH + 0.506HO2';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'IDHNDOO2'; Gstr{i,2} = 'HO2'; 
fIDHNDOO2(i)=fIDHNDOO2(i)-1; fHO2(i)=fHO2(i)-0.494; fIDHPN(i)=fIDHPN(i)+0.494; fHAC(i)=fHAC(i)+0.441; fETHLN(i)=fETHLN(i)+0.441; fMVK3OH4N(i)=fMVK3OH4N(i)+0.065; fHCHO(i)=fHCHO(i)+0.065; fOH(i)=fOH(i)+0.506; 

i=i+1;
Rnames{i} = 'IDHNDOO2 + NO = 0.858HAC + 0.858ETHLN + 0.142MVK3OH4N + 0.142HCHO + HO2 + NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,2.258,11,1,0);
Gstr{i,1} = 'IDHNDOO2'; Gstr{i,2} = 'NO'; 
fIDHNDOO2(i)=fIDHNDOO2(i)-1; fNO(i)=fNO(i)-1; fHAC(i)=fHAC(i)+0.858; fETHLN(i)=fETHLN(i)+0.858; fMVK3OH4N(i)=fMVK3OH4N(i)+0.142; fHCHO(i)=fHCHO(i)+0.142; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'IDHNDOO2 + NO = IDHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,2.258,11,1,0);
Gstr{i,1} = 'IDHNDOO2'; Gstr{i,2} = 'NO'; 
fIDHNDOO2(i)=fIDHNDOO2(i)-1; fNO(i)=fNO(i)-1; fIDHDN(i)=fIDHDN(i)+1; 

i=i+1;
Rnames{i} = 'IDHNDOO2 = ICHNP + HO2';
k(:,i) = 5.092E12*exp(-10000./T);
Gstr{i,1} = 'IDHNDOO2'; 
fIDHNDOO2(i)=fIDHNDOO2(i)-1; fICHNP(i)=fICHNP(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IDHNBOO + HO2 = 0.379HO2 + 0.379OH + 0.621IDHPN + 0.094MACR2OH3N + 0.242GLYC + 0.242PROPNN + 0.010MVK3OH4N + 0.033HAC + 0.033ETHLN + 0.104HCHO';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'IDHNBOO'; Gstr{i,2} = 'HO2'; 
fIDHNBOO(i)=fIDHNBOO(i)-1; fHO2(i)=fHO2(i)-0.621; fOH(i)=fOH(i)+0.379; fIDHPN(i)=fIDHPN(i)+0.621; fMACR2OH3N(i)=fMACR2OH3N(i)+0.094; fGLYC(i)=fGLYC(i)+0.242; fPROPNN(i)=fPROPNN(i)+0.242; fMVK3OH4N(i)=fMVK3OH4N(i)+0.01; fHAC(i)=fHAC(i)+0.033; fETHLN(i)=fETHLN(i)+0.033; fHCHO(i)=fHCHO(i)+0.104; 

i=i+1;
Rnames{i} = 'IDHNBOO + NO = 0.355MACR2OH3N + 0.546PROPNN + 0.546GLYC + 0.028MVK3OH4N + 0.071ETHLN + 0.071HAC + HO2 + NO2 + 0.383HCHO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,1.851,11,1,0);
Gstr{i,1} = 'IDHNBOO'; Gstr{i,2} = 'NO'; 
fIDHNBOO(i)=fIDHNBOO(i)-1; fNO(i)=fNO(i)-1; fMACR2OH3N(i)=fMACR2OH3N(i)+0.355; fPROPNN(i)=fPROPNN(i)+0.546; fGLYC(i)=fGLYC(i)+0.546; fMVK3OH4N(i)=fMVK3OH4N(i)+0.0287; fETHLN(i)=fETHLN(i)+0.071; fHAC(i)=fHAC(i)+0.071; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; fHCHO(i)=fHCHO(i)+0.383; 

i=i+1;
Rnames{i} = 'IDHNBOO + NO = IDHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,1.851,11,1,0.);
Gstr{i,1} = 'IDHNBOO'; Gstr{i,2} = 'NO'; 
fIDHNBOO(i)=fIDHNBOO(i)-1; fNO(i)=fNO(i)-1; fIDHDN(i)=fIDHDN(i)+1; 

i=i+1;
Rnames{i} = 'INPB + OH = IHPNBOO';
k(:,i) = 0.5102*8.72E-12*exp(390./T) ; 
Gstr{i,1} = 'INPB'; Gstr{i,2} = 'OH'; 
fINPB(i)=fINPB(i)-1; fOH(i)=fOH(i)-1; fIHPNBOO(i)=fIHPNBOO(i)+1; 

i=i+1;
Rnames{i} = 'INPD + OH = IHPNDOO';
k(:,i) = 0.6778*2.37E-11*exp(390./T) ; 
Gstr{i,1} = 'INPD'; Gstr{i,2} = 'OH'; 
fINPD(i)=fINPD(i)-1; fOH(i)=fOH(i)-1; fIHPNDOO(i)=fIHPNDOO(i)+1; 

i=i+1;
Rnames{i} = 'INPB + OH = OH + IHNE';
k(:,i) = F0AM_isop_EPO(T,M,6.673E-12,390,2.28E-20);
Gstr{i,1} = 'INPB'; Gstr{i,2} = 'OH'; 
fINPB(i)=fINPB(i)-1; fIHNE(i)=fIHNE(i)+1;

i=i+1;
Rnames{i} = 'INPD + OH = OH + IHNE';
k(:,i) = F0AM_isop_EPO(T,M,8.77E-12,390,2.185E-20);
Gstr{i,1} = 'INPD'; Gstr{i,2} = 'OH'; 
fINPD(i)=fINPD(i)-1; fIHNE(i)=fIHNE(i)+1; 

i=i+1;
Rnames{i} = 'INPD + OH = NO2 + ICHE';
k(:,i) = F0AM_isop_EPO(T,M,1.493E-11,390,2.715E-19);
Gstr{i,1} = 'INPD'; Gstr{i,2} = 'OH'; 
fINPD(i)=fINPD(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1; fICHE(i)=fICHE(i)+1; 

i=i+1;
Rnames{i} = 'INPB + OH = INO2B';
k(:,i) = 3.4E-12*exp(200./T);
Gstr{i,1} = 'INPB'; Gstr{i,2} = 'OH'; 
fINPB(i)=fINPB(i)-1; fOH(i)=fOH(i)-1; fINO2B(i)=fINO2B(i)+1; 

i=i+1;
Rnames{i} = 'INPD + OH = INO2D';
k(:,i) = 3.4E-12*exp(200./T);
Gstr{i,1} = 'INPD'; Gstr{i,2} = 'OH'; 
fINPD(i)=fINPD(i)-1; fOH(i)=fOH(i)-1; fINO2D(i)=fINO2D(i)+1; 

i=i+1;
Rnames{i} = 'INPD + OH = 0.843ISOP1N4CO + 0.157ISOP1CO4N + OH';
k(:,i) = 7.5E-12*exp(20./T);
Gstr{i,1} = 'INPD'; Gstr{i,2} = 'OH'; 
fINPD(i)=fINPD(i)-1; fISOP1N4CO(i)=fISOP1N4CO(i)+0.843; fISOP1CO4N(i)=fISOP1CO4N(i)+0.157; 

i=i+1;
Rnames{i} = 'ISOP1N4OH = NO2 + 0.45HO2 + 0.55MVK3OOH4OH + 0.55CO + 0.55OH + 0.45ISOP1CO4OH';
k(:,i) = J53 ; 
Gstr{i,1} = 'ISOP1N4OH'; 
fISOP1N4OH(i)=fISOP1N4OH(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.45; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.55; fCO(i)=fCO(i)+0.55; fOH(i)=fOH(i)+0.55; fISOP1CO4OH(i)=fISOP1CO4OH(i)+0.45; 

i=i+1;
Rnames{i} = 'ISOP1OH4N = NO2 + 0.45HO2 + 0.45ISOP1OH4CO + 0.55CO + 0.55OH + 0.55MACR2OOH3OH';
k(:,i) = J53 ;  
Gstr{i,1} = 'ISOP1OH4N'; 
fISOP1OH4N(i)=fISOP1OH4N(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.45; fISOP1OH4CO(i)=fISOP1OH4CO(i)+0.45; fCO(i)=fCO(i)+0.55; fOH(i)=fOH(i)+0.55; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.55; 

i=i+1;
Rnames{i} = 'ISOP1OH2N = NO2 + MVK + HO2 + HCHO';
k(:,i) = J55 ; 
Gstr{i,1} = 'ISOP1OH2N'; 
fISOP1OH2N(i)=fISOP1OH2N(i)-1; fNO2(i)=fNO2(i)+1; fMVK(i)=fMVK(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3N4OH = NO2 + MACR + HO2 + HCHO';
k(:,i) = J54 ; 
Gstr{i,1} = 'ISOP3N4OH'; 
fISOP3N4OH(i)=fISOP3N4OH(i)-1; fNO2(i)=fNO2(i)+1; fMACR(i)=fMACR(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'INPB = NO2 + HCHO + 0.097MACR + 0.903MVK + 0.67OH + 0.33HO2';
k(:,i) = J55 + J41 ; % SUN*8.3E-6 ;
Gstr{i,1} = 'INPB'; 
fINPB(i)=fINPB(i)-1; fNO2(i)=fNO2(i)+1; fHCHO(i)=fHCHO(i)+1; fMACR(i)=fMACR(i)+0.097; fMVK(i)=fMVK(i)+0.903; fOH(i)=fOH(i)+0.67; fHO2(i)=fHO2(i)+0.33; 

i=i+1;
Rnames{i} = 'INPD = 0.217NO2 + 0.182IHOO1 + 0.034IHOO4 + 0.063ISOP1CO4N + 0.062ISOP1N4CO + 0.125HO2 + 0.659INO + 0.783OH ';
k(:,i) = J53 + J41 ; % SUN*8.3E-6 ;
Gstr{i,1} = 'INPD'; 
fINPD(i)=fINPD(i)-1; fNO2(i)=fNO2(i)+1; fIHOO1(i)=fIHOO1(i)+0.182; fIHOO4(i)=fIHOO4(i)+0.034; fISOP1CO4N(i)=fISOP1CO4N(i)+0.063; fISOP1N4CO(i)=fISOP1N4CO(i)+0.062; fHO2(i)=fHO2(i)+0.125; fINO(i)=fINO(i)+0.659; fOH(i)=fOH(i)+0.783; 

i=i+1;
Rnames{i} = 'IHNB = NO2 + HO2 + HCHO + 0.5MVK + 0.5MACR';
k(:,i) = J55 ; % SUN*1.8E-6 ; 
Gstr{i,1} = 'IHNB'; 
fIHNB(i)=fIHNB(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fMVK(i)=fMVK(i)+0.5; fMACR(i)=fMACR(i)+0.5; 

i=i+1;
Rnames{i} = 'IHPNDOO = OH + ICHNP';
k(:,i) = 6.55E12*exp(-10000./T);
Gstr{i,1} = 'IHPNDOO'; 
fIHPNDOO(i)=fIHPNDOO(i)-1; fOH(i)=fOH(i)+1; fICHNP(i)=fICHNP(i)+1; 

i=i+1;
Rnames{i} = 'IHPNBOO = OH + 0.5ICHNP + 0.5IHNPE';
k(:,i) = 8.72E12*exp(-10000./T);
Gstr{i,1} = 'IHPNBOO'; 
fIHPNBOO(i)=fIHPNBOO(i)-1; fOH(i)=fOH(i)+1; fICHNP(i)=fICHNP(i)+0.5; fIHNPE(i)=fIHNPE(i)+0.5; 

i=i+1;
Rnames{i} = 'IHPNBOO + HO2 = 0.234IHNDP + 0.060MACR2OOH3N + 0.340GLYC + 0.249HPETHNL + 0.042HAC + 0.008MVK3OOH4N + 0.009HPAC + 0.054MVK3OOH4OH + 0.004MACR2OOH3OH + 1.147OH + 0.326HO2 + 0.126HCHO + 0.589PROPNN + 0.051ETHLN + 0.058NO2';
k(:,i) = 2.64E-13*exp(1300./T); 
Gstr{i,1} = 'IHPNBOO'; Gstr{i,2} = 'HO2'; 
fIHPNBOO(i)=fIHPNBOO(i)-1; fHO2(i)=fHO2(i)-0.674; fIHNDP(i)=fIHNDP(i)+0.234; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.06; fGLYC(i)=fGLYC(i)+0.34; fHPETHNL(i)=fHPETHNL(i)+0.249; fHAC(i)=fHAC(i)+0.042; fMVK3OOH4N(i)=fMVK3OOH4N(i)+0.008; fHPAC(i)=fHPAC(i)+0.009; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.054; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.004; fOH(i)=fOH(i)+1.147; fHCHO(i)=fHCHO(i)+0.126; fPROPNN(i)=fPROPNN(i)+0.589; fETHLN(i)=fETHLN(i)+0.051; fNO2(i)=fNO2(i)+0.058;

i=i+1;
Rnames{i} = 'IHPNDOO + HO2 = 0.387IHNDP + 0.050MACR2OOH3N + 0.471HPETHNL + 0.023MACR2OH3N + 0.005MVK3OOH4N + 0.010MVK3OH4N + 0.054HPAC + 0.646OH + 0.580HO2 + 0.088HCHO + 0.471PROPNN + 0.054ETHLN';
k(:,i) = 2.64E-13*exp(1300./T); 
Gstr{i,1} = 'IHPNDOO'; Gstr{i,2} = 'HO2'; 
fIHPNDOO(i)=fIHPNDOO(i)-1; fHO2(i)=fHO2(i)-0.42; fIHNDP(i)=fIHNDP(i)+0.387; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.05; fHPETHNL(i)=fHPETHNL(i)+0.471; fMACR2OH3N(i)=fMACR2OH3N(i)+0.023; fMVK3OH4N(i)=fMVK3OH4N(i)+0.005; fMVK3OH4N(i)=fMVK3OH4N(i)+0.01; fHPAC(i)=fHPAC(i)+0.054; fOH(i)=fOH(i)+0.646; fHCHO(i)=fHCHO(i)+0.088; fPROPNN(i)=fPROPNN(i)+0.471; fETHLN(i)=fETHLN(i)+0.054; 

i=i+1;
Rnames{i} = 'IHPNBOO + NO = 0.384GLYC + 0.170MACR2OOH3N + 0.303HPETHNL + 0.014MVK3OOH4N + 0.051HAC + 0.013HPAC + 0.059MVK3OOH4OH + 0.006MACR2OOH3OH + 0.687PROPNN + 0.064ETHLN + 0.249HCHO + 1.065NO2 + 0.500HO2 + 0.435OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,6.092,12,1,0);
Gstr{i,1} = 'IHPNBOO'; Gstr{i,2} = 'NO'; 
fIHPNBOO(i)=fIHPNBOO(i)-1; fNO(i)=fNO(i)-1; fGLYC(i)=fGLYC(i)+0.384; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.17; fHPETHNL(i)=fHPETHNL(i)+0.303; fMVK3OOH4N(i)=fMVK3OOH4N(i)+0.014; fHAC(i)=fHAC(i)+0.051; fHPAC(i)=fHPAC(i)+0.013; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.059; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.006; fPROPNN(i)=fPROPNN(i)+0.687; fETHLN(i)=fETHLN(i)+0.064; fHCHO(i)=fHCHO(i)+0.249; fNO2(i)=fNO2(i)+1.064; fHO2(i)=fHO2(i)+0.5; fOH(i)=fOH(i)+0.435; 

i=i+1;
Rnames{i} = 'IHPNBOO + NO = IHPDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,6.092,12,1,0);
Gstr{i,1} = 'IHPNBOO'; Gstr{i,2} = 'NO'; 
fIHPNBOO(i)=fIHPNBOO(i)-1; fNO(i)=fNO(i)-1; fIHPDN(i)=fIHPDN(i)+1; 

i=i+1;
Rnames{i} = 'IHPNDOO + NO = 0.221MACR2OOH3N + 0.070MACR2OH3N + 0.590HPETHNL + 0.023MVK3OOH4N + 0.070HPAC + 0.026MVK3OH4N + 0.590PROPNN + 0.070ETHLN + 0.340HCHO + 1.000NO2 + 0.904HO2 + 0.096OH';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,4.383,12,1,0);
Gstr{i,1} = 'IHPNDOO'; Gstr{i,2} = 'NO'; 
fIHPNDOO(i)=fIHPNDOO(i)-1; fNO(i)=fNO(i)-1; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.221; fMACR2OH3N(i)=fMACR2OH3N(i)+0.07; fHPETHNL(i)=fHPETHNL(i)+0.59; fMVK3OOH4N(i)=fMVK3OOH4N(i)+0.023; fHPAC(i)=fHPAC(i)+0.07; fMVK3OH4N(i)=fMVK3OH4N(i)+0.026; fPROPNN(i)=fPROPNN(i)+0.59; fETHLN(i)=fETHLN(i)+0.07; fHCHO(i)=fHCHO(i)+0.34; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+0.904; fOH(i)=fOH(i)+0.096; 

i=i+1;
Rnames{i} = 'IHPNDOO + NO = IHPDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,4.383,12,1,0);
Gstr{i,1} = 'IHPNDOO'; Gstr{i,2} = 'NO'; 
fIHPNDOO(i)=fIHPNDOO(i)-1; fNO(i)=fNO(i)-1; fIHPDN(i)=fIHPDN(i)+1;

i=i+1;
Rnames{i} = 'IHNE + OH = 0.201IHNEOO + 0.388ICN1OO + 0.026ICN2OO + 0.172ICN3OO + 0.026MCRENOL + 0.187MVKENOL + 0.213NO2 + 0.213HCHO';
k(:,i) = F0AM_isop_EPO(T,M,3.22E-11,-400,1.014E-20);
Gstr{i,1} = 'IHNE'; Gstr{i,2} = 'OH'; 
fIHNE(i)=fIHNE(i)-1; fOH(i)=fOH(i)-1; fIHNEOO(i)=fIHNEOO(i)+0.201; fICN1OO(i)=fICN1OO(i)+0.388; fICN2OO(i)=fICN2OO(i)+0.026; fICN3OO(i)=fICN3OO(i)+0.172; fMCRENOL(i)=fMCRENOL(i)+0.026; fMVKENOL(i)=fMVKENOL(i)+0.187; fNO2(i)=fNO2(i)+0.213; fHCHO(i)=fHCHO(i)+0.213; 

i=i+1;
Rnames{i} = 'IHNE + OH = 0.201IHNEOO + 0.200ICN1OO + 0.599ICNE';
k(:,i) = 7.77E-12*exp(-400./T); 
Gstr{i,1} = 'IHNE'; Gstr{i,2} = 'OH'; 
fIHNE(i)=fIHNE(i)-1; fOH(i)=fOH(i)-1; fIHNEOO(i)=fIHNEOO(i)+0.201; fICN1OO(i)=fICN1OO(i)+0.2; fICNE(i)=fICNE(i)+0.599; 

i=i+1;
Rnames{i} = 'IHNEOO + HO2 = 0.8ICHNP + 0.2PROPNN + 0.2CO + 0.2HO2 + 0.2HCHO + 0.2OH';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'IHNEOO'; Gstr{i,2} = 'HO2'; 
fIHNEOO(i)=fIHNEOO(i)-1; fHO2(i)=fHO2(i)-0.8; fICHNP(i)=fICHNP(i)+0.8; fPROPNN(i)=fPROPNN(i)+0.2; fCO(i)=fCO(i)+0.2; fHCHO(i)=fHCHO(i)+0.2; fOH(i)=fOH(i)+0.2; 

i=i+1;
Rnames{i} = 'IHNEOO + NO = ICHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,15.73,11,1,0);
Gstr{i,1} = 'IHNEOO'; Gstr{i,2} = 'NO'; 
fIHNEOO(i)=fIHNEOO(i)-1; fNO(i)=fNO(i)-1; fICHDN(i)=fICHDN(i)+1;

i=i+1;
Rnames{i} = 'IHNEOO + NO = PROPNN + CO + HO2 + HCHO + NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,15.73,11,1,0);
Gstr{i,1} = 'IHNEOO'; Gstr{i,2} = 'NO'; 
fIHNEOO(i)=fIHNEOO(i)-1; fNO(i)=fNO(i)-1; fPROPNN(i)=fPROPNN(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1N4CO + OH = NO2 + ICHE';
k(:,i) = F0AM_isop_EPO(T,M,3.35E-12,390,2.715E-19);
Gstr{i,1} = 'ISOP1N4CO'; Gstr{i,2} = 'OH'; 
fISOP1N4CO(i)=fISOP1N4CO(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1; fICHE(i)=fICHE(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1N4CO + OH = 0.660CO + 0.660HO2 + 0.660MACR2OOH3N + 0.340ICN1OO';
k(:,i) = 4.42E-12*exp(390./T); 
Gstr{i,1} = 'ISOP1N4CO'; Gstr{i,2} = 'OH'; 
fISOP1N4CO(i)=fISOP1N4CO(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+0.66; fHO2(i)=fHO2(i)+0.66; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.66; fICN1OO(i)=fICN1OO(i)+0.34; 

i=i+1;
Rnames{i} = 'ISOP1CO4N + OH = 0.281CO + 0.281HO2 + 0.281MVK3OOH4N + 0.719ICN2OO';
k(:,i) = 6.48E-12*exp(390./T); 
Gstr{i,1} = 'ISOP1CO4N'; Gstr{i,2} = 'OH'; 
fISOP1CO4N(i)=fISOP1CO4N(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+0.281; fHO2(i)=fHO2(i)+0.281; fMVK3OOH4N(i)=fMVK3OOH4N(i)+0.281; fICN2OO(i)=fICN2OO(i)+0.719; 

i=i+1;
Rnames{i} = 'ISOP1CO4N + OH = NO2 + ICHE';
k(:,i) = F0AM_isop_EPO(T,M,2.09E-12,390,2.715E-19);
Gstr{i,1} = 'ISOP1CO4N'; Gstr{i,2} = 'OH'; 
fISOP1CO4N(i)=fISOP1CO4N(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1; fICHE(i)=fICHE(i)+1; 

i=i+1;
Rnames{i} = 'ISOP3CO4N + OH = ICN3OO';
k(:,i) = 2.7E-12*exp(390./T); 
Gstr{i,1} = 'ISOP3CO4N'; Gstr{i,2} = 'OH'; 
fISOP3CO4N(i)=fISOP3CO4N(i)-1; fOH(i)=fOH(i)-1; fICN3OO(i)=fICN3OO(i)+1;

i=i+1;
Rnames{i} = 'ICN1OO = MACR2OH3N + CO + OH';
k(:,i) = 1.0E7*exp(-5000./T); 
Gstr{i,1} = 'ICN1OO';
fICN1OO(i)=fICN1OO(i)-1; fMACR2OH3N(i)=fMACR2OH3N(i)+1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+1; 

i=i+1;
Rnames{i} = 'ICN1OO + HO2 = 0.250ICHNP + 0.563MACR2OH3N + 0.563CO + 0.187PROPNN + 0.187GLYC + 0.75OH + 0.75HO2';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'ICN1OO'; Gstr{i,2} = 'HO2'; 
fICN1OO(i)=fICN1OO(i)-1; fHO2(i)=fHO2(i)-0.25; fICHNP(i)=fICHNP(i)+0.25; fMACR2OH3N(i)=fMACR2OH3N(i)+0.563; fCO(i)=fCO(i)+0.563; fPROPNN(i)=fPROPNN(i)+0.187; fGLYC(i)=fGLYC(i)+0.187; fOH(i)=fOH(i)+0.75; 

i=i+1;
Rnames{i} = 'ICN1OO + NO = NO2 + HO2 + 0.25PROPNN + 0.25GLYX + 0.75MACR2OH3N + 0.75CO';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,22.837,11,1,0);
Gstr{i,1} = 'ICN1OO'; Gstr{i,2} = 'NO'; 
fICN1OO(i)=fICN1OO(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fPROPNN(i)=fPROPNN(i)+0.25; fGLYX(i)=fGLYX(i)+0.25; fMACR2OH3N(i)=fMACR2OH3N(i)+0.75; fCO(i)=fCO(i)+0.75; 

i=i+1;
Rnames{i} = 'ICN1OO + NO = ICHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,22.837,11,1,0.);
Gstr{i,1} = 'ICN1OO'; Gstr{i,2} = 'NO'; 
fICN1OO(i)=fICN1OO(i)-1; fNO(i)=fNO(i)-1; fICHDN(i)=fICHDN(i)+1;

i=i+1;
Rnames{i} = 'ICN2OO = MVK3OH4N + CO + OH';
k(:,i) = 1.0E7*exp(-5000./T); 
Gstr{i,1} = 'ICN2OO';
fICN2OO(i)=fICN2OO(i)-1; fMVK3OH4N(i)=fMVK3OH4N(i)+1; fCO(i)=fCO(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ICN2OO + HO2 = 0.15ICHNP + 0.638MVK3OH4N + 0.638CO + 0.212MGLY + 0.212ETHLN + 0.85OH + 0.85HO2';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'ICN2OO'; Gstr{i,2} = 'HO2'; 
fICN2OO(i)=fICN2OO(i)-1; fHO2(i)=fHO2(i)-0.15; fICHNP(i)=fICHNP(i)+0.15; fMVK3OH4N(i)=fMVK3OH4N(i)+0.638; fCO(i)=fCO(i)+0.638; fMGLY(i)=fMGLY(i)+0.212; fETHLN(i)=fETHLN(i)+0.212; fOH(i)=fOH(i)+0.85; 

i=i+1;
Rnames{i} = 'ICN2OO + NO = 0.75MVK3OH4N + NO2 + HO2 + 0.75CO + 0.25MGLY + 0.25ETHLN';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,18.181,11,1,0);
Gstr{i,1} = 'ICN2OO'; Gstr{i,2} = 'NO'; 
fICN2OO(i)=fICN2OO(i)-1; fNO(i)=fNO(i)-1; fMVK3OH4N(i)=fMVK3OH4N(i)+0.75; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fCO(i)=fCO(i)+0.75; fMGLY(i)=fMGLY(i)+0.25; fETHLN(i)=fETHLN(i)+0.25; 

i=i+1;
Rnames{i} = 'ICN2OO + NO = ICHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,18.181,11,1,0);
Gstr{i,1} = 'ICN2OO'; Gstr{i,2} = 'NO'; 
fICN2OO(i)=fICN2OO(i)-1; fNO(i)=fNO(i)-1; fICHDN(i)=fICHDN(i)+1; 

i=i+1;
Rnames{i} = 'ICN3OO + HO2 = 0.15ICHNP + 0.638HAC + 0.85OH + 0.638CH3CO3 + 0.638NO2 + 0.212HO2 + 0.212HCHO + 0.212MVK3CO4N';
k(:,i) = 2.6E-13*exp(1300./T); 
Gstr{i,1} = 'ICN3OO'; Gstr{i,2} = 'HO2'; 
fICN3OO(i)=fICN3OO(i)-1; fHO2(i)=fHO2(i)-0.788; fICHNP(i)=fICHNP(i)+0.15; fHAC(i)=fHAC(i)+0.638; fOH(i)=fOH(i)+0.85; fCH3CO3(i)=fCH3CO3(i)+0.638; fNO2(i)=fNO2(i)+0.638; fHCHO(i)=fHCHO(i)+0.212; fMVK3CO4N(i)=fMVK3CO4N(i)+0.212; 

i=i+1;
Rnames{i} = 'ICN3OO + NO = ICHDN';
k(:,i) = F0AM_isop_NIT(T,M,2.7E-12,350,18.181,11,1,0);
Gstr{i,1} = 'ICN3OO'; Gstr{i,2} = 'NO'; 
fICN3OO(i)=fICN3OO(i)-1; fNO(i)=fNO(i)-1; fICHDN(i)=fICHDN(i)+1;

i=i+1;
Rnames{i} = 'ICN3OO + NO = 0.25CH3CO3 + 0.25HAC + 0.75MVK3CO4N + 0.75HO2 + 0.75HCHO + 1.25NO2';
k(:,i) = F0AM_isop_ALK(T,M,2.7E-12,350,18.181,11,1,0);
Gstr{i,1} = 'ICN3OO'; Gstr{i,2} = 'NO'; 
fICN3OO(i)=fICN3OO(i)-1; fNO(i)=fNO(i)-1; fCH3CO3(i)=fCH3CO3(i)+0.25; fHAC(i)=fHAC(i)+0.25; fMVK3CO4N(i)=fMVK3CO4N(i)+0.75; fHO2(i)=fHO2(i)+0.75; fHCHO(i)=fHCHO(i)+0.75; fNO2(i)=fNO2(i)+1.25; 

i=i+1;
Rnames{i} = 'ISOP1N4CO = NO2 + 0.818CO + 0.637OH + 0.182HO2 + 0.182ISOP1CO4CO + 0.182MVK3OOH4CO + 0.455MCRENOL + 0.182C4HVP2';
k(:,i) = J56 ; % SUN*2.2E-4;  
Gstr{i,1} = 'ISOP1N4CO'; 
fISOP1N4CO(i)=fISOP1N4CO(i)-1; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.818; fOH(i)=fOH(i)+0.637; fHO2(i)=fHO2(i)+0.182; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.182; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.182; fMCRENOL(i)=fMCRENOL(i)+0.455; fC4HVP2(i)=fC4HVP2(i)+0.182; 

i=i+1;
Rnames{i} = 'ISOP1CO4N = NO2 + 0.888CO + 0.662OH + 0.112HO2 + 0.112ISOP1CO4CO + 0.112MVK3OOH4CO + 0.552MVKENOL + 0.224C4HVP1';
k(:,i) = J56 ; % SUN*2.32E-4;  
Gstr{i,1} = 'ISOP1CO4N'; 
fISOP1CO4N(i)=fISOP1CO4N(i)-1; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.888; fOH(i)=fOH(i)+0.662; fHO2(i)=fHO2(i)+0.112; fISOP1CO4CO(i)=fISOP1CO4CO(i)+0.112; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.112; fMVKENOL(i)=fMVKENOL(i)+0.552; fC4HVP1(i)=fC4HVP1(i)+0.224; 

i=i+1;
Rnames{i} = 'ISOP3CO4N = NO2 + HCHO + MACR1OO';
k(:,i) = J56 ; % SUN*2.24E-4;  
Gstr{i,1} = 'ISOP3CO4N'; 
fISOP3CO4N(i)=fISOP3CO4N(i)-1; fNO2(i)=fNO2(i)+1; fHCHO(i)=fHCHO(i)+1; fMACR1OO(i)=fMACR1OO(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1N4CO + OH = ICN4OO';
k(:,i) = 3.3E-12*exp(470./T) ; 
Gstr{i,1} = 'ISOP1N4CO'; Gstr{i,2} = 'OH'; 
fISOP1N4CO(i)=fISOP1N4CO(i)-1; fOH(i)=fOH(i)-1; fICN4OO(i)=fICN4OO(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1CO4N + OH = ICN5OO';
k(:,i) = 3.3E-12*exp(470./T) ; 
Gstr{i,1} = 'ISOP1CO4N'; Gstr{i,2} = 'OH'; 
fISOP1CO4N(i)=fISOP1CO4N(i)-1; fOH(i)=fOH(i)-1; fICN5OO(i)=fICN5OO(i)+1;  

i=i+1;
Rnames{i} = 'ICN4OO + NO = 0.67ICN4OO + 0.33CO2 + 0.33CO + 0.33HO2 + 0.33PROPNN + NO2';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'ICN4OO'; Gstr{i,2} = 'NO'; 
fICN4OO(i)=fICN4OO(i)-0.33; fNO(i)=fNO(i)-1; fCO(i)=fCO(i)+0.33; fHO2(i)=fHO2(i)+0.33; fPROPNN(i)=fPROPNN(i)+0.33; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ICN4OO + HO2 = 0.67ICN4OO + 0.33CO2 + 0.33CO + 0.33HO2 + 0.33PROPNN + OH';
k(:,i) = 2.54E-13*exp(1300./T); 
Gstr{i,1} = 'ICN4OO'; Gstr{i,2} = 'HO2'; 
fICN4OO(i)=fICN4OO(i)-0.33; fHO2(i)=fHO2(i)-0.67; fCO(i)=fCO(i)+0.33; fPROPNN(i)=fPROPNN(i)+0.33; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ICN5OO + NO = 0.67ICN5OO + 0.33CO2 + 0.33CH3CO3 + 0.33ETHLN + NO2';
k(:,i) = 2.7E-12*exp(350./T); 
Gstr{i,1} = 'ICN5OO'; Gstr{i,2} = 'NO'; 
fICN5OO(i)=fICN5OO(i)-0.33; fNO(i)=fNO(i)-1; fCH3CO3(i)=fCH3CO3(i)+0.33; fETHLN(i)=fETHLN(i)+0.33; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ICN5OO + HO2 = 0.67ICN5OO + 0.33CO2 + 0.33CH3CO3 + 0.33ETHLN + OH';
k(:,i) = 2.54E-13*exp(1300./T); 
Gstr{i,1} = 'ICN5OO'; Gstr{i,2} = 'HO2'; 
fICN5OO(i)=fICN5OO(i)-0.33; fHO2(i)=fHO2(i)-1; fCH3CO3(i)=fCH3CO3(i)+0.33; fETHLN(i)=fETHLN(i)+0.33; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'ICHE + OH = OH + 1.5CO + 0.5HCHO + 0.5MGLY + 0.5HAC';
k(:,i) = 9.85E-12*exp(410./T) ;
Gstr{i,1} = 'ICHE'; Gstr{i,2} = 'OH'; 
fICHE(i)=fICHE(i)-1; fCO(i)=fCO(i)+1.5; fHCHO(i)=fHCHO(i)+0.5; fMGLY(i)=fMGLY(i)+0.5; fHAC(i)=fHAC(i)+0.5;  

i=i+1;
Rnames{i} = 'ICNE + OH = OH + 2.000CO + 0.350PROPNN + 0.650MGLY + 0.650HO2 + 0.650NO2';
k(:,i) = 9.85E-12*exp(410./T) ;
Gstr{i,1} = 'ICNE'; Gstr{i,2} = 'OH'; 
fICNE(i)=fICNE(i)-1; fCO(i)=fCO(i)+2; fPROPNN(i)=fPROPNN(i)+0.35; fMGLY(i)=fMGLY(i)+0.65; fHO2(i)=fHO2(i)+0.65; fNO2(i)=fNO2(i)+0.65; 

i=i+1;
Rnames{i} = 'MVK3N4OH = CH3CO3 + GLYC + NO2';
k(:,i) = J56 ; % SUN*6.46E-5; 
Gstr{i,1} = 'MVK3N4OH'; 
fMVK3N4OH(i)=fMVK3N4OH(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fGLYC(i)=fGLYC(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3N4OH + OH = 0.650HCOOH + NO3 + 0.650MGLY + 0.350HCHO + 0.350PYRAC';
k(:,i) = 1.5E-12*exp(380/T) ;
Gstr{i,1} = 'MVK3N4OH'; Gstr{i,2} = 'OH'; 
fMVK3N4OH(i)=fMVK3N4OH(i)-1; fOH(i)=fOH(i)-1; fHCOOH(i)=fHCOOH(i)+0.65; fNO3(i)=fNO3(i)+1; fHCHO(i)=fHCHO(i)+0.35; fMGLY(i)=fMGLY(i)+0.65; fPYRAC(i)=fPYRAC(i)+0.35; 

i=i+1;
Rnames{i} = 'MVK3OH4N = CH3CO3 + ETHLN + HO2';
k(:,i) = J22+J53 ; % SUN*4.21E-5; 
Gstr{i,1} = 'MVK3OH4N'; 
fMVK3OH4N(i)=fMVK3OH4N(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fETHLN(i)=fETHLN(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OH4N + OH = MVK3OH4CO + NO2';
k(:,i) = 2.23E-12 ;
Gstr{i,1} = 'MVK3OH4N'; Gstr{i,1} = 'OH'; 
fMVK3OH4N(i)=fMVK3OH4N(i)-1; fOH(i)=fOH(i)-1; fMVK3OH4CO(i)=fMVK3OH4CO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OOH4OH = CH3CO3 + GLYC + OH';
k(:,i) = J22 + J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MVK3OOH4OH'; 
fMVK3OOH4OH(i)=fMVK3OOH4OH(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fGLYC(i)=fGLYC(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OOH4OH + OH = MVK3CO4OH + OH';
k(:,i) = 5.77E-11 ;
Gstr{i,1} = 'MVK3OOH4OH'; Gstr{i,2} = 'OH'; 
fMVK3OOH4OH(i)=fMVK3OOH4OH(i)-1; fMVK3CO4OH(i)=fMVK3CO4OH(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OH4OOH + OH = MVK3OH4CO + OH';
k(:,i) = 5.77E-11 ;
Gstr{i,1} = 'MVK3OH4OOH'; Gstr{i,2} = 'OH'; 
fMVK3OH4OOH(i)=fMVK3OH4OOH(i)-1; fMVK3OH4CO(i)=fMVK3OH4CO(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OH4OOH = OH + HO2 + HCHO + MGLY';
k(:,i) = J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MVK3OH4OOH'; 
fMVK3OH4OOH(i)=fMVK3OH4OOH(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fMGLY(i)=fMGLY(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OOH3OH = OH + CO + HO2 + HAC';
k(:,i) = J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MACR2OOH3OH';  
fMACR2OOH3OH(i)=fMACR2OOH3OH(i)-1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fHAC(i)=fHAC(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OOH3OH + OH = CO + OH + HAC';
k(:,i) = 2.7E-12*exp(470./T); 
Gstr{i,1} = 'MACR2OOH3OH'; Gstr{i,2} = 'OH'; 
fMACR2OOH3OH(i)=fMACR2OOH3OH(i)-1; fCO(i)=fCO(i)+1; fHAC(i)=fHAC(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OH3OOH = OH + HCHO + HO2 + MGLY';
k(:,i) = J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MACR2OH3OOH'; 
fMACR2OH3OOH(i)=fMACR2OH3OOH(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fMGLY(i)=fMGLY(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OH3OOH + OH = CO2 + OH + HPAC';
k(:,i) = 2.7E-12*exp(470./T); 
Gstr{i,1} = 'MACR2OH3OOH'; Gstr{i,2} = 'OH'; 
fMACR2OH3OOH(i)=fMACR2OH3OOH(i)-1; fHPAC(i)=fHPAC(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OH3N + OH = CO2 + OH + PROPNN';
k(:,i) = 2.7E-12*exp(470./T); 
Gstr{i,1} = 'MACR2OH3N'; Gstr{i,2} = 'OH'; 
fMACR2OH3N(i)=fMACR2OH3N(i)-1; fPROPNN(i)=fPROPNN(i)+1; 

i=i+1;
Rnames{i} = 'MACR2N3OH = HAC + CO + HO2 + NO2';
k(:,i) = J56 ; % SUN*6.46E-5; 
Gstr{i,1} = 'MACR2N3OH';
fMACR2N3OH(i)=fMACR2N3OH(i)-1; fHAC(i)=fHAC(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MACR2N3OH + OH = MACRNO2';
k(:,i) = 1.39E-11*exp(380/T);
Gstr{i,1} = 'MACR2N3OH'; Gstr{i,2} = 'OH'; 
fMACR2N3OH(i)=fMACR2N3OH(i)-1; fOH(i)=fOH(i)-1; fMACRNO2(i)=fMACRNO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OOH4N = CH3CO3 + OH + ETHLN';
k(:,i) = J53 + J22 ; % SUN*4.21E-5; 
Gstr{i,1} = 'MVK3OOH4N'; 
fMVK3OOH4N(i)=fMVK3OOH4N(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fOH(i)=fOH(i)+1; fETHLN(i)=fETHLN(i)+1; 

i=i+1;
Rnames{i} = 'MVK3OOH4N + OH = OH + MVK3CO4N';
k(:,i) = 5.77E-11 ;
Gstr{i,1} = 'MVK3OOH4N'; Gstr{i,2} = 'OH'; 
fMVK3OOH4N(i)=fMVK3OOH4N(i)-1; fMVK3CO4N(i)=fMVK3CO4N(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OOH3N + OH = CO + OH + PROPNN';
k(:,i) = 2.7E-12*exp(470./T); 
Gstr{i,1} = 'MACR2OOH3N'; Gstr{i,2} = 'OH'; 
fMACR2OOH3N(i)=fMACR2OOH3N(i)-1; fCO(i)=fCO(i)+1; fPROPNN(i)=fPROPNN(i)+1; 

i=i+1;
Rnames{i} = 'MACR2OOH3N = PROPNN + OH + CO + HO2';
k(:,i) =  J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MACR2OOH3N'; 
fMACR2OOH3N(i)=fMACR2OOH3N(i)-1; fPROPNN(i)=fPROPNN(i)+1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'MVK3CO4N = 2CH3CO3 + NO2';
k(:,i) = J56 + J35 ; % 2.5E-4*SUN;
Gstr{i,1} = 'MVK3CO4N'; 
fMVK3CO4N(i)=fMVK3CO4N(i)-1; fCH3CO3(i)=fCH3CO3(i)+2; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ISOP1CO4CO + OH = CO + HO2 + MVK3OOH4CO';
k(:,i) = 3.0E-12*exp(650./T) ;
Gstr{i,1} = 'ISOP1CO4CO'; Gstr{i,2} = 'OH'; 
fISOP1CO4CO(i)=fISOP1CO4CO(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+1; 

i=i+1;
Rnames{i} = 'HPETHNL + OH = CO + OH + HCHO';
k(:,i) = 1.55E-12*exp(340./T) ;
Gstr{i,1} = 'HPETHNL'; Gstr{i,2} = 'OH'; 
fHPETHNL(i)=fHPETHNL(i)-1; fCO(i)=fCO(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'HPETHNL + OH = GLYX + OH';
k(:,i) = 2.91E-11 ;
Gstr{i,1} = 'HPETHNL'; Gstr{i,2} = 'OH'; 
fHPETHNL(i)=fHPETHNL(i)-1; fGLYX(i)=fGLYX(i)+1; 

i=i+1;
Rnames{i} = 'HPETHNL = OH + CO + HO2 + HCHO';
k(:,i) = J22+J41 ; % 2.5E-4*SUN;
Gstr{i,1} = 'HPETHNL'; 
fHPETHNL(i)=fHPETHNL(i)-1; fOH(i)=fOH(i)+1; fCO(i)=fCO(i)+1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HMML + OH = 0.7MGLY + 0.7OH + 0.3CH3CO3 + 0.3HCOOH';
k(:,i) = 4.33E-12 ;
Gstr{i,1} = 'HMML'; Gstr{i,2} = 'OH'; 
fHMML(i)=fHMML(i)-1; fOH(i)=fOH(i)-0.3; fMGLY(i)=fMGLY(i)+0.7; fCH3CO3(i)=fCH3CO3(i)+0.3; fHCOOH(i)=fHCOOH(i)+0.3; 

i=i+1;
Rnames{i} = 'MACR1OOH = OH + 1.65CO2 + 0.65CH3OO + HCHO + 0.35CH3CO3';
k(:,i) = J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'MACR1OOH'; 
fMACR1OOH(i)=fMACR1OOH(i)-1; fOH(i)=fOH(i)+1; fCH3OO(i)=fCH3OO(i)+0.65; fHCHO(i)=fHCHO(i)+1; fCH3CO3(i)=fCH3CO3(i)+0.35; 

i=i+1;
Rnames{i} = 'MACR1OH + OH = 1.65CO2 + HCHO + 0.35CH3CO3 + 0.65CH3OO';
k(:,i) = 1.51E-11 ;
Gstr{i,1} = 'MACR1OH'; Gstr{i,2} = 'OH'; 
fMACR1OH(i)=fMACR1OH(i)-1; fOH(i)=fOH(i)-1; fCH3OO(i)=fCH3OO(i)+0.65; fHCHO(i)=fHCHO(i)+1; fCH3CO3(i)=fCH3CO3(i)+0.35; 

i=i+1;
Rnames{i} = 'HPAC = CH3CO3 + HCHO + OH';
k(:,i) = J22+J41 ; % SUN*3.0E-5;
Gstr{i,1} = 'HPAC'; 
fHPAC(i)=fHPAC(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fHCHO(i)=fHCHO(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'HPAC + OH = MGLY + OH';
k(:,i) = 8.39E-12 ;
Gstr{i,1} = 'HPAC'; Gstr{i,2} = 'OH'; 
fHPAC(i)=fHPAC(i)-1; fMGLY(i)=fMGLY(i)+1; 

i=i+1;
Rnames{i} = 'IDN = 0.1MVK + 0.01MACR + 0.11HCHO + 1.11NO2 + 0.455INO + 0.455ISOP1CO4N + 0.455HO2';
k(:,i) = J53 + J55 ; % SUN*5.0E-5;
Gstr{i,1} = 'IDN'; 
fIDN(i)=fIDN(i)-1; fMVK(i)=fMVK(i)+0.1; fMACR(i)=fMACR(i)+0.01; fHCHO(i)=fHCHO(i)+0.11; fNO2(i)=fNO2(i)+1.11; fINO(i)=fINO(i)+0.455; fISOP1CO4N(i)=fISOP1CO4N(i)+0.455; fHO2(i)=fHO2(i)+0.455; 

i=i+1;
Rnames{i} = 'IDN + OH = NO2 + IHNE';
k(:,i) = F0AM_isop_EPO(T,M,2.37E-11,390,2.715E-19);
Gstr{i,1} = 'IDN'; Gstr{i,2} = 'OH'; 
fIDN(i)=fIDN(i)-1; fOH(i)=fOH(i)-1; fNO2(i)=fNO2(i)+1; fIHNE(i)=fIHNE(i)+1; 

i=i+1;
Rnames{i} = 'IDN + OH = IDNOO';
k(:,i) = 0.74*2.37E-11*exp(390./T) ;
Gstr{i,1} = 'IDN'; Gstr{i,2} = 'OH'; 
fIDN(i)=fIDN(i)-1; fOH(i)=fOH(i)-1; fIDNOO(i)=fIDNOO(i)+1; 

i=i+1;
Rnames{i} = 'IDNOO + NO = PROPNN + 1.11NO2 + 0.11GLYC + 0.89ETHLN + 0.89HO2';
k(:,i) = 2.7E-12*exp(350./T);
Gstr{i,1} = 'IDNOO'; Gstr{i,2} = 'NO'; 
fIDNOO(i)=fIDNOO(i)-1; fNO(i)=fNO(i)-1; fPROPNN(i)=fPROPNN(i)+1; fNO2(i)=fNO2(i)+1.11; fGLYC(i)=fGLYC(i)+0.11; fETHLN(i)=fETHLN(i)+0.89; fHO2(i)=fHO2(i)+0.89; 

i=i+1;
Rnames{i} = 'IDNOO + HO2 = 0.18IHPDN + 0.82PROPNN + 0.82OH + 0.09GLYC + 0.09NO2 + 0.73ETHLN + 0.73HO2';
k(:,i) = 2.71E-13*exp(1300./T); 
Gstr{i,1} = 'IDNOO'; Gstr{i,2} = 'HO2'; 
fIDNOO(i)=fIDNOO(i)-1; fHO2(i)=fHO2(i)-0.27; fIHPDN(i)=fIHPDN(i)+0.18; fPROPNN(i)=fPROPNN(i)+0.82; fOH(i)=fOH(i)+0.82; fGLYC(i)=fGLYC(i)+0.09; fETHLN(i)=fETHLN(i)+0.73; fNO2(i)=fNO2(i)+0.09; 

i=i+1;
Rnames{i} = 'MACRNO2 + HO2 = 0.5HAC + 0.5OH + 0.5NO2 + 0.5CO2 + 0.13O3 + 0.13MACRNOH + 0.37MACRNOOH';
k(:,i) = 3.14E-12*exp(580./T) ;
Gstr{i,1} = 'MACRNO2'; Gstr{i,2} = 'HO2'; 
fMACRNO2(i)=fMACRNO2(i)-1; fHO2(i)=fHO2(i)-1; fHAC(i)=fHAC(i)+0.5; fOH(i)=fOH(i)+0.5; fNO2(i)=fNO2(i)+0.5; fO3(i)=fO3(i)+0.13; fMACRNOH(i)=fMACRNOH(i)+0.13; fMACRNOOH(i)=fMACRNOOH(i)+0.37; 

i=i+1;
Rnames{i} = 'MACRNO2 + NO = HAC + 2NO2 + CO2';
k(:,i) = 7.5E-12*exp(290./T) ; 
Gstr{i,1} = 'MACRNO2'; Gstr{i,2} = 'NO'; 
fMACRNO2(i)=fMACRNO2(i)-1; fNO(i)=fNO(i)-1; fHAC(i)=fHAC(i)+1; fNO2(i)=fNO2(i)+2; 

i=i+1;
Rnames{i} = 'MACRNO2 + NO2 = MPANHN';
k(:,i) = F0AM_isop_TROE(2.591E-28,0,-6.87,1.125E-11,0,-1.105,0.3,T,M);
Gstr{i,1} = 'MACRNO2'; Gstr{i,2} = 'NO2'; 
fMACRNO2(i)=fMACRNO2(i)-1; fNO2(i)=fNO2(i)-1; fMPANHN(i)=fMPANHN(i)+1; 

i=i+1;
Rnames{i} = 'MACRNO2 + NO3 = HAC + 2NO2 + CO2';
k(:,i) = 4.0E-12 ; 
Gstr{i,1} = 'MACRNO2'; Gstr{i,2} = 'NO3'; 
fMACRNO2(i)=fMACRNO2(i)-1; fNO3(i)=fNO3(i)-1; fHAC(i)=fHAC(i)+1; fNO2(i)=fNO2(i)+2; 

i=i+1;
Rnames{i} = 'MACRNO2 + CH3OO = 0.7HAC + 0.7CO2 + 0.7NO2 + 0.7HO2 + HCHO + 0.3MACRNOH';
k(:,i) = 2.9E-12*exp(500./T) ;
Gstr{i,1} = 'MACRNO2'; Gstr{i,2} = 'CH3OO'; 
fMACRNO2(i)=fMACRNO2(i)-1; fCH3OO(i)=fCH3OO(i)-1; fHAC(i)=fHAC(i)+0.7; fNO2(i)=fNO2(i)+0.7; fHO2(i)=fHO2(i)+0.7; fHCHO(i)=fHCHO(i)+1; fMACRNOH(i)=fMACRNOH(i)+0.3; 

i=i+1;
Rnames{i} = 'MACRNOH + OH = HAC + NO2 + CO2';
k(:,i) = 1.34E-12 ;
Gstr{i,1} = 'MACRNOH'; Gstr{i,2} = 'OH'; 
fMACRNOH(i)=fMACRNOH(i)-1; fOH(i)=fOH(i)-1; fHAC(i)=fHAC(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'MACRNOOH + OH = MACRNO2';
k(:,i) = 4.42E-12 ;
Gstr{i,1} = 'MACRNOOH'; Gstr{i,2} = 'OH'; 
fMACRNOOH(i)=fMACRNOOH(i)-1; fOH(i)=fOH(i)-1; fMACRNO2(i)=fMACRNO2(i)+1;

i=i+1;
Rnames{i} = 'MACRNOOH = HAC + NO2 + OH + CO2';
k(:,i) = J22+J41+J51; % SUN*6.49E-6 ;
Gstr{i,1} = 'MACRNOOH'; 
fMACRNOOH(i)=fMACRNOOH(i)-1; fHAC(i)=fHAC(i)+1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'MPANHN = MACRNO2 + NO2';
k(:,i) = 1.58E16*exp(-13500./T);
Gstr{i,1} = 'MPANHN';
fMPANHN(i)=fMPANHN(i)-1; fMACRNO2(i)=fMACRNO2(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ETHLN = NO2 + HCHO + CO + HO2';
k(:,i) = 4.3*J56 ; % SUN*6.46E-5 ;
Gstr{i,1} = 'ETHLN';
fETHLN(i)=fETHLN(i)-1; fNO2(i)=fNO2(i)+1; fHCHO(i)=fHCHO(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'ETHLN + OH = NCH2CO3';
k(:,i) = 3.4E-12 ;
Gstr{i,1} = 'ETHLN'; Gstr{i,2} = 'OH'; 
fETHLN(i)=fETHLN(i)-1; fOH(i)=fOH(i)-1; fNCH2CO3(i)=fNCH2CO3(i)+1; 

i=i+1;
Rnames{i} = 'ETHLN + NO3 = HNO3 + NCH2CO3';
k(:,i) = 1.4E-12*exp(-1860./T);
Gstr{i,1} = 'ETHLN'; Gstr{i,2} = 'NO3'; 
fETHLN(i)=fETHLN(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fNCH2CO3(i)=fNCH2CO3(i)+1; 

i=i+1;
Rnames{i} = 'GLYC = HCHO + 2HO2 + CO';
k(:,i) = J15 ; % SUN*1.66E-5 ;
Gstr{i,1} = 'GLYC'; 
fGLYC(i)=fGLYC(i)-1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+2; fCO(i)=fCO(i)+1; 

i=i+1;
Rnames{i} = 'GLYC + OH = 0.2GLYX + 0.2HO2 + 0.8HOCH2CO3';
k(:,i) = 8.0E-12 ;
Gstr{i,1} = 'GLYC'; Gstr{i,2} = 'OH'; 
fGLYC(i)=fGLYC(i)-1; fOH(i)=fOH(i)-1; fGLYX(i)=fGLYX(i)+0.2; fHO2(i)=fHO2(i)+0.2; fHOCH2CO3(i)=fHOCH2CO3(i)+0.8; 

i=i+1;
Rnames{i} = 'GLYC + NO3 = HNO3 + HOCH2CO3';
k(:,i) = 1.4E-12*exp(-1860./T);
Gstr{i,1} = 'GLYC'; Gstr{i,2} = 'NO3'; 
fGLYC(i)=fGLYC(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fHOCH2CO3(i)=fHOCH2CO3(i)+1; 

i=i+1;
Rnames{i} = 'PYRAC + OH = CH3CO3 + CO2';
k(:,i) = 8.0E-13 ;
Gstr{i,1} = 'PYRAC'; Gstr{i,2} = 'OH'; 
fPYRAC(i)=fPYRAC(i)-1; fOH(i)=fOH(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'PYRAC = CH3CO3 + CO2 + HO2';
k(:,i) = J34 ; % SUN*3.1E-3 ;
Gstr{i,1} = 'PYRAC'; 
fPYRAC(i)=fPYRAC(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'PROPNN + OH = MGLY + NO2';
k(:,i) = 1.0E-13 ;
Gstr{i,1} = 'PROPNN'; Gstr{i,2} = 'OH'; 
fPROPNN(i)=fPROPNN(i)-1; fOH(i)=fOH(i)-1; fMGLY(i)=fMGLY(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'PROPNN = NO2 + HCHO + CH3CO3';
k(:,i) = J56 ; % SUN*8.44E-6;
Gstr{i,1} = 'PROPNN'; 
fPROPNN(i)=fPROPNN(i)-1; fNO2(i)=fNO2(i)+1; fHCHO(i)=fHCHO(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'HCOOH + OH = HO2 + CO2';
k(:,i) = 4.5E-13 ;
Gstr{i,1} = 'HCOOH'; Gstr{i,2} = 'OH'; 
fHCOOH(i)=fHCOOH(i)-1; fOH(i)=fOH(i)-1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HAC + OH = MGLY + HO2';
k(:,i) = 1.6E-12*exp(305./T);
Gstr{i,1} = 'HAC'; Gstr{i,2} = 'OH'; 
fHAC(i)=fHAC(i)-1; fOH(i)=fOH(i)-1; fMGLY(i)=fMGLY(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HAC = CH3CO3 + HO2 + HCHO';
k(:,i) = J22 ; % SUN*2.7E-6;
Gstr{i,1} = 'HAC'; 
fHAC(i)=fHAC(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'MGLY + OH = CH3CO3 + CO';
k(:,i) = 1.9E-12*exp(575./T) ;
Gstr{i,1} = 'MGLY'; Gstr{i,2} = 'OH'; 
fMGLY(i)=fMGLY(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'MGLY + NO3 = CH3CO3 + CO + HNO3';
k(:,i) = 2.4*1.4E-12*exp(-1860./T) ;
Gstr{i,1} = 'MGLY'; Gstr{i,2} = 'NO3'; 
fMGLY(i)=fMGLY(i)-1; fNO3(i)=fNO3(i)-1; fCO(i)=fCO(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; fHNO3(i)=fHNO3(i)+1;

i=i+1;
Rnames{i} = 'MGLY = CH3CO3 + CO + HO2';
k(:,i) = J34 ; % SUN*3.1E-3;
Gstr{i,1} = 'MGLY'; 
fMGLY(i)=fMGLY(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'GLYX = CO + CO + H2';
k(:,i) = J31 ; % SUN*8.61E-6 ; 
Gstr{i,1} = 'GLYX'; 
fGLYX(i)=fGLYX(i)-1; fCO(i)=fCO(i)+2; fH2(i)=fH2(i)+1; 

i=i+1;
Rnames{i} = 'GLYX = 2CO + 2HO2';
k(:,i) = J33 ; % SUN*8.36E-5 ; 
Gstr{i,1} = 'GLYX'; 
fGLYX(i)=fGLYX(i)-1; fCO(i)=fCO(i)+2; fHO2(i)=fHO2(i)+2; 

i=i+1;
Rnames{i} = 'GLYX = HCHO + CO';
k(:,i) = J32 ; % SUN*3.71E-5 ; 
Gstr{i,1} = 'GLYX'; 
fGLYX(i)=fGLYX(i)-1; fCO(i)=fCO(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'GLYX + NO3 = HNO3 + HCOCO';
k(:,i) = 1.4E-12*exp(-1860./T) ;
Gstr{i,1} = 'GLYX'; Gstr{i,2} = 'NO3'; 
fGLYX(i)=fGLYX(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+1; fHCOCO(i)=fHCOCO(i)+1; 

i=i+1;
Rnames{i} = 'GLYX + OH = HCOCO';
k(:,i) = 3.1E-12*exp(340./T);
Gstr{i,1} = 'GLYX'; Gstr{i,2} = 'OH'; 
fGLYX(i)=fGLYX(i)-1; fOH(i)=fOH(i)-1; fHCOCO(i)=fHCOCO(i)+1; 

i=i+1;
Rnames{i} = 'HCOCO = 2CO + HO2';
k(:,i) = 7.0E11*exp(-3160./T) ;
Gstr{i,1} = 'HCOCO'; 
fHCOCO(i)=fHCOCO(i)-1; fCO(i)=fCO(i)+2; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HCOCO + O2 = 0.76CO2 + 0.76OH + 0.24HO2 + 1.24CO';
k(:,i) =  5.0E-12*4.2*M.*0.79 ;
Gstr{i,1} = 'HCOCO';
fHCOCO(i)=fHCOCO(i)-1; fCO(i)=fCO(i)+1.24; fHO2(i)=fHO2(i)+0.24; fOH(i)=fOH(i)+0.76; 

i=i+1;
Rnames{i} = 'CH3CO3 + HO2 = 0.13O3 + 0.13CH3CO2H + 0.37CH3CO3H + 0.5CH3OO + 0.5CO2 + 0.5OH';
k(:,i) = 3.14E-12*exp(580./T) ;
Gstr{i,1} = 'CH3CO3'; Gstr{i,2} = 'HO2'; 
fCH3CO3(i)=fCH3CO3(i)-1; fHO2(i)=fHO2(i)-1; fO3(i)=fO3(i)+0.13; fCH3CO2H(i)=fCH3CO2H(i)+0.13; fCH3CO3H(i)=fCH3CO3H(i)+0.37; fCH3OO(i)=fCH3OO(i)+0.5; fOH(i)=fOH(i)+0.5; 

i=i+1;
Rnames{i} = 'CH3CO3 + NO = NO2 + CH3OO + CO2';
k(:,i) = 7.5E-12*exp(290./T) ; 
Gstr{i,1} = 'CH3CO3'; Gstr{i,2} = 'NO'; 
fCH3CO3(i)=fCH3CO3(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fCH3OO(i)=fCH3OO(i)+1; 

i=i+1;
Rnames{i} = 'CH3CO3 + NO2 = PAN';
k(:,i) = F0AM_isop_TROE(2.591E-28,0,-6.87,1.125E-11,0,-1.105,0.3,T,M);
Gstr{i,1} = 'CH3CO3'; Gstr{i,2} = 'NO2'; 
fCH3CO3(i)=fCH3CO3(i)-1; fNO2(i)=fNO2(i)-1; fPAN(i)=fPAN(i)+1; 

i=i+1;
Rnames{i} = 'PAN = CH3CO3 + NO2';
k(:,i) = F0AM_isop_TROE(3.871E-3,-12100,0,5.4E16,-13830,0,0.3,T,M);
Gstr{i,1} = 'PAN'; 
fPAN(i)=fPAN(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'PAN + OH = HCHO + CO + NO2';
k(:,i) = 3.0E-14 ;
Gstr{i,1} = 'PAN'; Gstr{i,2} = 'OH'; 
fPAN(i)=fPAN(i)-1; fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)+1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3CO3 + NO3 = NO2 + CH3OO + CO2';
k(:,i) = 4.0E-12 ; 
Gstr{i,1} = 'CH3CO3'; Gstr{i,2} = 'NO3'; 
fCH3CO3(i)=fCH3CO3(i)-1; fNO3(i)=fNO3(i)-1; fCH3OO(i)=fCH3OO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'CH3CO3 + MO2 = 0.7HO2 + HCHO + 0.7CO2 + 0.7CH3OO + 0.3CH3CO2H';
k(:,i) = 2.9E-12*exp(500./T) ;
Gstr{i,1} = 'CH3CO3'; Gstr{i,2} = 'CH3OO'; 
fCH3CO3(i)=fCH3CO3(i)-1; fCH3OO(i)=fCH3OO(i)-0.3; fHO2(i)=fHO2(i)+0.7; fHCHO(i)=fHCHO(i)+1; fCH3CO2H(i)=fCH3CO2H(i)+0.3; 

i=i+1;
Rnames{i} = 'CH3CO2H + OH = CH3OO + CO2';
k(:,i) = 8.0E-13 ;
Gstr{i,1} = 'CH3CO2H'; Gstr{i,2} = 'OH'; 
fCH3CO2H(i)=fCH3CO2H(i)-1; fOH(i)=fOH(i)-1; fCH3OO(i)=fCH3OO(i)+1; 

i=i+1;
Rnames{i} = 'CH3CO3H + OH = CH3CO3';
k(:,i) = 3.7E-12 ;
Gstr{i,1} = 'CH3CO3H'; Gstr{i,2} = 'OH'; 
fCH3CO3H(i)=fCH3CO3H(i)-1; fOH(i)=fOH(i)-1; fCH3CO3(i)=fCH3CO3(i)+1; 

i=i+1;
Rnames{i} = 'CH3CO3H = CH3OO + OH + CO2';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'CH3CO3H'; 
fCH3CO3H(i)=fCH3CO3H(i)-1; fCH3OO(i)=fCH3OO(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'NCH2CO3 + HO2 = 0.5HCHO + 0.5NO2 + 0.5OH + 0.5CO2 + 0.13O3 + 0.13NCH2CO2H + 0.37NCH2CO3H';
k(:,i) = 3.14E-12*exp(580./T);
Gstr{i,1} = 'NCH2CO3'; Gstr{i,2} = 'HO2'; 
fNCH2CO3(i)=fNCH2CO3(i)-1; fHO2(i)=fHO2(i)-1; fHCHO(i)=fHCHO(i)+0.5; fNO2(i)=fNO2(i)+0.5; fOH(i)=fOH(i)+0.5; fO3(i)=fO3(i)+0.5; fNCH2CO2H(i)=fNCH2CO2H(i)+0.13; fNCH2CO3H(i)=fNCH2CO3H(i)+0.37;  

i=i+1;
Rnames{i} = 'NCH2CO3 + NO = CO2 + HCHO + 2NO2';
k(:,i) = 7.5E-12*exp(290./T) ;
Gstr{i,1} = 'NCH2CO3'; Gstr{i,2} = 'NO'; 
fNCH2CO3(i)=fNCH2CO3(i)-1; fNO(i)=fNO(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+2; 

i=i+1;
Rnames{i} = 'NCH2CO3 + NO2 = PNAN';
k(:,i) = F0AM_isop_TROE(2.591E-28,0,-6.87,1.125E-11,0,-1.105,0.3,T,M);
Gstr{i,1} = 'NCH2CO3'; Gstr{i,2} = 'NO2'; 
fNCH2CO3(i)=fNCH2CO3(i)-1; fNO2(i)=fNO2(i)-1; fPNAN(i)=fPNAN(i)+1;

i=i+1;
Rnames{i} = 'NCH2CO3 + NO3 = 2NO2 + HCHO + CO2';
k(:,i) = 2.3E-12*1.74 ;
Gstr{i,1} = 'NCH2CO3'; Gstr{i,2} = 'NO3'; 
fNCH2CO3(i)=fNCH2CO3(i)-1; fNO3(i)=fNO3(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+2; 

i=i+1;
Rnames{i} = 'NCH2CO3 + CH3OO = 1.7HCHO + 0.7NO2 + 0.7CO2 + 0.7HO2 + 0.3NCH2CO2H';
k(:,i) = 1.0E-11 ;
Gstr{i,1} = 'NCH2CO3'; Gstr{i,2} = 'CH3OO'; 
fNCH2CO3(i)=fNCH2CO3(i)-1; fCH3OO(i)=fCH3OO(i)-1; fHCHO(i)=fHCHO(i)+1.7; fNO2(i)=fNO2(i)+0.7; fHO2(i)=fHO2(i)+0.7; fNCH2CO2H(i)=fNCH2CO2H(i)+0.3; 

i=i+1;
Rnames{i} = 'HOCH2CO3 + HO2 = 0.5HO2 + 0.5CO2 + 0.5HCHO + 0.5OH + 0.13O3 + 0.13HOCH2CO2H + 0.37HOCH2CO3H';
k(:,i) = 3.14E-12*exp(580./T) ;
Gstr{i,1} = 'HOCH2CO3'; Gstr{i,2} = 'HO2'; 
fHOCH2CO3(i)=fHOCH2CO3(i)-1; fHO2(i)=fHO2(i)-0.5; fHCHO(i)=fHCHO(i)+0.5; fOH(i)=fOH(i)+0.5; fO3(i)=fO3(i)+0.5; fHOCH2CO2H(i)=fHOCH2CO2H(i)+0.13; fHOCH2CO3H(i)=fHOCH2CO3H(i)+0.37; 

i=i+1;
Rnames{i} = 'HOCH2CO3 + NO = NO2 + HO2 + HCHO + CO2';
k(:,i) = 7.5E-12*exp(290./T) ;
Gstr{i,1} = 'HOCH2CO3'; Gstr{i,2} = 'NO'; 
fHOCH2CO3(i)=fHOCH2CO3(i)-1; fNO(i)=fNO(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'HOCH2CO3 + NO2 = PHAN';
k(:,i) = F0AM_isop_TROE(2.591E-28,0,-6.87,1.125E-11,0,-1.105,0.3,T,M);
Gstr{i,1} = 'HOCH2CO3'; Gstr{i,2} = 'NO2'; 
fHOCH2CO3(i)=fHOCH2CO3(i)-1; fNO2(i)=fNO2(i)-1; fPHAN(i)=fPHAN(i)+1; 

i=i+1;
Rnames{i} = 'HOCH2CO3 + NO3 = NO2 + HO2 + HCHO + CO2';
k(:,i) = 2.3E-12*1.74 ;
Gstr{i,1} = 'HOCH2CO3'; Gstr{i,2} = 'NO3'; 
fHOCH2CO3(i)=fHOCH2CO3(i)-1; fNO3(i)=fNO3(i)-1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; 

i=i+1;
Rnames{i} = 'HOCH2CO3 + CH3OO = 1.7HCHO + 0.7HO2 + 0.7CO2 + 0.3HOCH2CO2H';
k(:,i) = 1.0E-11 ;
Gstr{i,1} = 'HOCH2CO3'; Gstr{i,2} = 'CH3OO'; 
fHOCH2CO3(i)=fHOCH2CO3(i)-1; fCH3OO(i)=fCH3OO(i)-1; fHCHO(i)=fHCHO(i)+1.7; fHO2(i)=fHO2(i)+0.7; fHOCH2CO2H(i)=fHOCH2CO2H(i)+0.3; 

i=i+1;
Rnames{i} = 'HOCH2CO2H + OH = HCHO + CO2 + HO2';
k(:,i) = 2.73E-12 ;
Gstr{i,1} = 'HOCH2CO2H'; Gstr{i,2} = 'OH'; 
fHOCH2CO2H(i)=fHOCH2CO2H(i)-1; fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'HOCH2CO3H + OH = HOCH2CO3';
k(:,i) = 6.19E-12 ;
Gstr{i,1} = 'HOCH2CO3H'; Gstr{i,2} = 'OH'; 
fHOCH2CO3H(i)=fHOCH2CO3H(i)-1; fOH(i)=fOH(i)-1; fHOCH2CO3(i)=fHOCH2CO3(i)+1; 

i=i+1;
Rnames{i} = 'HOCH2CO3H = HCHO + CO2 + HO2 + OH';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'HOCH2CO3H'; 
fHOCH2CO3H(i)=fHOCH2CO3H(i)-1; fHCHO(i)=fHCHO(i)+1; fHO2(i)=fHO2(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'PHAN = NO2 + HOCH2CO3';
k(:,i) = F0AM_isop_TROE(3.871E-3,-12100,0,5.4E16,-13830,0,0.3,T,M);
Gstr{i,1} = 'PHAN'; 
fPHAN(i)=fPHAN(i)-1; fNO2(i)=fNO2(i)+1; fHOCH2CO3(i)=fHOCH2CO3(i)+1; 

i=i+1;
Rnames{i} = 'PHAN + OH = HCHO + CO + NO2';
k(:,i) = 1.12E-12 ;
Gstr{i,1} = 'PHAN'; Gstr{i,2} = 'OH'; 
fPHAN(i)=fPHAN(i)-1; fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)+1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NCH2CO2H + OH = HCHO + NO2 + CO2';
k(:,i) = 1.68E-13 ;
Gstr{i,1} = 'NCH2CO2H'; Gstr{i,2} = 'OH'; 
fNCH2CO2H(i)=fNCH2CO2H(i)-1; fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'NCH2CO3H + OH = NCH2CO3';
k(:,i) = 3.63E-12 ;
Gstr{i,1} = 'NCH2CO3H'; Gstr{i,2} = 'OH'; 
fNCH2CO3H(i)=fNCH2CO3H(i)-1; fOH(i)=fOH(i)-1; fNCH2CO3(i)=fNCH2CO3(i)+1;

i=i+1;
Rnames{i} = 'NCH2CO3H = HCHO + NO2 + OH + CO2';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'NCH2CO3H';
fNCH2CO3H(i)=fNCH2CO3H(i)-1; fHCHO(i)=fHCHO(i)+1; fNO2(i)=fNO2(i)+1; fOH(i)=fOH(i)+1; 

i=i+1;
Rnames{i} = 'PNAN = NCH2CO3 + NO2';
k(:,i) = F0AM_isop_TROE(3.871E-3,-12100,0,5.4E16,-13830,0,0.3,T,M);
Gstr{i,1} = 'PNAN'; 
fPNAN(i)=fPNAN(i)-1; fNCH2CO3(i)=fNCH2CO3(i)+1; fNO2(i)=fNO2(i)+1;

i=i+1;
Rnames{i} = 'PNAN + OH = HCHO + CO + 2NO2';
k(:,i) = 1.12E-14 ;
Gstr{i,1} = 'PNAN'; Gstr{i,2} = 'OH'; 
fPNAN(i)=fPNAN(i)-1; fOH(i)=fOH(i)-1; fHCHO(i)=fHCHO(i)+1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+2; 

i=i+1;
Rnames{i} = 'MVK + O3 = 0.545MGLY + 0.5SCI + 0.6HCHO + 0.38CH3CO3 + 0.1HO2 + 0.08OH + 0.18CO + 0.075PYRAC + 0.045H2O2';
k(:,i) = 8.5E-16*exp(-1520./T) ;
Gstr{i,1} = 'MVK'; Gstr{i,2} = 'O3'; 
fMVK(i)=fMVK(i)-1; fO3(i)=fO3(i)-1; fMGLY(i)=fMGLY(i)+0.545; fSCI(i)=fSCI(i)+0.5; fHCHO(i)=fHCHO(i)+0.6; fCH3CO3(i)=fCH3CO3(i)+0.38; fHO2(i)=fHO2(i)+0.1; fOH(i)=fOH(i)+0.08; fCO(i)=fCO(i)+0.18; fPYRAC(i)=fPYRAC(i)+0.075; fH2O2(i)=fH2O2(i)+0.045; 

i=i+1;
Rnames{i} = 'MACR + O3 = 0.88MGLY + 0.88SCI + 0.12HCHO + 0.12OH + 0.12CO + 0.12CH3CO3';
k(:,i) = 1.4E-15*exp(-2100./T) ;
Gstr{i,1} = 'MACR'; Gstr{i,2} = 'O3'; 
fMACR(i)=fMACR(i)-1; fO3(i)=fO3(i)-1; fMGLY(i)=fMGLY(i)+0.88; fSCI(i)=fSCI(i)+0.88; fHCHO(i)=fHCHO(i)+0.12; fOH(i)=fOH(i)+0.12; fCO(i)=fCO(i)+0.12; fCH3CO3(i)=fCH3CO3(i)+0.12; 

i=i+1;
Rnames{i} = 'MACR + NO3 = 0.32HNO3 + 0.32MACR1OO + 0.68OH + 0.68CO + 0.68PROPNN';
k(:,i) = 1.8E-13*exp(-1190./T) ;
Gstr{i,1} = 'MACR'; Gstr{i,2} = 'NO3'; 
fMACR(i)=fMACR(i)-1; fNO3(i)=fNO3(i)-1; fHNO3(i)=fHNO3(i)+0.32; fMACR1OO(i)=fMACR1OO(i)+0.32; fOH(i)=fOH(i)+0.68; fCO(i)=fCO(i)+0.68; fPROPNN(i)=fPROPNN(i)+0.68; 

i=i+1;
Rnames{i} = 'MACR = 1.65CO + HO2 + HCHO + 0.65CH3OO + 0.35CH3CO3';
k(:,i) = J18;
Gstr{i,1} = 'MACR'; 
fMACR(i)=fMACR(i)-1; fCO(i)=fCO(i)+1.65; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fCH3OO(i)=fCH3OO(i)+0.65; fCH3CO3(i)=fCH3CO3(i)+0.35; 

i=i+1;
Rnames{i} = 'MACR = HO2 + MACR1OO';
k(:,i) = J19;
Gstr{i,1} = 'MACR'; 
fMACR(i)=fMACR(i)-1; fHO2(i)=fHO2(i)+1; fMACR1OO(i)=fMACR1OO(i)+1; 

i=i+1;
Rnames{i} = 'MVK = HO2 + HCHO + CH3CO3 + CO';
k(:,i) = J23+J24;
Gstr{i,1} = 'MVK'; 
fMVK(i)=fMVK(i)-1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+1; fCH3CO3(i)=fCH3CO3(i)+1; fCO(i)=fCO(i)+1; 

i=i+1;
Rnames{i} = 'ICHNP + OH = CO + NO2 + 0.75MVK3OOH4OH + 0.25MACR2OOH3OH';
k(:,i) = 1.0E-11 ;
Gstr{i,1} = 'ICHNP'; Gstr{i,2} = 'OH'; 
fICHNP(i)=fICHNP(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+1; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+1; 

i=i+1;
Rnames{i} = 'ICHNP = MGLY + OH + NO2 + GLYC';
k(:,i) = J41 + J54 ; % SUN*(1.15E-5) ;
Gstr{i,1} = 'ICHNP'; 
fICHNP(i)=fICHNP(i)-1; fMGLY(i)=fMGLY(i)+1; fGLYC(i)=fGLYC(i)+1; fOH(i)=fOH(i)+1; fNO2(i)=fNO2(i)+1; 

i=i+1;
Rnames{i} = 'ICHNP = 0.5MVK3OOH4OH + 0.5MACR2OOH3OH + CO + NO2 + HO2';
k(:,i) = J22+J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'ICHNP'; 
fICHNP(i)=fICHNP(i)-1; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.5; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.5; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; fHO2(i)=fHO2(i)+1; 

i=i+1;
Rnames{i} = 'IHNDP + OH = OH + ICHNP';
k(:,i) = 2.0E-12 ;
Gstr{i,1} = 'IHNDP'; Gstr{i,2} = 'OH'; 
fIHNDP(i)=fIHNDP(i)-1; fICHNP(i)=fICHNP(i)+1; 

i=i+1;
Rnames{i} = 'IHNDP = NO2 + IHPOO1';
k(:,i) = J54 ; % SUN*5.0E-6 ;
Gstr{i,1} = 'IHNDP'; 
fIHNDP(i)=fIHNDP(i)-1; fNO2(i)=fNO2(i)+1; fIHPOO1(i)=fIHPOO1(i)+1; 

i=i+1;
Rnames{i} = 'IHNDP = OH + 0.5IDHNDOO1 + 0.5IDHNBOO';
k(:,i) = 2*J41 ; % SUN*1.3E-5 ;
Gstr{i,1} = 'IHNDP';
fIHNDP(i)=fIHNDP(i)-1; fOH(i)=fOH(i)+1; fIDHNDOO1(i)=fIDHNDOO1(i)+0.5; fIDHNBOO(i)=fIDHNBOO(i)+0.5; 

i=i+1;
Rnames{i} = 'IDHDN + OH = ICHDN';
k(:,i) = 2.0E-12 ;
Gstr{i,1} = 'IDHDN'; Gstr{i,2} = 'OH'; 
fIDHDN(i)=fIDHDN(i)-1; fOH(i)=fOH(i)-1; fICHDN(i)=fICHDN(i)+1; 

i=i+1;
Rnames{i} = 'IDHDN = 2NO2 + GLYC + HAC';
k(:,i) = J53 + J55 ; % SUN*5.0E-6*2 ;
Gstr{i,1} = 'IDHDN'; 
fIDHDN(i)=fIDHDN(i)-1; fNO2(i)=fNO2(i)+2; fGLYC(i)=fGLYC(i)+1; fHAC(i)=fHAC(i)+1; 

i=i+1;
Rnames{i} = 'IDHPN + OH = 0.33OH + 0.67HO2 + 0.33IDHCN + 0.67ICHNP';
k(:,i) = 3.0E-12 ;
Gstr{i,1} = 'IDHPN'; Gstr{i,2} = 'OH'; 
fIDHPN(i)=fIDHPN(i)-1; fOH(i)=fOH(i)-0.67; fHO2(i)=fHO2(i)+0.67; fIDHCN(i)=fIDHCN(i)+0.33; fICHNP(i)=fICHNP(i)+0.67; 

i=i+1;
Rnames{i} = 'IDHPN = OH + 0.7HO2 + 0.55HCHO + 0.5MACR2N3OH + 0.3GLYC + 0.45HAC + 0.3NO2 + 0.15ETHLN + 0.05MVK3N4OH';
k(:,i) = J41+J51 ; % SUN*(6.5E-6) ;
Gstr{i,1} = 'IDHPN'; 
fIDHPN(i)=fIDHPN(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+0.7; fHCHO(i)=fHCHO(i)+0.55; fMACR2N3OH(i)=fMACR2N3OH(i)+0.5; fGLYC(i)=fGLYC(i)+0.3; fHAC(i)=fHAC(i)+0.45; fNO2(i)=fNO2(i)+0.3; fETHLN(i)=fETHLN(i)+0.15; fMVK3N4OH(i)=fMVK3N4OH(i)+0.05; 

i=i+1;
Rnames{i} = 'IDHPN = NO2 + 0.8HAC + 0.7HO2 + 0.5HPETHNL + 0.3OH + 0.35GLYC + 0.15HCHO + 0.15MACR2OOH3OH + 0.05HPAC';
k(:,i) = J41+J51 ; % SUN*(6.5E-6) ;
Gstr{i,1} = 'IDHPN'; 
fIDHPN(i)=fIDHPN(i)-1; fNO2(i)=fNO2(i)+1; fHAC(i)=fHAC(i)+0.8; fHO2(i)=fHO2(i)+0.7; fHPETHNL(i)=fHPETHNL(i)+0.5; fOH(i)=fOH(i)+0.3; fGLYC(i)=fGLYC(i)+0.35; fHCHO(i)=fHCHO(i)+0.15; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.15; fHPAC(i)=fHPAC(i)+0.05; 

i=i+1;
Rnames{i} = 'IDHCN + OH = CO + NO2 + 0.75MVK3OH4OH + 0.25MACR2OH3OH';
k(:,i) = 1.0E-11 ;
Gstr{i,1} = 'IDHCN'; Gstr{i,2} = 'OH'; 
fIDHCN(i)=fIDHCN(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.75; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.25; 

i=i+1;
Rnames{i} = 'IDHCN = 0.5HO2 + NO2 + 0.5CO + 0.375MVK3OH4OH + 0.125MACR2OH3OH + 0.5HOCH2CO3 + 0.5HAC';
k(:,i) = J22 + J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'IDHCN'; 
fIDHCN(i)=fIDHCN(i)-1; fHO2(i)=fHO2(i)+0.5; fNO2(i)=fNO2(i)+1; fCO(i)=fCO(i)+0.5; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.375; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.125; fHOCH2CO3(i)=fHOCH2CO3(i)+0.5; fHAC(i)=fHAC(i)+0.5; 

i=i+1;
Rnames{i} = 'IDHCN = NO2 + 0.5HO2 + 0.375GLYC + 0.375MGLY + 0.125GLYX + 0.625HAC + 0.5HOCH2CO3';
k(:,i) = J51 ; % SUN*5.0E-6 ;
Gstr{i,1} = 'IDHCN'; 
fIDHCN(i)=fIDHCN(i)-1; fHO2(i)=fHO2(i)+0.5; fNO2(i)=fNO2(i)+1; fGLYC(i)=fGLYC(i)+0.375; fMGLY(i)=fMGLY(i)+0.375; fGLYX(i)=fGLYX(i)+0.125; fHAC(i)=fHAC(i)+0.625; fHOCH2CO3(i)=fHOCH2CO3(i)+0.5; 

i=i+1;
Rnames{i} = 'ICPDH + OH = CO + 0.5HO2 + 0.5OH + 0.5MACR2OOH3OH + 0.35MVK3OH4OH + 0.15MACR2OH3OH';
k(:,i) = 1.0E-11 ;
Gstr{i,1} = 'ICPDH'; Gstr{i,2} = 'OH'; 
fICPDH(i)=fICPDH(i)-1; fOH(i)=fOH(i)-0.5; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+0.5; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.5; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.35; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.15; 

i=i+1;
Rnames{i} = 'ICPDH = CO + 1.5HO2 + 0.5OH + 0.5MACR2OOH3OH + 0.35MVK3OH4OH + 0.15MACR2OH3OH';
k(:,i) = J22 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'ICPDH'; 
fICPDH(i)=fICPDH(i)-1; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1.5; fOH(i)=fOH(i)+0.5; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.5; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.35; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.15; 

i=i+1;
Rnames{i} = 'ICPDH = OH + HO2 + 0.1HCHO + 0.1MVK3OH4CO + 0.438HAC + 0.438GLYX + 0.088GLYC + 0.088MGLY + 0.122CO + 0.122MACR2OH3OH';
k(:,i) = J41 ;
Gstr{i,1} = 'ICPDH'; 
fICPDH(i)=fICPDH(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fHCHO(i)=fHCHO(i)+0.1; fMVK3OH4CO(i)=fMVK3OH4CO(i)+0.1; fHAC(i)=fHAC(i)+0.438; fGLYX(i)=fGLYX(i)+0.438; fGLYC(i)=fGLYC(i)+0.088; fMGLY(i)=fMGLY(i)+0.088; fCO(i)=fCO(i)+0.122; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.122; 

i=i+1;
Rnames{i} = 'IDHDP + OH = OH + 0.333ICPDH + 0.667IDHPE';
k(:,i) = 3.0E-12 ;
Gstr{i,1} = 'IDHDP'; Gstr{i,2} = 'OH'; 
fIDHDP(i)=fIDHDP(i)-1; fICPDH(i)=fICPDH(i)+0.333; fIDHPE(i)=fIDHPE(i)+0.667; 

i=i+1;
Rnames{i} = 'IDHDP = 1.25OH + 0.25GLYC + 0.25HAC + 0.75ICPDH + 0.75HO2';
k(:,i) = 2*J41 ;
Gstr{i,1} = 'IDHDP';
fIDHDP(i)=fIDHDP(i)-1; fOH(i)=fOH(i)+1.25; fGLYC(i)=fGLYC(i)+0.25; fHAC(i)=fHAC(i)+0.25; fICPDH(i)=fICPDH(i)+0.75; fHO2(i)=fHO2(i)+0.75; 

i=i+1;
Rnames{i} = 'IDHPE + OH = OH + CO2 + 0.571MACR2OOH3OH + 0.429MVK3OOH4OH';
k(:,i) = 3.0E-12 ;
Gstr{i,1} = 'IDHPE'; Gstr{i,2} = 'OH'; 
fIDHPE(i)=fIDHPE(i)-1; fMACR2OOH3OH(i)=fMACR2OOH3OH(i)+0.571; fMVK3OOH4OH(i)=fMVK3OOH4OH(i)+0.429; 

i=i+1;
Rnames{i} = 'IDHPE = OH + HO2 + 0.429MGLY + 0.429GLYC + 0.571GLYX + 0.571HAC';
k(:,i) = J41 ;
Gstr{i,1} = 'IDHPE'; 
fIDHPE(i)=fIDHPE(i)-1; fOH(i)=fOH(i)+1; fHO2(i)=fHO2(i)+1; fMGLY(i)=fMGLY(i)+0.429; fGLYC(i)=fGLYC(i)+0.429; fGLYX(i)=fGLYX(i)+0.571; fHAC(i)=fHAC(i)+0.571; 

i=i+1;
Rnames{i} = 'IHNPE + OH = 0.5HO2 + CO + OH + 0.5PROPNN + 0.5CO2 + 0.5MACR2OOH3N';
k(:,i) = 3.0E-12*exp(20./T);
Gstr{i,1} = 'IHNPE'; Gstr{i,2} = 'OH'; 
fIHNPE(i)=fIHNPE(i)-1; fHO2(i)=fHO2(i)+0.5; fCO(i)=fCO(i)+1; fPROPNN(i)=fPROPNN(i)+0.5; fMACR2OOH3N(i)=fMACR2OOH3N(i)+0.5; 

i=i+1;
Rnames{i} = 'IHNPE = OH + PROPNN + HO2 + GLYX';
k(:,i) = J41 + J51 ; 
Gstr{i,1} = 'IHNPE'; 
fIHNPE(i)=fIHNPE(i)-1; fOH(i)=fOH(i)+1; fPROPNN(i)=fPROPNN(i)+1; fHO2(i)=fHO2(i)+1; fGLYX(i)=fGLYX(i)+1; 

i=i+1;
Rnames{i} = 'IHPDN + OH = 0.875IDNOO + 0.125OH + 0.125ICHDN';
k(:,i) = 1.0E-12 ;
Gstr{i,1} = 'IHPDN'; Gstr{i,2} = 'OH'; 
fIHPDN(i)=fIHPDN(i)-1; fOH(i)=fOH(i)-0.875; fIDNOO(i)=fIDNOO(i)+1; fICHDN(i)=fICHDN(i)+1; 

i=i+1;
Rnames{i} = 'IHPDN = OH + 0.82PROPNN + 0.68HO2 + 0.32GLYC + 0.32NO2 + 0.18ICHDN + 0.5ETHLN';
k(:,i) = J41 ; 
Gstr{i,1} = 'IHPDN'; 
fIHPDN(i)=fIHPDN(i)-1; fOH(i)=fOH(i)+1; fPROPNN(i)=fPROPNN(i)+0.82; fHO2(i)=fHO2(i)+0.68; fGLYC(i)=fGLYC(i)+0.32; fNO2(i)=fNO2(i)+0.32; fICHDN(i)=fICHDN(i)+0.18; fETHLN(i)=fETHLN(i)+0.5; 

i=i+1;
Rnames{i} = 'IHPDN = NO2 + 0.25IDHNDOO1 + 0.125IDHNBOO + 0.25PROPNN + 0.25MVK3OH4N + 0.125MVK3N4OH + 0.375HCHO + 0.535OH + 0.09HO2 + 0.09HPETHNL + 0.16GLYC';
k(:,i) = 2.0*J51 ;
Gstr{i,1} = 'IHPDN'; 
fIHPDN(i)=fIHPDN(i)-1; fNO2(i)=fNO2(i)+1; fIDHNDOO1(i)=fIDHNDOO1(i)+0.25; fIDHNBOO(i)=fIDHNBOO(i)+0.125; fPROPNN(i)=fPROPNN(i)+0.25; fMVK3OH4N(i)=fMVK3OH4N(i)+0.25; fMVK3N4OH(i)=fMVK3N4OH(i)+0.125; fHCHO(i)=fHCHO(i)+0.375; fOH(i)=fOH(i)+0.535; fHO2(i)=fHO2(i)+0.09; fHPETHNL(i)=fHPETHNL(i)+0.09; fGLYC(i)=fGLYC(i)+0.16; 

i=i+1;
Rnames{i} = 'ICHDN + OH = CO + NO2 + 0.5MACR2OH3N + 0.5MVK3OH4N';
k(:,i) = 1.0E-11*0.4 ;
Gstr{i,1} = 'ICHDN'; Gstr{i,2} = 'OH'; 
fICHDN(i)=fICHDN(i)-1; fOH(i)=fOH(i)-1; fCO(i)=fCO(i)+1; fNO2(i)=fNO2(i)+1; fMACR2OH3N(i)=fMACR2OH3N(i)+0.5; fMVK3OH4N(i)=fMVK3OH4N(i)+0.5; 

i=i+1;
Rnames{i} = 'ICHDN = 1.5NO2 + 0.65CO + 0.1MACR2OH3OH + 0.1MVK3OH4OH + 0.075MVK3OH4N + 0.075MACR2OH3N + 0.275HO2 + 0.025GLYX + 0.025PROPNN + 0.025MGLY + 0.025ETHLN + 0.525HAC + 0.375HCHO + 0.075MVK3CO4N + 0.225NCH2CO3';
k(:,i) = 2.0*J51 ;
Gstr{i,1} = 'ICHDN'; 
fICHDN(i)=fICHDN(i)-1; fNO2(i)=fNO2(i)+1.5; fCO(i)=fCO(i)+0.65; fMACR2OH3OH(i)=fMACR2OH3OH(i)+0.1; fMVK3OH4OH(i)=fMVK3OH4OH(i)+0.1; fMVK3OH4N(i)=fMVK3OH4N(i)+0.075; fMACR2OH3N(i)=fMACR2OH3N(i)+0.075; fHO2(i)=fHO2(i)+0.275; fGLYX(i)=fGLYX(i)+0.025; fPROPNN(i)=fPROPNN(i)+0.025; fMGLY(i)=fMGLY(i)+0.025; fETHLN(i)=fETHLN(i)+0.025; fHAC(i)=fHAC(i)+0.525; fHCHO(i)=fHCHO(i)+0.375; fMVK3CO4N(i)=fMVK3CO4N(i)+0.075; fNCH2CO3(i)=fNCH2CO3(i)+0.225; 

i=i+1;
Rnames{i} = 'ICHDN = NO2 + 0.6HAC + 0.6NCH2CO3 + 0.4HO2 + 0.4CO + 0.2MACR2OH3N + 0.2MVK3OH4N';
k(:,i) = J22 + J41 ; % SUN*3.0E-5 ;
Gstr{i,1} = 'ICHDN'; 
fICHDN(i)=fICHDN(i)-1; fNO2(i)=fNO2(i)+1; fHAC(i)=fHAC(i)+0.6; fNCH2CO3(i)=fNCH2CO3(i)+0.6; fHO2(i)=fHO2(i)+0.4; fCO(i)=fCO(i)+0.4; fMACR2OH3N(i)=fMACR2OH3N(i)+0.2; fMVK3OH4N(i)=fMVK3OH4N(i)+0.2; 

i=i+1;
Rnames{i} = 'IDCHP + OH = 0.888CO + 0.444OH + 0.444HO2 + 0.318MVK3OH4CO + 0.08IEPOXAOO + 0.126MACR2OH3CO + 0.444MVK3OOH4CO + 0.032IEPOXBOO';
k(:,i) = 2.25E-11 ;
Gstr{i,1} = 'IDCHP'; Gstr{i,2} = 'OH'; 
fIDCHP(i)=fIDCHP(i)-1; fOH(i)=fOH(i)-0.556; fCO(i)=fCO(i)+0.888; fHO2(i)=fHO2(i)+0.444; fMVK3OH4CO(i)=fMVK3OH4CO(i)+0.318; fIEPOXAOO(i)=fIEPOXAOO(i)+0.081; fMACR2OH3CO(i)=fMACR2OH3CO(i)+0.126; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.444; fIEPOXBOO(i)=fIEPOXBOO(i)+0.032; 

i=i+1;
Rnames{i} = 'IDCHP = 0.546OH + CO + 1.454HO2 + 0.391MVK3OH4CO + 0.155MACR2OH3CO + 0.454MVK3OOH4CO';
k(:,i) = 2*J22 + J41 ;
Gstr{i,1} = 'IDCHP'; 
fIDCHP(i)=fIDCHP(i)-1; fOH(i)=fOH(i)+0.546; fCO(i)=fCO(i)+1; fHO2(i)=fHO2(i)+1.454; fMVK3OH4CO(i)=fMVK3OH4CO(i)+0.391; fMACR2OH3CO(i)=fMACR2OH3CO(i)+0.155; fMVK3OOH4CO(i)=fMVK3OOH4CO(i)+0.454; 


%END OF REACTION LIST