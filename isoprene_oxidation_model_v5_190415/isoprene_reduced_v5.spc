#include atoms

#DEFVAR
H2 = IGNORE ;
CO = IGNORE ;
CO2 = IGNORE ;
HCHO = IGNORE ;
HCOOH = IGNORE ;
HAC = IGNORE ; {hydroxyacetone}
MGLY = IGNORE ; {methylglyoxal}
HPETHNL = IGNORE ; {hydroperoxy ethanal}
CH3CO3 = IGNORE ; {acetylperoxy radical}
GLYC = IGNORE ; {glycolaldehyde}
GLYX = IGNORE ; {glyoxal}
HPAC = IGNORE ; {hydroperoxyacetone}
PYRAC = IGNORE ; {pyruvic acid}
SCI = IGNORE ; {C1 stabilized Criegee intermediate}
H2O2 = IGNORE ; {hydrogen peroxide}
HMHP = IGNORE ; {hydroxymethyl hydroperoxide}
CH3OH = IGNORE ; {methanol}
PROPNN = IGNORE ; {propanone nitrate}
ETHLN = IGNORE ; {ethanal nitrate}
CH3OO = IGNORE ;
HMML = IGNORE ;
ISOP = IGNORE ; {isoprene}
IHOO1 = IGNORE ; {isoprene-(1-OH)-OO adducts}
IHOO4 = IGNORE ; {isoprene-(4-OH)-OO adducts}
ISOP1CO4OOH = IGNORE ;
ISOP1OOH4CO = IGNORE ;
ISOP1CO4OH = IGNORE ;
ISOP1OH4CO = IGNORE ;
ISOP3CO4OH = IGNORE ;
ISOP1CO2OOH = IGNORE ;
ISOP3OOH4CO = IGNORE ;
ICHOO = IGNORE ; {ISOP1OH2OO3CO4OH}
ISOP1OH2N = IGNORE ;
ISOP1OH4N = IGNORE ;
ISOP1N4OH = IGNORE ;
ISOP3N4OH = IGNORE ;
ISOP1OH2OOH = IGNORE ;
ISOP3OOH4OH = IGNORE ;
ISOP1OH4OOH = IGNORE ;
ISOP1OOH4OH = IGNORE ;
IHPOO1 = IGNORE ; {lumped ISOP1OH2OOH3OH4OO and ISOP1OH2OO3OH4OOH}
IHPOO2 = IGNORE ; {lumped ISOP1OO2OH3OOH4OH and ISOP1OOH2OH3OO4OH}
IHPOO3 = IGNORE ; {lumped ISOP1OH2OO3OOH4OH and ISOP1OH2OOH3OO4OH}
IEPOXt = IGNORE ; {trans-beta-IEPOX}
IEPOXc = IGNORE ; {cis-beta-IEPOX}
IEPOXD = IGNORE ; {lumped delta-IEPOX}
ICHE = IGNORE ; {isoprene-carbonyl-hydroxy-epoxide}
IEPOXAOO = IGNORE ; {ISOP1CO2OO3OH4OH}
IEPOXBOO = IGNORE ; {ISOP1OH2OH3OO4CO}
MPAN = IGNORE ;
ISOP1CO4CO = IGNORE ;
C4HVP1 = IGNORE ; 
C4HVP2 = IGNORE ;
HPALD1OO = IGNORE ; {lumped ISOP1CO3OO4OOH and ISOP121CO3OO4N}
HPALD2OO = IGNORE ; {lumped ISOP1OOH2OO4CO and ISOP1N2OO344CO}
ISOPNOO1 = IGNORE ; {lumped ISOP1OH2N3OO4OH and ISOP1OH2N3OH4OO}
ISOPNOO2 = IGNORE ; {lumped ISOP1OH2OO3N4OH and ISOP1OO2OH3N4OH}
INO2B = IGNORE ; {beta-nitrooxy peroxy radicals from isoprene + NO3}
INO2D = IGNORE ; {delta-nitrooxy peroxy radicals from isoprene + NO3}
INPB = IGNORE ; {beta-nitrooxy hydroperoxides from isoprene + NO3}
INPD = IGNORE ; {delta-nitrooxy hydroperoxides from isoprene + NO3}
ISOP1N4CO = IGNORE ;
ISOP1CO4N = IGNORE ;
ISOP3CO4N = IGNORE ;
IHNB = IGNORE ; {beta-nitrooxy alcohols from isoprene + NO3}
INO = IGNORE ; {ISOP1N4O}
IDN = IGNORE ; {dinitrates from isoprene + NO3}
IDHNBOO = IGNORE ; {lumped hydroxy-peroxy radicals from IHNB}
IDHNDOO1 = IGNORE ; {lumped ISOP1N2OH3OO4OH and ISOP1N2OO3OH4OH}
IDHNDOO2 = IGNORE ; {lumped ISOP1OH2OO3OH4N and ISOP1OH2OH3OO4N}
IHPNBOO = IGNORE ; {lumped hydroxy-peroxy radicals from IPNB}
IHPNDOO = IGNORE ; {lumped hydroxy-peroxy radicals from IPND}
IHNE = IGNORE ; {isoprene-hydroxy-nitrooxy-epoxide}
ICNE = IGNORE ; {isoprene-carbonyl-nitrooxy-epoxide}
IHNEOO = IGNORE ; {lumped peroxy radicals from IHNE}
ICN1OO = IGNORE ; {ISOP1N2OH3OO4CO}
ICN2OO = IGNORE ; {ISOP1CO2OO3OH4N}
ICN3OO = IGNORE ; {ISOP1OH2OO3CO4N}
ICN4OO = IGNORE ; {ISOP1N4CO4OO}
ICN5OO = IGNORE ; {ISOP1CO1OO4N}
IHNPE = IGNORE ; {lumped tetrafunctionalized isoprene: hydroxy-nitrooxy-hydroperoxy-epoxide}
ICHNP = IGNORE ; {lumped tetrafunctionalized isoprene: carbonyl-hydroxy-nitrooxy-hydroperoxide}
IHNDP = IGNORE ; {lumped tetrafunctionalized isoprene: hydroxy-nitrooxy-dihydroperoxide}
IDHDN = IGNORE ; {lumped tetrafunctionalized isoprene: dihydroxy-dinitrate}
IDHPN = IGNORE ; {lumped tetrafunctionalized isoprene: dihydroxy-hydroperoxy-nitrate}
IDHCN = IGNORE ; {lumped tetrafunctionalized isoprene: dihydroxy-carbonyl-nitrate}
ICPDH = IGNORE ; {lumped tetrafunctionalized isoprene: carbonyl-hydroperoxy-diol}
IDHDP = IGNORE ; {lumped tetrafunctionalized isoprene: dihydroxy-dihydroperoxide}
IDHPE = IGNORE ; {lumped tetrafunctionalized isoprene: dihydroxy-hydroperxy-epoxide}
IHPDN = IGNORE ; {lumped tetrafunctionalized isoprene: hydroxy-hydroperoxy-dinitrate}
ICHDN = IGNORE ; {lumped tetrafunctionalized isoprene: carbonyl-hydroxy-dinitrate}
IDCHP = IGNORE ; {lumped tetrafunctionalized isoprene: hydroxy-hydroperoxy-dialdehyde}
MVK = IGNORE ; {methyl vinyl ketone}
MACR = IGNORE ; {methacrolein}
MVK3CO4N = IGNORE ; 
MVK3OH4OH = IGNORE ;
MACR2OH3OH = IGNORE ;
MVK3OOH4OH = IGNORE ;
MACR2OOH3OH = IGNORE ;
MVK3CO4OH = IGNORE ;
MVK3OH4CO = IGNORE ;
MACR2OH3CO = IGNORE ;
MVKOHOO = IGNORE ; {lumped MVK-OH-OO radicals}
MVK3OH4OOH = IGNORE ;
MVK3OH4N = IGNORE ;
MVK3N4OH = IGNORE ;
MCROHOO = IGNORE ; {MACR2OO3OH}
MACR1OO = IGNORE ;
MACR1OOH = IGNORE ;
MACR1OH = IGNORE ;
MACR2N3OH = IGNORE ;
MVK3OOH4CO = IGNORE ;
MVKENOL = IGNORE ;
MCRENOL = IGNORE ;
MVK3OOH4N = IGNORE ;
MACR2OOH3N = IGNORE ;
MACR2OH3N = IGNORE ;
MACR2OH3OOH = IGNORE ;
HNO3 = IGNORE ;

#DEFFIX
M = IGNORE ;
O2 = IGNORE ;
OH = IGNORE ;
NO = IGNORE ;
HO2 = IGNORE ;
NO2 = IGNORE ;
NO3 = IGNORE ;
MO2 = IGNORE ;
H2O = IGNORE ;
O3 = IGNORE ;