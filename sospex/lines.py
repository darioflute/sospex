import collections

"""
ATOMIC AND IONIC FINE-STRUCTURE LINES (2-205 UM) FROM HTTPS://WWW.MPE.MPG.DE/IR/ISO/LINELISTS/INDEX.HTML
H2, HD, H2O, OH LINES FROM HTTP://WWW.MPE-GARCHING.MPG.DE/ISO/LINELISTS/MOLECULAR.HTML
CO LINES FROM WW.MPE-GARCHING.MPG.DE/ISO/LINELISTS/CO.HTML
[CIII] line from https://physics.nist.gov/PhysRefData/ASD/lines_form.html
"""
def vac2air(vac):
    # Formula from Morton (1991, ApJS, 77, 119) to transform vacuum into air wavelengths
    # Check also: https://classic.sdss.org/dr7/products/spectra/vacwavelength.html
    # vac in Angstrom
    vac *= 1.e4 # from um to A
    air = vac / (1.0 + 2.735182E-4 + 131.4182 / vac**2 + 2.76249E8 / vac**4)
    return air*1e-4

def define_lines(reference='vacuum'):
    _two = '$_2$'
    _three = '$_3'
    _two18 = '$_2$$^{18}$'
    _1h3h = '$_{1/2+3/2}$'
    _1h1h = '$_{1/2+1/2}$'
    _3h3h = '$_{3/2+3/2}$'
    _3h1h = '$_{3/2+1/2}$'
    alpha = u'\u03B1'
    beta = u'\u03B2'
    gamma = u'\u03B3'
    delta = u'\u03B4'
    eps = u'\u03B5'
    if reference == 'vacuum':
        return collections.OrderedDict([
            ('[SiVI]   2P1/2-2P3/2',['[SiVI]',  +1.963410]),
            ('[SiVII] 3P1-3P2',['[SiVII]', +2.483340]),
            ('[SiIX]   3P2-3P1',['[SiIX]',  +2.584230]),
            ('[AlV]    2P1/2-2P3/2',['[AlV]',   +2.905190]),
            ('[CoII]   a5F5-a3F4',['[CoII]',  +2.984600]),
            ('[MgVIII] 2P3/2-2P1/2',['[MgVIII]',+3.027950]),
            ('[NiI]    a1D2-a3D3',['[NiI]',   +3.120000]),
            ('[KVII]   2P3/2-2P1/2',['[KVII]',  +3.190530]),
            ('[CaIV]   2P1/2-2P3/2',['[CaIV]',  +3.206710]),
            ('[AlVI]   3P1-3P2',['[AlVI]',  +3.659710]),
            ('[AlVIII] 3P2-3P1',['[AlVIII]',+3.690000]),
            ('[SiIX]   3P1-3P0',['[SiIX]',  +3.935700]),
            ('[NiI]    a1D2-a3D2',['[NiI]',   +3.952400]),
            ('[FeII]   a4F5/2-a6D7/2',['[FeII]',  +4.076319]),
            ('[FeII]   a4F3/2-a6D5/2',['[FeII]',  +4.081906]),
            ('[CaVII]  3P2-3P1',['[CaVII]', +4.085800]),
            ('[FeII]   a4F7/2-a6D9/2',['[FeII]',  +4.114990]),
            ('[CaV]    3P1-3P2',['[CaV]',   +4.159370]),
            ('[FeII]   a4F3/2-a6D3/2',['[FeII]',  +4.434839]),
            ('[MgIV]   2P1/2-2P3/2',['[MgIV]',  +4.486680]),
            ('[ArVI]   2P3/2-2P1/2',['[ArVI]',  +4.529520]),
            ('[FeII]   a4F5/2-a6D5/2',['[FeII]',  +4.607664]),
            ('[KIII]   2P1/2-2P3/2',['[KIII]',  +4.618020]),
            ('[FeII]   a4F3/2-a6D1/2',['[FeII]',  +4.671945]),
            ('[NaVII]  2P3/2-2P1/2',['[NaVII]', +4.684720]),
            ('[FeII]   a4F7/2-a6D7/2',['[FeII]',  +4.889137]),
            ('[FeII]   a4F5/2-a6D3/2',['[FeII]',  +5.062350]),
            ('[FeII]   a4F9/2-a6D9/2',['[FeII]',  +5.340169]),
            ('[MgVII]  3P2-3P1',['[MgVII]', +5.503200]),
            ('[KVI]    3P2-3P1',['[KVI]',   +5.575000]),
            ('[MgV]    3P1-3P2',['[MgV]',   +5.609850]),
            ('[FeII]   a4F7/2-a6D5/2',['[FeII]',  +5.673905]),
            ('[AlVIII] 3P1-3P0',['[AlVIII]',+5.850000]),
            ('[NiI]    a1D2-a3D1',['[NiI]',   +5.893300]),
            ('[KIV]    3P1-3P2',['[KIV]',   +5.982000]),
            ('[CaVII]  3P1-3P0',['[CaVII]', +6.154000]),
            ('[SiVII]  3P0-3P1',['[SiVII]', +6.492200]),
            ('[NiII]   2D3/2-2D5/2',['[NiII]',  +6.636000]),
            ('[ClV]    2P3/2-2P1/2',['[ClV]',   +6.706670]),
            ('[FeII]   a4F9/2-a6D7/2',['[FeII]',  +6.721283]),
            ('[ArII]   2P1/2-2P3/2',['[ArII]',  +6.985274]),
            ('[NaIII]  2P1/2-2P3/2',['[NaIII]', +7.317700]),
            ('[NiI]    a3F3-a3F4',['[NiI]',   +7.506700]),
            ('[NeVI]   2P3/2-2P1/2',['[NeVI]',  +7.652400]),
            ('[FeVII]  3F4-3F3',['[FeVII]', +7.814500]),
            ('[ArV]    3P2-3P1',['[ArV]',   +7.901600]),
            ('[NaVI]   3P2-3P1',['[NaVI]',  +8.610590]),
            ('[KVI]    3P1-3P0',['[KVI]',   +8.829900]),
            ('[ArIII]  3P1-3P2',['[ArIII]', +8.991380]),
            ('[MgVII]  3P1-3P0',['[MgVII]', +9.009000]),
            ('[NaIV]   3P1-3P2',['[NaIV]',  +9.041000]),
            ('[AlVI]   3P0-3P1',['[AlVI]',  +9.116000]),
            ('[FeVII]  3F3-3F2',['[FeVII]', +9.526700]),
            ('[SIV]    2P3/2-2P1/2',['[SIV]',  +10.510500]),
            ('[CoII]   a3F3-a3F4',['[CoII]', +10.521000]),
            ('[NiII]   4F7/2-4F9/2',['[NiII]', +10.682200]),
            ('[NiI]    a3F2-a3F3',['[NiI]',  +11.307500]),
            ('[ClI]    2P1/2-2P3/2',['[ClI]',  +11.333347]),
            ('[CaV]    3P0-3P1',['[CaV]',  +11.482000]),
            ('[ClIV]   3P2-3P1',['[ClIV]', +11.761900]),
            ('[CoIII]  a4F7/2-a4F9/2',['[CoIII]',+11.890000]),
            ('[NiI]    a3D1-a3D2',['[NiI]',  +12.001000]),
            ('[CoI]    a4F7/2-a4F9/2',['[CoI]',  +12.255000]),
            ('[NiII]   4F5/2-4F7/2',['[NiII]', +12.728800]),
            ('[NeII]   2P1/2-2P3/2',['[NeII]', +12.813550]),
            ('[ArV]    3P1-3P0',['[ArV]',  +13.102200]),
            ('[FV]     2P3/2-2P1/2',['[FV]',   +13.432000]),
            ('[MgV]    3P0-3P1',['[MgV]',  +13.521300]),
            ('[NeV]    3P2-3P1',['[NeV]',  +14.321700]),
            ('[ClII]   3P1-3P2',['[ClII]', +14.367800]),
            ('[NaVI]   3P1-3P0',['[NaVI]', +14.396400]),
            ('[CoII]   a5F4-a5F5',['[CoII]', +14.740000]),
            ('[NiI]    a3D2-a3D3',['[NiI]',  +14.814200]),
            ('[CoI]    b4F7/2-b4F9/2',['[CoI]',  +15.155000]),
            ('[KIV]    3P0-3P1',['[KIV]',  +15.390000]),
            ('[CoII]   a3F2-a3F3',['[CoII]', +15.460000]),
            ('[NeIII]  3P1-3P2',['[NeIII]',+15.555100]),
            ('[CoIII]  a4F5/2-a4F7/2',['[CoIII]',+16.390000]),
            ('[CoI]    a4F5/2-a4F7/2',['[CoI]',  +16.925000]),
            ('[PIII]   2P3/2-2P1/2',['[PIII]', +17.885000]),
            ('[FeII]   a4F7/2-a4F9/2',['[FeII]', +17.935950]),
            ('[NiII]   4F3/2-4F5/2',['[NiII]', +18.240500]),
            ('[CoI]    b4F5/2-b4F7/2',['[CoI]',  +18.264000]),
            ('[SIII]   3P2-3P1',['[SIII]', +18.713000]),
            ('[CoII]   a5F3-a5F4',['[CoII]', +18.804000]),
            ('[ClIV]   3P1-3P0',['[ClIV]', +20.310700]),
            ('[NaIV]   3P0-3P1',['[NaIV]', +21.290000]),
            ('[ArIII]  3P0-3P1',['[ArIII]',+21.830200]),
            ('[FeIII]  5D3-5D4',['[FeIII]',+22.925000]),
            ('[FeI]    a5D3-a5D4',['[FeI]',  +24.042330]),
            ('[CoIII]  a4F3/2-a4F5/2',['[CoIII]',+24.070000]),
            ('[NeV]    3P1-3P0',['[NeV]',  +24.317500]),
            ('[FeII]   a4F5/2-a4F7/2',['[FeII]', +24.519250]),
            ('[FI]     2P1/2-2P3/2',['[FI]',   +24.747500]),
            ('[CoI]    a4F3/2-a4F5/2',['[CoI]',  +24.845000]),
            ('[SI]     3P1-3P2',['[SI]',   +25.249000]),
            ('[CoII]   a5F2-a5F3',['[CoII]', +25.681000]),
            ('[FIV]    3P2-3P1',['[FIV]',  +25.830000]),
            ('[OIV]    2P3/2-2P1/2',['[OIV]',  +25.890300]),
            ('[CoI]    b4F3/2-b4F5/2',['[CoI]',  +25.930000]),
            ('[FeII]   a6D7/2-a6D9/2',['[FeII]', +25.988290]),
            ('[FII]    3P1-3P2',['[FII]',  +29.330000]),
            ('[PII]    3P2-3P1',['[PII]',  +32.870000]),
            ('[FeIII]  5D2-5D3',['[FeIII]',+33.038400]),
            ('[ClII]   3P0-3P1',['[ClII]', +33.281000]),
            ('[SIII]   3P1-3P0',['[SIII]', +33.481000]),
            ('[FeI]    a5D2-a5D3',['[FeI]',  +34.713300]),
            ('[SiII]   2P3/2-2P1/2',['[SiII]', +34.815200]),
            ('[FeII]   a6D5/2-a6D7/2',['[FeII]', +35.348650]),
            ('[FeII]   a4F3/2-a4F5/2',['[FeII]', +35.777400]),
            ('[NeIII]  3P0-3P1',['[NeIII]',+36.013500]),
            ('[CoII]   a5F1-a5F2',['[CoII]', +39.274000]),
            ('[FIV]    3P1-3P0',['[FIV]',  +44.070000]),
            ('[FeII]   a6D3/2-a6D5/2',['[FeII]', +51.300440]),
            ('[FeIII]  5D1-5D2',['[FeIII]',+51.680000]),
            ('[OIII]   3P2-3P1',['[OIII]', +51.814500]),
            ('[FeI]    a5D1-a5D2',['[FeI]',  +54.310930]),
            ('[SI]     3P0-3P1',['[SI]',   +56.311000]),
            ('[NIII]   2P3/2-2P1/2',['[NIII]', +57.317000]),
            ('[PII]    3P1-3P0',['[PII]',  +60.640000]),
            ('[OI]     3P1-3P2',['[OI]',   +63.183705]),
            ('[FII]    3P0-3P1',['[FII]',  +67.200000]),
            ('[SiI]    3P2-3P1',['[SiI]',  +68.473000]),
            ('[FeII]   a6D1/2-a6D3/2',['[FeII]', +87.384400]),
            ('[OIII]   3P1-3P0',['[OIII]', +88.356000]),
            ('[AlI]    2P3/2-2P1/2',['[AlI]',  +89.237000]),
            ('[FeIII]  5D0-5D1',['[FeIII]',105.370000]),
            ('[FeI]    a5D0-a5D1',['[FeI]',  111.182800]),
            ('[NII]    3P2-3P1',['[NII]',  121.897570]),
            ('[SiI]    3P1-3P0',['[SiI]',  129.681730]),
            ('[OI]     3P0-3P1',['[OI]',   145.525439]),
            ('[CII]    2P3/2-2P1/2',['[CII]',  157.740900]),
            ('[CIII]   1s22s2p/3P0/1',['[CIII]', 177.4308]),
            ('[NII]    3P1-3P0',['[NII]',  205.178230]),
            ('CO_J=13->12',['CO',  200.272476]),	         
            ('CO_J=14->13',['CO',  185.999313]),	         
            ('CO_J=15->14',['CO',  173.631434]),	         
            ('CO_J=16->15',['CO',  162.811630]),	         
            ('CO_J=11->16',['CO',  153.266708]),	         
            ('CO_J=18->17',['CO',  144.784195]),	         
            ('CO_J=19->18',['CO',  137.196332]),	         
            ('CO_J=20->19',['CO',  130.368927]),	         
            ('CO_J=21->20',['CO',  124.193352]),	         
            ('CO_J=22->21',['CO',  118.580719]),	         
            ('CO_J=23->22',['CO',  113.457603]),	         
            ('CO_J=24->23',['CO',  108.762810]),	         
            ('CO_J=25->24',['CO',  104.444952]),	         
            ('CO_J=26->25',['CO',  100.460533]),                   
            ('CO_J=27->26',['CO',   96.772514]),	         
            ('CO_J=28->27',['CO',   93.349123]),	         
            ('CO_J=29->28',['CO',   90.163002]),	         
            ('CO_J=30->29',['CO',   87.190422]),	         
            ('CO_J=31->30',['CO',   84.410721]),	         
            ('CO_J=32->31',['CO',   81.805809]),	         
            ('CO_J=33->32',['CO',   79.359810]),	         
            ('CO_J=34->33',['CO',   77.058693]),	         
            ('CO_J=35->34',['CO',   74.890053]),	         
            ('CO_J=36->35',['CO',   72.84272]),                
            ('CO_J=37->36',['CO',   70.90710]),                
            ('CO_J=38->37',['CO',   69.07426]),                
            ('CO_J=39->38',['CO',   67.33630]),                
            ('CO_J=40->39',['CO',   65.68611]),
            ('12CO(1-0)', ['$^{12}$CO$_{1-0}$', 2600.7576]),
            ('12CO(2-1)', ['$^{12}$CO$_{2-1}$', 1300.4036]),
            ('13CO(1-0)', ['$^{13}$CO$_{1-0}$', 2720.4063]),        
            ('13CO(2-1)', ['$^{13}$CO$_{2-1}$', 1360.2280]),        
            ('H2;(1,0)Q(1)',['H'+_two,2.4065914]),                                  
            ('H2;(2,1)Q(1)',['H'+_two,2.5510]),                          
            ('H2;(2,1)Q(3)',['H'+_two,2.5698]),                          
            ('H2;(2,1)Q(5)',['H'+_two,2.6040]),                          
            ('H2;(1,0)O(2)',['H'+_two,2.6269]),                          
            ('H2;(1,0)O(3)',['H'+_two,2.8025]),                          
            ('H2;(1,0)O(5)',['H'+_two,3.2350]),                          
            ('H2;(0,0)S(13)',['H'+_two,3.8468]),                          
            ('H2;(0,0)S(11)',['H'+_two,4.1813]),                          
            ('H2;(0,0)S(10)',['H'+_two,4.4099]),                          
            ('H2;(0,0)S(9)',['H'+_two,4.69461]),                          
            ('H2;(0,0)S(8)',['H'+_two,5.05303]),                          
            ('H2;(0,0)S(7)',['H'+_two,5.51116]),                          
            ('H2;(0,0)S(6)',['H'+_two,6.10856]),                          
            ('H2;(0,0)S(5)',['H'+_two,6.90952]),                          
            ('H2;(0,0)S(4)',['H'+_two,8.02505]),                          
            ('H2;(0,0)S(3)',['H'+_two,9.66491]),                          
            ('HD;(0,0)R(10)',['HD',11.57346]),                          
            ('H2;(0,0)S(2)',['H'+_two,12.27861]),                          
            ('HD;(0,0)R(9)',['HD',12.47181]),                          
            ('HD;(0,0)R(8)',['HD',13.59265]),                          
            ('HD;(0,0)R(7)',['HD',15.25104]),                          
            ('HD;(0,0)R(6)',['HD',16.89381]),                          
            ('H2;(0,0)S(1)',['H'+_two,17.03483]),                          
            ('p-H2O;5_51-4_04',['p-H'+_two+'O',19.2300]),                          
            ('HD;(0,0)R(5)',['HD',19.43100]),                          
            ('o-H2O;5_50-4_23',['o-H'+_two+'O',22.6391]),                          
            ('HD;(0,0)R(4)',['HD',23.03376]),                          
            ('OH1/2-3/2;9/2-7/2',['OH'+_1h3h,24.614]),                   
            ('OH1/2-3/2;9/2-7/2',['OH'+_1h3h,24.642]),                   
            ('o-H2O;5_41-4_14',['o-H'+_two+'O',25.9402]),                          
            ('H2;(0,0)S(0)',['H'+_two,28.21883]),                          
            ('HD;(0,0)R(3)',['HD',28.50197]),                          
            ('p-H2O;4_40-3_13',['p-H'+_two+'O',28.9138]),                          
            ('OH1/2-3/2;7/2-5/2',['OH'+_1h3h,28.939]),                   
            ('OH1/2-3/2;7/2-5/2',['OH'+_1h3h,28.940]),                   
            ('o-H2O;7_25-6_16',['o-H'+_two+'O',29.8363]),                          
            ('p-H2O;5_42-4_13',['p-H'+_two+'O',29.8849]),                          
            ('o-H2O;4_41-3_12',['o-H'+_two+'O',31.7715]),                          
            ('OH1/2-3/2;5/2--3/2+3-2',['OH'+_1h3h,34.6034]),                          
            ('OH1/2-3/2;5/2--3/2+2-1',['OH'+_1h3h,34.6034]),                          
            ('OH1/2-3/2;5/2+-3/2-3-2',['OH'+_1h3h,34.6294]),                          
            ('OH1/2-3/2;5/2+-3/2-2-1',['OH'+_1h3h,34.6293]),                          
            ('p-H2O;5_33-4_04',['p-H'+_two+'O',35.4710]),                          
            ('HD;(0,0)R(2)',['HD',37.70155]),                          
            ('o-H2O;4_41-4_14',['o-H'+_two+'O',37.9839]),                          
            ('OH1/2-3/2;9/2--9/2+5-5',['OH'+_1h3h,39.5139]),                          
            ('OH1/2-3/2;9/2--9/2+4-4',['OH'+_1h3h,39.5144]),                          
            ('OH1/2-3/2;9/2+-9/2-5-5',['OH'+_1h3h,39.6378]),                          
            ('OH1/2-3/2;9/2+-9/2-4-4',['OH'+_1h3h,39.6380]),                          
            ('o-H2O;6_43-5_32',['o-H'+_two+'O',40.3367]),                          
            ('o-H2O;4_32-3_03',['o-H'+_two+'O',40.6904]),                          
            ('OH1/2-3/2;7/2+-7/2-4-4',['OH'+_1h3h,43.9497]),                          
            ('OH1/2-3/2;7/2+-7/2-3-3',['OH'+_1h3h,43.9503]),                          
            ('OH1/2-3/2;7/2--7/2+4-4',['OH'+_1h3h,44.0723]),                          
            ('OH1/2-3/2;7/2--7/2+3-3',['OH'+_1h3h,44.0724]),                          
            ('o-H2O;5_23-4_14',['o-H'+_two+'O',45.1112]),                          
            ('p-H2O;3_31-2_02',['p-H'+_two+'O',46.4835]),                          
            ('p-H2O;4_40-4_13',['p-H'+_two+'O',47.0285]),                          
            ('o-H2O;5_32_4_23',['o-H'+_two+'O',47.9723]),                          
            ('OH1/2-3/2;5/2--5/2+3-3',['OH'+_1h3h,48.7040]),                          
            ('OH1/2-3/2;5/2--5/2+2-2',['OH'+_1h3h,48.7044]),                          
            ('OH1/2-3/2;5/2+-5/2-3-3',['OH'+_1h3h,48.8168]),                          
            ('OH1/2-3/2;5/2+-5/2-2-2',['OH'+_1h3h,48.8168]),                          
            ('p-H2O;4_40-3_31',['p-H'+_two+'O',49.2808]),                          
            ('o-H2O;4_41-3_30',['o-H'+_two+'O',49.3357]),                          
            ('OH3/2-3/2;11/2+-9/2-6-5',['OH'+_3h3h,52.9339]),                          
            ('OH3/2-3/2;11/2+-9/2-5-4',['OH'+_3h3h,52.9338]),                          
            ('OH3/2-3/2;11/2--9/2+6-5',['OH'+_3h3h,53.0572]),                          
            ('OH3/2-3/2;11/2--9/2+5-4',['OH'+_3h3h,53.0571]),                          
            ('OH1/2-3/2;3/2--3/2+2-2',['OH'+_1h3h,53.3513]),                          
            ('OH1/2-3/2;3/2--3/2+1-1',['OH'+_1h3h,53.3509]),                          
            ('OH1/2-3/2;3/2+-3/2-2-2',['OH'+_1h3h,53.2614]),                          
            ('OH1/2-3/2;3/2+-3/2-1-1',['OH'+_1h3h,53.2616]),                          
            ('OH1/2-1/2;9/2--7/2+5-4',['OH'+_1h1h,55.9497]),                          
            ('OH1/2-1/2;9/2--7/2+4-3',['OH'+_1h1h,55.9497]),                          
            ('OH1/2-1/2;9/2+-7/2-5-4',['OH'+_1h1h,55.8908]),                          
            ('OH1/2-1/2;9/2+-7/2-4-3',['OH'+_1h1h,55.8909]),                          
            ('HD;(0,0)R(1)',['HD',56.22975]),                          
            ('p-H2O;4_31-3_22',['p-H'+_two+'O',56.3242]),                          
            ('p-H2O;4_22-3_13',['p-H'+_two+'O',57.6361]),                          
            ('o-H2O;4_32-3_21',['o-H'+_two+'O',58.6982]),                          
            ('p-H2O;4_31-4_04',['p-H'+_two+'O',61.8084]),                          
            ('p-H2O;8_08-7_17',['p-H'+_two+'O',63.4570]),                          
            ('OH3/2-3/2;9/2--7/2+5-4',['OH'+_3h3h,65.1315]),                          
            ('OH3/2-3/2;9/2--7/2+4-3',['OH'+_3h3h,65.1314]),                          
            ('OH3/2-3/2;9/2+-7/2-5-4',['OH'+_3h3h,65.2788]),                          
            ('OH3/2-3/2;9/2+-7/2-4-3',['OH'+_3h3h,65.2786]),                          
            ('o-H2O;3_30-2_21',['o-H'+_two+'O',66.4372]),                          
            ('p-H2O;3_31-2_20',['p-H'+_two+'O',67.0886]),                          
            ('o-H2O;3_30-3_03',['o-H'+_two+'O',67.2689]),                          
            ('p-H2O;5_24-4_13',['p-H'+_two+'O',71.0662]),                          
            ('OH1/2-1/2;7/2--5/2+4-3',['OH'+_1h1h,71.170752]),                          
            ('OH1/2-1/2;7/2--5/2+3-2',['OH'+_1h1h,71.170850]),                          
            ('OH1/2-1/2;7/2+-5/2-4-3',['OH'+_1h1h,71.215827]),                          
            ('OH1/2-1/2;7/2+-5/2-3-2',['OH'+_1h1h,71.215869]),                          
            ('o-H2O;7_07-6_16',['o-H'+_two+'O',71.9460]),                          
            ('o-H2O;3_21-2_12',['o-H'+_two+'O',75.3804]),                          
            ('o-H2O;4_23-3_12',['o-H'+_two+'O',78.7414]),                          
            ('OH1/2-3/2;1/2--3/2+1-2',['OH'+_1h3h,79.117296]),                          
            ('OH1/2-3/2;1/2--3/2+0-1',['OH'+_1h3h,79.118034]),                          
            ('OH1/2-3/2;1/2+-3/2-1-2',['OH'+_1h3h,79.181726]),                          
            ('OH1/2-3/2;1/2+-3/2-0-1',['OH'+_1h3h,79.180928]),                          
            ('o-H2O;6_16-5_05',['o-H'+_two+'O',82.0304]),                          
            ('p-H2O;6_06-5_15',['p-H'+_two+'O',83.2831]),                          
            ('OH3/2-3/2;7/2+-5/2-4-3',['OH'+_3h3h,84.420391]),                          
            ('OH3/2-3/2;7/2+-5/2-3-2',['OH'+_3h3h,84.419938]),                          
            ('OH3/2-3/2;7/2--5/2+4-3',['OH'+_3h3h,84.596825]),                          
            ('OH3/2-3/2;7/2--5/2+3-2',['OH'+_3h3h,84.596312]),                          
            ('p-H2O;3_22-2_11',['p-H'+_two+'O',89.9878]),                          
            ('p-H2O;5_15-4_04',['p-H'+_two+'O',95.6263]),                          
            ('OH3/2-1/2;3/2+-5/2-2-3',['OH'+_3h1h,96.312205]),                          
            ('OH3/2-1/2;3/2+-5/2-1-2',['OH'+_3h1h,96.313812]),                          
            ('OH3/2-1/2;3/2--5/2+2-3',['OH'+_3h1h,96.367469]),                          
            ('OH3/2-1/2;3/2--5/2+1-2',['OH'+_3h1h,96.367404]),                          
            ('OH1/2-1/2;5/2--3/2+3-2',['OH'+_1h1h,98.724877]),                          
            ('OH1/2-1/2;5/2--3/2+2-1',['OH'+_1h1h,98.724920]),                          
            ('OH1/2-1/2;5/2+-3/2-3-2',['OH'+_1h1h,98.736890]),                          
            ('OH1/2-1/2;5/2+-3/2-2-1',['OH'+_1h1h,98.737085]),                          
            ('o-H2O;5_05-4_14',['o-H'+_two+'O',99.4924]),                          
            ('o-H2O;5_14-4_23',['o-H'+_two+'O',100.9127]),                          
            ('p-H2O;2_20-1_11',['p-H'+_two+'O',100.9828]),  
            ('OH3_3--2_2',['OH 3$_3-$_2$_2$', 101.70]),         # Gonzales-Alonso, Smith et al.,ApJ 613:247               
            ('OH3_4--2_3',['OH 3$_4-$_2$_3$', 101.92]),         # Gonzales-Alonso, Smith et al.,ApJ 613:247                           
            ('p-H218O;2_20-1_11',['p-$^{18}$H$_2$O', 102.004]),                   
            ('o-H2O;2_21-1_10',['o-H'+_two+'O',108.0730]),                          
            ('o-H218O;2_21-1_10',['o-$^{18}$H$_2$O',109.347]),                   
            ('HD;(0,0)R(0)',['HD',112.07251]),                          
            ('o-H2O;4_14-3_03',['o-H'+_two+'O',113.5366]),                          
            ('OH1/2-3/2;5/2+-7/2-3-4',['OH'+_1h3h,115.1530]),                          
            ('OH1/2-3/2;5/2+-7/2-2-3',['OH'+_1h3h,115.1541]),                          
            ('OH1/2-3/2;5/2--7/2+3-4',['OH'+_1h3h,115.3858]),                          
            ('OH1/2-3/2;5/2--7/2+2-3',['OH'+_1h3h,115.3890]),                          
            ('o-H2O;7_34-6_43',['o-H'+_two+'O',116.7836]),                          
            ('OH3/2-3/2;5/2--3/2+3-2',['OH'+_3h3h,119.234182]),                          
            ('OH3/2-3/2;5/2--3/2+2-1',['OH'+_3h3h,119.232438]),                          
            ('OH3/2-3/2;5/2+-3/2-3-2',['OH'+_3h3h,119.441669]),                          
            ('OH3/2-3/2;5/2+-3/2-2-1',['OH'+_3h3h,119.439806]),                          
            ('18OH3/2-3/2;5/2+-3/2-3-2',['$^{18}$OH'+_3h3h,119.9659]),                  
            ('18OH3/2-3/2;5/2+-3/2-2-1',['$^{18}$OH'+_3h3h,119.9641]),                  
            ('18OH3/2-3/2;5/2--3/2+3-2',['$^{18}$OH'+_3h3h,120.1726]),                  
            ('18OH3/2-3/2;5/2--3/2+2-1',['$^{18}$OH'+_3h3h,120.1707]),                  
            ('o-H2O;4_32-4_23',['o-H'+_two+'O',121.7191]),                          
            ('p-H2O;4_04-3_13',['p-H'+_two+'O',125.3534]),                          
            ('p-H2O;3_31-3_22',['p-H'+_two+'O',126.7126]),                          
            ('o-H2O;4_23-4_14',['o-H'+_two+'O',132.4070]),                          
            ('OH1/2-3/2;7/2--9/2+4-5',['OH'+_1h3h,134.8449]),                          
            ('OH1/2-3/2;7/2--9/2+3-4',['OH'+_1h3h,134.8476]),                          
            ('OH1/2-3/2;7/2+-9/2-4-5',['OH'+_1h3h,134.9642]),                          
            ('OH1/2-3/2;7/2+-9/2-3-4',['OH'+_1h3h,134.9696]),                          
            ('o-H2O;5_14-5_05',['o-H'+_two+'O',134.9346]),                          
            ('o-H2O;3_30-3_21',['o-H'+_two+'O',136.4944]),                          
            ('p-H2O;3_13-2_02',['p-H'+_two+'O',138.5272]),                          
            ('p-H2O;4_13-3_22',['p-H'+_two+'O',144.5181]),                          
            ('13CO;18-17',['13CO',151.4312]),         
            ('p-H2O;3_22-3_13',['p-H'+_two+'O',156.1930]),                          
            ('p-H2O;3_31-4_04',['p-H'+_two+'O',158.3090]),                          
            ('OH1/2-1/2;3/2+-1/2-2-1',['OH'+_1h1h,163.124275]),                          
            ('OH1/2-1/2;3/2+-1/2-1-0',['OH'+_1h1h,163.122481]),                          
            ('OH1/2-1/2;3/2--1/2+2-1',['OH'+_1h1h,163.397176]),                          
            ('OH1/2-1/2;3/2--1/2+1-0',['OH'+_1h1h,163.396902]),                          
            ('o-H2O;3_03-2_12',['o-H'+_two+'O',174.6264]),                          
            ('o-H2O;2_12-1_01',['o-H'+_two+'O',179.5265]),                          
            ('o-H2O;2_21-2_12',['o-H'+_two+'O',180.4880]),                          
            ('o-H218O;2_12-1_01',['o-H'+_two18+'O',181.051]),                   
            ('p-H2O;4_13-4_04',['p-H'+_two+'O',187.1104]),
            # Young et al. (2009, http://www.physics.arizona.edu/~young/binaries/TDS/paper.pdf) - absorption line
            ('NH3 j1<-0', ['NH'+_three, 171.6306]),
            ('NH3 j2<-1', ['NH'+_three, 364.2034]),
            ('NH3 j3<-2', ['NH'+_three, 529.13]),
            ('PAH  5.27', ['PAH',  5.27]),  # FRAC_FWHM: 0.034
            ('PAH  5.70', ['PAH',  5.70]),  # FRAC_FWHM: 0.035
            ('PAH  6.22', ['PAH',  6.22]),  # FRAC_FWHM: 0.030
            ('PAH  6.69', ['PAH',  6.69]),  # FRAC_FWHM: 0.07
            ('PAH  7.42', ['PAH',  7.42]),  # FRAC_FWHM: 0.126
            ('PAH  7.60', ['PAH',  7.60]),  # FRAC_FWHM: 0.044
            ('PAH  7.85', ['PAH',  7.85]),  # FRAC_FWHM: 0.053
            ('PAH  8.33', ['PAH',  8.33]),  # FRAC_FWHM: 0.05
            ('PAH  8.61', ['PAH',  8.61]),  # FRAC_FWHM: 0.039
            ('PAH 10.68', ['PAH', 10.68]),  # FRAC_FWHM: 0.02
            ('PAH 11.23', ['PAH', 11.23]),  # FRAC_FWHM: 0.012
            ('PAH 11.33', ['PAH', 11.33]),  # FRAC_FWHM: 0.032
            ('PAH 11.99', ['PAH', 11.99]),  # FRAC_FWHM: 0.045
            ('PAH 12.62', ['PAH', 12.62]),  # FRAC_FWHM: 0.042
            ('PAH 12.69', ['PAH', 12.69]),  # FRAC_FWHM: 0.013
            ('PAH 13.48', ['PAH', 13.48]),  # FRAC_FWHM: 0.04
            ('PAH 14.04', ['PAH', 14.04]),  # FRAC_FWHM: 0.016
            ('PAH 14.19', ['PAH', 14.19]),  # FRAC_FWHM: 0.025
            ('PAH 15.90', ['PAH', 15.90]),  # FRAC_FWHM: 0.02
            ('PAH 16.45', ['PAH', 16.45]),  # FRAC_FWHM: 0.014
            ('PAH 17.04', ['PAH', 17.04]),  # FRAC_FWHM: 0.065
            ('PAH 17.37', ['PAH', 17.37]),  # FRAC_FWHM: 0.012
            ('PAH 17.87', ['PAH', 17.87]),  # FRAC_FWHM: 0.016
            ('PAH 18.92', ['PAH', 18.92]),  # FRAC_FWHM: 0.019
            ('PAH 33.10', ['PAH', 33.10]),  # FRAC_FWHM: 0.05    
            ('HI 21cm', ['HI', 211061.140542]), # 21 cm
            # SDSS lines in vacuum
            ('OVI 1033', ['OVI', 1033.82e-4]),
            ('Ly-alpha 1215', ['Ly' + alpha, 1215.24e-4]),
            ('N-V 1240', ['N-V', 1240.81e-4]),
            ('OI 1304', ['OI', 1304.53e-4]),
            ('Si-IV 1397', ['Si-IV', 1397.61e-4]),
            ('C-IV 1549', ['C-IV', 1549.48e-4]),
            ('He-II 1640', ['He-II', 1640.40e-4]),
            ('O-III 1666', ['O-III', 1665.85e-4]),
            ('Al-III 1857', ['Al-III', 1857.40e-4]),
            ('C-III 1908', ['C-III', 1908.734e-4]),
            ('C-II 2326', ['C-II', 2326.00e-4]),
            ('Ne-IV 2439', ['Ne-IV', 2439.50e-4]),
            ('Mg-II 2799', ['Mg-II', 2799.117e-4]),  # SDSS
            ('Ne-V 3346', ['Ne-V', 3346.79e-4]),
            ('Ne-VI 3426', ['Ne-VI', 3426.85e-4]),
            ('[OII] 3728', ['[OII]', 3728.48e-4]),
            ('[NeIII] 3869', ['[NeIII]', 3869.85e-4]),
            ('H-8 3890', ['H-8', 3890.15e-4]),
            ('H-eps 3971', ['H' + eps, 3971.20e-4]),
            ('H-delta 4102', ['H' + delta, 4102.89e-4]),
            ('H-gamma 4341', ['H' + gamma, 4341.68e-4]),
            ('[OIII] 4363', ['[OIII]', 4363.2e-4]),
            ('H-beta 4862', ['H' + beta, 4862.721e-4]),
            ('[OIII] 4960', ['[OIII]', 4960.295e-4]),
            ('[OIII] 5008', ['[OIII]', 5008.239e-4]),
            ('[NI] 5199', ['[NI]', 5199e-4]),
            ('[NII] 5755', ['[NII]', 5754.6e-4]),
            ('HeI 5877', ['HeI', 5877.29e-4]),
            ('[OI] 6302', ['[OI]', 6302.05e-4]),
            ('[NII] 6549', ['[NII]', 6549.86e-4]),
            ('H-alpha 6564', ['H' + alpha, 6564.614e-4]),
            ('[NII] 6585', ['[NII]', 6585.27e-4]),
            ('SII 6718', ['SII', 6718.29e-4]),
            ('SII 6732', ['SII', 6732.68e-4])
            #('A:Ca(H) 3934', ['A:Ca(H)', 3934.78e-4]),
            #('A:Ca(K) 3969', ['A:Ca(K)', 3969.59e-4]),
            #('A:G-band 4300', ['A:G-band', 4300.4e-4]),
            #('A:Mg-1 5167', ['A:Mg-1', 5167.3222e-4]),
            #('A:Mg-2 5172', ['A:Mg-2', 5172.6847e-4]),
            #('A:Mg-3 5183', ['A:Mg-3', 5183.6046e-4]),
            #('A:Na 5894', ['A:Na', 5894.57e-4]),
            #('A:H-delta 4102', ['A:H' + delta, 4102.89e-4]),
            #('A:H-gamma 4341', ['A:H' + gamma, 4341.68e-4]),
            #('A:H-beta 4862', ['A:H' + beta, 4862.68e-4]),
            #('A:H-alpha 6564', ['A:H' + alpha, 6564.61e-4])
        ])
    else:
        return collections.OrderedDict([
            # SDSS lines in air
            #('OVI 1033', ['OVI', 1033.82e-4]),
            #('Ly-alpha 1215', ['Ly' + alpha, 1215.24e-4]),
            #('N-V 1240', ['N-V', 1240.81e-4]),
            #('OI 1304', ['OI', 1304.53e-4]),
            #('Si-IV 1397', ['Si-IV', 1397.61e-4]),
            #('C-IV 1549', ['C-IV', 1549.48e-4]),
            #('He-II 1640', ['He-II', 1640.40e-4]),
            #('O-III 1666', ['O-III', 1665.85e-4]),
            #('Al-III 1857', ['Al-III', 1857.40e-4]),
            #('C-III 1908', ['C-III', 1908.734e-4]),
            #('C-II 2326', ['C-II', 2326.00e-4]),
            #('Ne-IV 2439', ['Ne-IV', 2439.50e-4]),
            #('Mg-II 2799', ['Mg-II', 2799.117e-4]),  # SDSS
            #('Ne-V 3346', ['Ne-V', 3346.79e-4]),
            #('Ne-VI 3426', ['Ne-VI', 3426.85e-4]),
            #('[OII] 3728', ['[OII]', 3728.48e-4]),
            #('[NeIII] 3869', ['[NeIII]', 3869.85e-4]),
            #('H-8 3890', ['H-8', 3890.15e-4]),
            #('H-eps 3971', ['H' + eps, 3971.20e-4]),
            #('H-delta 4102', ['H' + delta, 4102.89e-4]),
            #('H-gamma 4341', ['H' + gamma, 4341.68e-4]),
            #('[OIII] 4363', ['[OIII]', 4363.2e-4]),
            #('H-beta 4862', ['H' + beta, 4862.721e-4]),
            #('[OIII] 4960', ['[OIII]', 4960.295e-4]),
            #('[OIII] 5008', ['[OIII]', 5008.239e-4]),
            #('[NI] 5199', ['[NI]', 5199e-4]),
            #('[NII] 5755', ['[NII]', 5754.6e-4]),
            #('HeI 5877', ['HeI', 5877.29e-4]),
            #('[OI] 6302', ['[OI]', 6302.05e-4]),
            #('[NII] 6549', ['[NII]', 6549.86e-4]),
            #('H-alpha 6564', ['H' + alpha, 6564.614e-4]),
            #('[NII] 6585', ['[NII]', 6585.27e-4]),
            #('SII 6718', ['SII', 6718.29e-4]),
            #('SII 6732', ['SII', 6732.68e-4]),
            #('A:Ca(H) 3934', ['A:Ca(H)', 3934.78e-4]),
            #('A:Ca(K) 3969', ['A:Ca(K)', 3969.59e-4]),
            #('A:G-band 4300', ['A:G-band', 4300.4e-4]),
            #('A:Mg-1 5167', ['A:Mg-1', 5167.3222e-4]),
            #('A:Mg-2 5172', ['A:Mg-2', 5172.6847e-4]),
            #('A:Mg-3 5183', ['A:Mg-3', 5183.6046e-4]),
            #('A:Na 5894', ['A:Na', 5894.57e-4]),
            #('A:H-delta 4102', ['A:H' + delta, 4102.89e-4]),
            #('A:H-gamma 4341', ['A:H' + gamma, 4341.68e-4]),
            #('A:H-beta 4862', ['A:H' + beta, 4862.68e-4]),
            #('A:H-alpha 6564', ['A:H' + alpha, 6564.61e-4])
            ('OVI 1033', ['OVI', vac2air(1033.82e-4)]),
            ('Ly-alpha 1215', ['Ly' + alpha, vac2air(1215.24e-4)]),
            ('N-V 1240', ['N-V', vac2air(1240.81e-4)]),
            ('OI 1304', ['OI', vac2air(1304.53e-4)]),
            ('CII 1335', ['CII', vac2air(1335.31e-4)]),
            ('Si-IV 1397', ['Si-IV', vac2air(1397.61e-4)]),
            ('C-IV 1549', ['C-IV', vac2air(1549.48e-4)]),
            ('He-II 1640', ['He-II', vac2air(1640.40e-4)]),
            ('O-III 1666', ['O-III', vac2air(1665.85e-4)]),
            ('Al-III 1857', ['Al-III', vac2air(1857.40e-4)]),
            ('C-III 1908', ['C-III', vac2air(1908.734e-4)]),
            ('C-II 2326', ['C-II', vac2air(2326.00e-4)]),
            ('Ne-IV 2439', ['Ne-IV', vac2air(2439.50e-4)]),
            ('Mg-II 2799', ['Mg-II', vac2air(2799.117e-4)]),  # SDSS
            ('Ne-V 3346', ['Ne-V', vac2air(3346.79e-4)]),
            ('Ne-VI 3426', ['Ne-VI', vac2air(3426.85e-4)]),
            ('[OII] 3728', ['[OII]', vac2air(3728.48e-4)]),
            ('[NeIII] 3869', ['[NeIII]', vac2air(3869.85e-4)]),
            ('H-8 3890', ['H-8', vac2air(3890.15e-4)]),
            ('H-eps 3971', ['H' + eps, vac2air(3971.20e-4)]),
            ('H-delta 4102', ['H' + delta, vac2air(4102.89e-4)]),
            ('H-gamma 4341', ['H' + gamma, vac2air(4341.68e-4)]),
            ('[OIII] 4363', ['[OIII]', vac2air(4363.2e-4)]),
            ('H-beta 4862', ['H' + beta, vac2air(4862.68e-4)]),
            ('[OIII] 4960', ['[OIII]', vac2air(4960.30e-4)]),
            ('[OIII] 5008', ['[OIII]', vac2air(5008.24e-4)]),
            ('[NI] 5199', ['[NI]', vac2air(5199e-4)]),
            ('[NII] 5755', ['[NII]', vac2air(5754.6e-4)]),
            ('HeI 5877', ['HeI', vac2air(5877.29e-4)]),
            ('[OI] 6302', ['[OI]', vac2air(6302.05e-4)]),
            ('[NII] 6549', ['[NII]', vac2air(6549.85e-4)]),
            ('H-alpha 6564', ['H' + alpha, vac2air(6564.61e-4)]),
            ('[NII] 6585', ['[NII]', vac2air(6585.28e-4)]),
            ('SII 6718', ['SII', vac2air(6718.29e-4)]),
            ('SII 6732', ['SII', vac2air(6732.67e-4)]),
            ('[SIII] 9068', ['[SIII]', vac2air(9068.600e-4)]),
            ('[SIII] 9531', ['[SIII]', vac2air(9531.100e-4)]),
            ('[CI] 9824',['[CI]', vac2air(9824.130e-4)]),
            ('[CI] 9850',['[CI]', vac2air(9850.260e-4)])
            ])

    