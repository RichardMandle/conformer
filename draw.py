from rdkit import Chem
from PIL import Image, ImageChops, ImageDraw
from rdkit.Chem import Draw, AllChem, rdAbbreviations, rdCoordGen

def abbrev_mol():
    '''
    Function that uses the rdAbbrevations module in rdkit.Chem
    to create some handy space-saving abbreviations.
    
    e.g. CCCCCCCC in smiles - octane - becomes H17C8- 
    
    abbreviations are a matched "definition" string and "smarts" code.   
    Useful SMARTS information: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
    
    This code is quite fragile; we get lots of errors when trying to put 'linkers' (e.g. 'O','S')
    into the alkyl chain, because we then need to know which way roudn the text goes (eg. OC5 vs C5O...)
    
    In the end, we just abbreviate the linear alkyl portions. Anything bigger than C3 and smaller than C100
    
        LinkList - a list of linking groups in SMARTS format. Currently we've stuck to things we are likely to use (O, S, N, C, alkene, alkyne, halogens)
        BCI - Bonus Carbon Index - a way of keeping track of 'extra' carbon atoms from 'C' linkers
        label - the text label we'll use, format is C(n)H(2m{+1 for terminus}). 
        sma - smarts code to abbreviate
        
        defns = a list of labels and smarts codes which is parsed to abbreviations
        defn2a/b - labels for "functional" group abbreviation (CN, NO2, TMS etc.)
        sma2 - smarts codes for defn2a/b, above
        
        returns(abbrevs) a list of abbreviations that RDkit can use.
        

    Generally, new link_list entries will be one of 3 possibilities:
        1st: [AtomSymbol] + ";D2"; this specifies we look for [AtomSymbol] with at least two explicit (non hydrogen) connections.
        2nd: [AtomSymbol] + ";h1"; this implies we have a terminal [AtomSymbol]
        3rd: [AtomSymbol]; this is monovalent, so a halogen (basically ONLY a halogen).
    Why do this? Well, they have different labelling styles!
    '''
    defns=[] # this is a list of our definitions; a smarts string and corresponding label
    
    Link_List = ['O;D2','O;h1',         # oxygen 
                 'S;D2','S;h1',         # silicon
                 'N;h2','N;h1',         # nitrogen
                 'C]=[C',               # alkene
                 'C]#[C',               # alkyne

                 'F','Cl','Br','I',]    # add your linkers here. Do we need more than these?
    #Link_List = [s + '' for s in Link_List] # append ';D2', makes sure we have a non-terminal group (which would be displayed incorrectly)

    # add carbons LAST. Special case!
     
    Link_List.insert(0,'C&h2') # methylene, -CH2-
    Link_List.insert(0,'C&h3&D1') # methyl, -CH3

    BCI=[0,1] #BCI=bonus carbon index; keeps track of when we need to add another carbon to the label:

    # nested loops construct the label
    for i in range(0,len(Link_List)): 
        for q in range(0,len(Link_List)):  
            for n in range(2,20,1): # this defines to what length chain we'll abbreviate.
                
                #our label is made on the fly. 
                # we need two; one for the 'left' and one for the 'right', depending how it's drawn...
                label1 = ((Link_List[i].replace(';D2','').replace(';h1','H')*(i not in BCI)) +  
                        'C<sub>' + 
                        str(n+(i in BCI) + (q in BCI)) + 
                        '</sub>H<sub>' +
                        str((2*(n+(i in BCI)+(q in BCI)))+(q == BCI[0])+(i == BCI[0]))  + 
                        '</sub>'
                        +(Link_List[q].replace(';D2','').replace(';h1','H')*(q not in BCI)))
                
                label2 = ((Link_List[q].replace(';D2','').replace(';h1','H')[::-1]*(q not in BCI)) +  
                        'C<sub>' + 
                        str(n+(i in BCI) + (q in BCI)) + 
                        '</sub>H<sub>' +
                        str((2*(n+(i in BCI)+(q in BCI)))+(q == BCI[0])+(i == BCI[0]))  + 
                        '</sub>'
                        +(Link_List[i].replace(';D2','').replace(';h1','H')*(i not in BCI)))
                #for label2 just flip [i] and [q] arround
                
                replace_dict = {
                                "[":"", # remove brackets for alkene
                                "]":"", # remove brackets for alkene
                                "#":"#", # replace hash with ... hash, because ≡ doesn't work.
                                "2h;":"H<sub>2</sub>",
                                ";h2":"H<sub>2</sub>"
                               }

                for symb in replace_dict:
                    label1 = label1.replace(symb,replace_dict[symb])
                    label2 = label2.replace(symb,replace_dict[symb])

                sma = ('label' + '\t' +  # sma is our smiles/smarts lookup list; for a given defn we have a matching smiles string.
                       '['+Link_List[i] + ']'+
                       'C'*n +  
                       '['+Link_List[q] + 
                       ']' + '\t' + 
                       label1 + '\t' + 
                       label2 + '\n')

                defns.append(sma)
                
    # misc functional groups
    # these are all hard coded, but common, functional groups.
    # hard-coding is not really ideal, but there is no obvious way to do it (names are not systematic, Bn, Bz etc)

    #defn2a is our definitions for the structure when group is on the right hand side...
    defn2a = ['CN',
            'NO<sub>2</sub>',
            'SF<sub>5</sub>',
            'CF<sub>3</sub>',
            'NCS',
            
            # reactive groups
            'OTf',
            'NTf<sub>2</sub>',
            'OTs',
            'B(pin)',
            'B(mida)',
            'B(OH)<sub>2</sub>',
              
             # protecting groups
             'TMS',
             'TBDMS',
             'TBDPS',
             'TIPS',
             'OTMS',
             'OTBDMS',
             'OTBDPS',
             'OTIPS',
             'Bn',
             'Bz',
             'Trt',
              
             'OAc',
             'O<sup>t</sup>Bu',
             'O<sup>i</sup>Pr,'
            ]

    # defn2b is definitions for our group whne on the left hand side...
    defn2b = ['NC',
            'O<sub>2</sub>N',
            'F<sub>5</sub>S',
            'F<sub>3</sub>C',
            'SCN',
            
            # reactive groups
            'TfO',
            'Tf<sub>2</sub>N',
            'TsO',
            'B(pin)',
            'B(mida)',
            '(HO)<sub>2</sub>B',
             
             # protecting groups
             'TMS',
             'TBDMS',
             'TBDPS',
             'TIPS',
             'TMSO',
             'TBDMSO',
             'TBDPSO',
             'TIPSO',
             'Bn',
             'Bz',
             'Trt',
              
             'AcO',
             '<sup>t</sup>BuO',
             '<sup>i</sup>PrO,'
            ]

    sma2 = ['*C#N',
            '*[N+](=O)[O-]',
            'S(F)(F)(F)(F)(F)',
            'C(F)(F)(F)',
            'N=C=S',
            
            #reactive groups
            '*OS(=O)(=O)C(F)(F)(F)',                        # triflate ester
            '*N(S(=O)(=O)C(F)(F)(F))(S(=O)(=O)C(F)(F)(F))', #triflate amide?? wording?
            '*OS(=O)(=O)c0ccc(C)cc0',                       #tosylate
            '*B0OC(C)(C)C(C)(C)O0',                         #pinacol boronate
            '*B0OC(=O)CN(C)CC(=O)O0',                       #mida boronate
            '*B(O)(O)',                                     #boronic acid
            
            # protecting grps
            '*[Si](C)(C)(C)',
            '*[Si](C)(C)C(C)(C)(C)',
            '*[Si](c1ccccc1)(c2ccccc2)C(C)(C)(C)',
            '*[Si](C(C)C)(C(C)C)(C(C)C)',
            '*O[Si](C)(C)(C)',
            '*O[Si](C)(C)C(C)(C)(C)',
            '*O[Si](c1ccccc1)(c2ccccc2)C(C)(C)(C)',
            '*O[Si](C(C)C)(C(C)C)(C(C)C)',
            '*Cc1ccccc1',
            '*OC(=O)c1ccccc1',
            '*C(c1ccccc1)(c2ccccc2)(c3ccccc3)',
            
            '*OC(=O)C',                                      #acetoxy
            '*OC(C)(C)C',                                    #tbutoxy
            '*OC(C)C',
           ]

    for n in range(len(defn2a)):
        defns.append(defn2a[n] + '\t' + sma2[n] + '\t' + defn2a[n] + '\t' + defn2b[n] +  '\n')
    #print(defns)
    abbrev_defns = '\t'.join([x for x in defns[::-1]])

    abbrevs = rdAbbreviations.ParseAbbreviations(abbrev_defns)
    return (abbrevs)
    

def align_mol(MolObject):
    '''
    Function for 'aligning' a molecule when drawn by RDkit.
    we have a list of templates that we look for. These are SMARTS strings, e.g.
    https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
    
    We look for a match in our structure, and align to this.
    
    The Templates are listed in terms of decreasing priority, i.e. [0] is the highest
    priority.
    
    We basically want the molecule to have its length parallel to the tape length; so
    mostly ring-to-ring bonds are good for alignment. Might be an argument for some
    odd-heterocycles to go at the top (e.g. oxadiazole).
    '''
    Templates = ['c-c',       # aromatic-to-aromatic bond (E.g. biphenyl)
                'c-[!c;R]',   # aromatic to cycic hydrocarbon
                'c-C',       # aromatic-to-alicyclic
                '[R]-[R]',    # any ring-to-ring bond
                '[C;R]!@[C;R]', # bond between saturated rings.
                '[c&H0aac&H0]',   # substituted 1,4-carbons in a benzene ring
                '[*;h1]~[C;R]@[C;R]@[C;R]@[C;R]',   # substituted 1,4-carbons in a saturated ring
                '*#*',  # any tripple bond
                ]

    for template in Templates[::-1]:
        tplt = Chem.MolFromSmarts(template)
        if MolObject.GetSubstructMatches(tplt):
            #print('match found for :' +  template)
            AllChem.Compute2DCoords(tplt) 
            AllChem.GenerateDepictionMatching2DStructure(MolObject,tplt)

    return(MolObject)
    
def crop_img(filename):
    '''
    Code for automatically removing whitepsace from image.
    
    Credit:
    https://stackoverflow.com/questions/10615901/trim-whitespace-using-pil
    '''
    im = Image.open(filename)
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)

def draw_mol(mol):
    '''
    Draws a structure
    
    aligns to "template", i.e. the first atom - makes it nice and linear
        
    Alignment is handled by the align_mol function (above!)
    
    Args:
    Smiles - smiles string to draw_mol
    
    Returns:
    img - a PIL image of the molecule in questeion
    '''
    
    #abbrevs = rdAbbreviations.GetDefaultAbbreviations()        # use the standard rdkit dict of abbreviations
    
    abbrevs = abbrev_mol()
    mol = rdAbbreviations.CondenseMolAbbreviations(mol,abbrevs,maxCoverage=1)
    
    ps = rdCoordGen.CoordGenParams()
    psminimizerPrecision = ps.sketcherBestPrecision
    rdCoordGen.AddCoords(mol,ps)
    
    mol = align_mol(mol)
    img = Draw.MolToImage(mol)
    img.save("temp_mol_img.png")
    return(img)
