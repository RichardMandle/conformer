from rdkit import Chem
from PIL import Image, ImageChops, ImageDraw
from rdkit.Chem import Draw, AllChem, rdCoordGen

'''
functions pertaining to drawing chemical structures.

Removed functions:
abbreviation engine - this is the same as the one used in ChemLabel; it makes no sense to abbreviate atom groups when
                      we care explicitly about connectivity.
                      
'''

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
    for n in range(1,20):
        Templates.insert(0, '[R]-' + ('[!R]-')*n + '-[R]')

    for template in Templates[::-1]:
        tplt = Chem.MolFromSmarts(template)
        if MolObject.GetSubstructMatches(tplt):
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

    ps = rdCoordGen.CoordGenParams()
    psminimizerPrecision = ps.sketcherBestPrecision
    rdCoordGen.AddCoords(mol,ps)
    
    mol = align_mol(mol)
    img = Draw.MolToImage(mol)
    img.save("temp_mol_img.png")
    return(img)
