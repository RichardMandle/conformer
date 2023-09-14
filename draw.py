import subprocess
from rdkit import Chem
from rdkit.Chem import PyMol

def pymol_draw(mol, path, name='NoName', style=True, ray=False, grid=True, start=0, end=-1):

    """
    Function for calling pymol and making a nice image
    having style = 'yes' will give a standard image in the UoL-SMP style
    raytacing defaults to off (ray='no') but can be called ('yes') or run externally in PyMol
    
    Make sure you update the pymol path in the os.system call
    
    you can specify start = i and end = m, where i and m are integers, to
    output only conformers between the ith and mth to Pymol. Python
    indexes from zero (i.e. hte first conformer is no 0, the 2nd is 1 etc.)
    
    This function is painfully slow, is there a better way?
    """
    
    #first up, launch PyMol with the -R flag so we can pass data to it:
    subprocess.call(path + ' -R')
    print('\nSending data to PyMol and drawing; this can be slow!')
    if end == -1:
        end = mol.GetNumConformers()
        
    if path is None or path == '':
        print('Error: PyMol path not set!')
        return
        
    v= PyMol.MolViewer()
    v.DeleteAll()
    
    for i in range(start,end):
        v.ShowMol(mol, confId=i,name=name+'_conf-%d'%i,showOnly=False)
    
    if style:
        # call grid mode; coloring options
        # putting multiple v.server.do calls on one line speeds it up a lot.
        v.server.do('set_color oxygen, [1.0,0.4,0.4];'+
                    'set_color nitrogen, [0.5,0.5,1.0];'+
                    'set_color carbon, [0.5,0.5,0.5];'+
                    'set_color hydrogen, [1.0,1.0,1.0];'+
                    'util.cbag;recolor;'+
                    'remove solvent;'+
                    'as spheres;'+
                    'bg white;'+
                    'set light_count,8;'+
                    'set spec_count,1;'+
                    'set shininess, 10;'+
                    'set specular, 0.25;'+
                    'set ambient,0;'+
                    'set direct,0;'+
                    'set reflect,1.5;'+
                    'set ray_shadow_decay_factor, 0.1;'+
                    'set ray_shadow_decay_range, 2;'+
                    'unset depth_cue;' + 
                    ('set grid_mode, on;' * grid))

        if ray:
            print('I am going to raytrace')
            v.server.do('ray') # you can raytrace if you want

    return() 