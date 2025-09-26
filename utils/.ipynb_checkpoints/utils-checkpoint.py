import gzip
import numpy as np

class PDB(object):
    '''
    Class object to load PDB structure files of pdb or cif format.

    Parameters
    ----------
    filename : str
        Path to file.

    ignore_waters : bool
        To ignore or not the waters in the file, default False.

    ignore_other_HETATM : bool
        To ignore or not the HETATMs in the file, excluding waters, default False.

    multi_models : bool
        To include in the data or not multiple models of the structure in the file, default False.

    Attributes
    ----------
    natoms : int
        Number of protein atoms of all the chains in the structure.

    water_coords : numpy array
        Coordinates of water atoms, note that here we only keep oxygen atoms, shape (W,3).

    wb : numpy array
        B factor of water coordinates, shape (W).

    HETATM_coords : numpy array
        Coordinates of HETATM atoms, excluding water molecules, shape (H,3).

    HETATM_name : numpy array
        Names of HETATM atoms, excluding water molecules, shape (H).

    HETATM_num : numpy array
        Number of HETATM instance, excluding water molecules, shape (H).

    HETATM_atomnum : numpy array
        Atom number of HETATM atoms, excluding water molecules, shape (H).

    HETATM_atomtype : numpy array
        Atom type of HETATM atoms, excluding water molecules, shape (H).

    SSE : numpy array
        Empty numpy array ready to keep SSE info, shape (N).

    resolution : float
        Resolution of structure.

    SSEraw : numpy array
        Info of ranges for SSE info, shape (S,3)

    atomnum : numpy array
        Atom number of protein atoms, shape (N).

    atomname : numpy array
        Atom name of protein atoms, shape (N).

    atomalt : numpy array
        Alternative atoms for protein atoms, shape (N).

    resname : numpy array
        Residue name that protein atom belongs, shape (N).

    atomtype : numpy array
        Atom type of protein atoms, shape (N).

    resnum : numpy array
        Residue number of protein residues, shape (N).

    resalt : numpy array
        Alternative residues for protein residues, shape (N).

    chain : numpy array
        Chain that atom belongs, shape (N).

    coords : numpy array
        Atom coordinates of protein atoms, shape (N,3).

    occupancy : numpy array
        Occupancy of protein atoms, shape (N).

    b : numpy array
        B factor of protein atoms, shape (N).

    self.cella : float
        Unit cell parameters, length of a axis.

    self.cellb : float
        Unit cell parameters, length of b axis.

    self.cellc : float
        Unit cell parameters, length of c axis.

    self.cellalpha : float
        Unit cell parameters, angle of a axis.

    self.cellbeta : float
        Unit cell parameters, angle of b axis.

    self.cellgamma : float
        Unit cell parameters, angle of c axis.
    '''
    def __init__(self, filename, ignore_waters=False, ignore_other_HETATM=False, multi_models=False):
        #Define lists and variables
        self.natoms, self.water_coords, self.wb = 0, [], []
        self.HETATM_coords, self.HETATM_name, self.HETATM_chain, self.HETATM_num, self.HETATM_atomnum, self.HETATM_atomtype = [],[],[],[],[],[]
        self.SSE, self.resolution, self.SSEraw = [], None, []
        self.atomnum, self.atomname, self.atomalt, self.resname, self.atomtype, self.resnum, self.resalt, self.chain, self.coords = [],[],[],[],[],[],[],[],[]
        self.occupancy, self.b  = [],[]
        #Check if file is in cif format
        cif = False
        if filename.split('.')[-1] == 'cif' or filename.split('.')[-2] == 'cif':
            cif = True
        #If file is compressed uncompress it
        if filename.split('.')[-1] == 'gz':
            with gzip.open(filename, 'rt', encoding='utf-8') as file:
                f = file.readlines()
        #Else just open it
        else:
            with open(filename, 'r') as file:
                f = file.readlines()
        #If the format is pdb
        if not cif:
            #Start reading the PDB file
            for i, line in enumerate(f):
                sline = line.split()
                #If line starts with 'ATOM'
                if line[:4]=='ATOM':
                    #Keep protein's atoms info.
                    self.atomnum.append(int(float(sline[1])))
                    self.atomname.append(line[12:16].strip())
                    self.atomalt.append(line[16])
                    self.resname.append(line[17:21].strip())
                    atomtype = sline[-1]
                    if len(atomtype)-1:
                        try:
                            int(atomtype[1])
                            atomtype = atomtype[0]
                        except:
                            atomtype = atomtype[0].upper() + atomtype[1].lower()
                    self.atomtype.append(atomtype)
                    self.resnum.append(int(float(line[22:26])))
                    self.resalt.append(line[26])
                    self.chain.append(line[21])
                    self.coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    self.occupancy.append(float(line[56:60]))
                    self.b.append(float(line[60:66]))
                    #self.charge[atom] = line[78:80].strip('\n')
                    #self.nelectrons[atom] = electrons.get(self.atomtype[atom].upper(),6)
                    self.natoms += 1
                    continue
                #If line starts with 'HETATM'
                if line[:6] == 'HETATM':
                    #Keep waters
                    if not ignore_waters and line[13] == 'O' and ((line[17:20]=='HOH') or (line[17:20]=='TIP') or (line[17:20]=='WAT')):
                        self.water_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        self.wb.append(float(line[60:66]))
                        continue
                    #Keep the rest hetatm
                    if not ignore_other_HETATM and sline[-1][0]!='H':
                        self.HETATM_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        self.HETATM_name.append(line[17:21].strip())
                        self.HETATM_chain.append(line[21])
                        self.HETATM_num.append(int(float(line[22:26])))
                        self.HETATM_atomnum.append(int(float(line[6:12])))
                        atomtype = sline[-1]
                        if len(atomtype)-1:
                            try:
                                int(atomtype[1])
                                atomtype = atomtype[0]
                            except:
                                atomtype = atomtype[0].upper() + atomtype[1].lower()
                        self.HETATM_atomtype.append(atomtype)
                        continue
                #Here you can add more self objects from Header.
                #Resolution info in pdb format files.
                if line[:22] == 'REMARK   2 RESOLUTION.':
                    self.resolution = sline[3]
                    continue
                #Helix info in pdb format files.
                if line[:5] == 'HELIX':
                    self.SSEraw.append([range(int(line[20:25]),int(line[32:37])+1), line[19], 'H'])
                    continue
                #Beta sheet info in pdb format files.
                if line[:5] == 'SHEET':
                    self.SSEraw.append([range(int(line[22:26]),int(line[33:37])+1), line[21], 'S'])
                    continue
                #Crystal info in pdb format files.
                if line[:6] == 'CRYST1':
                    self.cella, self.cellb, self.cellc = float(sline[1]), float(sline[2]), float(sline[3])
                    self.cellalpha, self.cellbeta, self.cellgamma = float(sline[4]), float(sline[5]), float(sline[6])
                    continue
                #Return structure only of Model 1 if there are multiple Models
                if sline[0] == 'MODEL' and not cif and not multi_models and sline[1] != '1':
                    break
        #If the format is cif
        else:
            b_strands = False
            #Start reading the PDB file
            for i, line in enumerate(f):
                sline = line.split()
                #If line starts with 'ATOM'
                if line[:4]=='ATOM':
                    #First check if you have multiple models and you dont want them
                    if int(sline[20])-1 and not multi_models:
                        break
                    self.atomnum.append(int(float(sline[1])))
                    self.atomname.append(sline[3])
                    self.atomalt.append(sline[4])
                    self.resname.append(sline[5])
                    atomtype = sline[2]
                    if len(atomtype) == 2:
                        atomtype = atomtype[0].upper() + atomtype[1].lower()
                    self.atomtype.append(atomtype)
                    self.resnum.append(int(float(sline[16])))
                    self.chain.append(sline[6])
                    self.coords.append([float(sline[10]), float(sline[11]), float(sline[12])])
                    self.occupancy.append(float(sline[13]))
                    self.b.append(float(sline[14]))
                    self.natoms += 1
                    continue
                #If line starts with 'HETATM'
                if line[:6] == 'HETATM':
                    if not ignore_waters and sline[2] == 'O' and ((sline[5]=='HOH') or (sline[5]=='TIP') or (sline[5]=='WAT')):
                        #Keep coordinates of water Oxygens.
                        self.water_coords.append([float(sline[10]), float(sline[11]), float(sline[12])])
                        self.wb.append(float(sline[14]))
                        continue
                    if not ignore_other_HETATM and sline[2]!='H':
                        self.HETATM_coords.append([float(sline[10]), float(sline[11]), float(sline[12])])
                        self.HETATM_name.append(sline[5])
                        self.HETATM_chain.append(sline[7])
                        self.HETATM_num.append(int(float(sline[16])))
                        self.HETATM_atomnum.append(int(float(sline[1])))
                        self.HETATM_atomtype.append(sline[2])
                        continue
                #Here you can add more self objects from Header.
                #Resolution info in cif format files.
                if line[:25] == '_reflns.d_resolution_high' or line[:33] == '_em_3d_reconstruction.resolution ':
                    self.resolution = sline[1]
                    continue
                #Helix info in cif format files.
                if line[:6] == 'HELX_P':
                    self.SSEraw.append([range(int(sline[-7]),int(sline[-4])+1), sline[-8], 'H'])
                    continue
                #Beta sheet info in cif format files.
                if line[:35] == '_struct_sheet_range.end_auth_seq_id':
                    b_strands = True
                    continue
                if b_strands:
                    if line[:1] == '#':
                        b_strands=False
                        continue
                    self.SSEraw.append([range(int(sline[-4]),int(sline[-1])+1), sline[-2], 'S'])
                    continue
                #Crystal info in cif format files.
                if line[:5] == '_cell':
                    if sline[0] == '_cell.length_a':
                        self.cella = float(sline[1])
                    if sline[0] == '_cell.length_b':
                        self.cellb = float(sline[1])
                    if sline[0] == '_cell.length_c':
                        self.cellc = float(sline[1])
                    if sline[0] == '_cell.angle_alpha':
                        self.cellalpha = float(sline[1])
                    if sline[0] == '_cell.angle_beta':
                        self.cellbeta = float(sline[1])
                    if sline[0] == '_cell.angle_gamma':
                        self.cellgamma = float(sline[1])
                    continue
        #Return structure when every atom and water atom have been searched.
        try:
            self.resolution = float(self.resolution)
        except:
            pass
        #Prepare SecondaryStructureElements info
        self.SSE = np.zeros((self.natoms), dtype=np.dtype((str,1)))
        self.SSEraw = np.array(self.SSEraw, dtype=object)
        #Turn every list to numpy array
        for attribute, value in vars(self).items():
            if type(value)==list:
                setattr(self, attribute, np.array(value))
        return

    def remove_waters(self):
        idx = np.where((self.resname=="HOH") | (self.resname=="TIP"))
        self.remove_atoms_from_object(idx)

    def remove_by_atomtype(self, atomtype):
        idx = np.where((self.atomtype==atomtype))
        self.remove_atoms_from_object(idx)

    def remove_by_atomname(self, atomname):
        idx = np.where((self.atomname==atomname))
        self.remove_atoms_from_object(idx)

    def remove_by_atomnum(self, atomnum):
        idx = np.where((self.atomnum==atomnum))
        self.remove_atoms_from_object(idx)

    def remove_by_resname(self, resname):
        idx = np.where((self.resname==resname))
        self.remove_atoms_from_object(idx)

    def remove_by_resnum(self, resnum):
        idx = np.where((self.resnum==resnum))
        self.remove_atoms_from_object(idx)

    def remove_by_chain(self, chain):
        idx = np.where((self.chain==chain))
        self.remove_atoms_from_object(idx)

    def remove_atoms_from_object(self, idx):
        mask = np.ones(self.natoms, dtype=bool)
        mask[idx] = False
        self.atomnum = self.atomnum[mask]
        self.atomname = self.atomname[mask]
        self.atomalt = self.atomalt[mask]
        self.resalt = self.resalt[mask]
        self.resname = self.resname[mask]
        self.resnum = self.resnum[mask]
        self.chain = self.chain[mask]
        self.coords = self.coords[mask]
        self.occupancy = self.occupancy[mask]
        self.b = self.b[mask]
        self.atomtype = self.atomtype[mask]
        self.SSE = self.SSE[mask]
        #self.charge = self.charge[mask]
        #self.nelectrons = self.nelectrons[mask]
        self.natoms = len(self.atomnum)

    def rearrange_resalt(self):
        #Keep indexes where you have added residues
        ind = np.where(self.resalt!=' ')[0]
        if not len(ind):
            return
        #Find chains of these added residues
        resalt_ch = self.chain[ind]
        #Find unique chains
        diff_chains = np.unique(resalt_ch)
        #For each unique chain id
        for ch in diff_chains:
            #Keep indexes of added residues only for your chain
            ind_resalt_ch = ind[np.where(resalt_ch==ch)]
            #Keep last index of your chain's residues
            last_ch_ind = np.where(self.chain==ch)[0][-1]+1
            #For each atom in chain in added residue
            for atom in ind_resalt_ch:
                #Find where added residue starts
                if self.resalt[atom-1]!=self.resalt[atom] and self.resnum[atom-1] == self.resnum[atom]:
                    #Add to all consquent residue +1 number
                    self.resnum[atom:last_ch_ind] += 1
                #If added residues are before residue number, change the numbering of the original residue +1
                if atom+1+1>self.natoms:
                    continue
                if self.resalt[atom+1] == ' ' and self.resnum[atom+1] == self.resnum[atom]:
                    self.resnum[atom+1:last_ch_ind] += 1
        #Remove alternative residues indexes
        self.resalt[ind] = ' '

def clean_pdb_alternative_atoms(pdb):
    '''
    Remove alterative atoms in pdb protein coordinates

    Parameters
    ----------
    pdb : class
        object of PDB class

    Returns
    -------
    pdb : class
        cleaned PDB class object from alternative atoms.
    '''
    delete = np.where(np.isin(pdb.atomalt, (' ', 'A', '.'))==False)[0]
    #If you have indexes to delete proceed
    if len(delete):
        #Remove selected atoms if any.
        pdb.remove_atoms_from_object(delete)
    return

def clean_pdb_atomname_and_resname(pdb):
    '''
    Remove non canonical atomnames and resnames.

    Parameters
    ----------
    pdb : class
        object of PDB class

    Returns
    -------
    pdb : class
        cleaned PDB class object from non canonical atomnames and resnames.
    '''
    #Find where you have non canonical atomnames and resnames
    non_canonical_atomnames = np.where(np.isin(pdb.atomname, atomnames)==False)[0]
    non_canonical_resnames = np.where(np.isin(pdb.resname, resnames)==False)[0]
    #If you have indexes
    if len(non_canonical_atomnames) or len(non_canonical_resnames):
        #Compute the union.
        delete = np.union1d(non_canonical_atomnames, non_canonical_resnames)
        #Remove selected atoms if any.
        pdb.remove_atoms_from_object(delete)
    return

def clean_pdb(pdb):
    '''
    Remove alterative atoms in pdb protein coordinates

    Parameters
    ----------
    pdb : class
        object of PDB class

    Returns
    -------
    pdb : class
        cleaned PDB class object from alternative atoms and non canonical atomnames and resnames.
    '''
    alt_atoms = np.where(np.isin(pdb.atomalt, (' ', 'A', '.'))==False)[0]
    #Find where you have non canonical atomnames and resnames
    non_canonical_atomnames = np.where(np.isin(pdb.atomname, atomnames)==False)[0]
    non_canonical_resnames = np.where(np.isin(pdb.resname, resnames)==False)[0]
    #If you have indexes to delete proceed
    if len(non_canonical_atomnames) or len(non_canonical_resnames) or len(alt_atoms):
        #Compute the union.
        delete = np.union1d(alt_atoms, np.union1d(non_canonical_atomnames, non_canonical_resnames))
        #Remove selected atoms if any.
        pdb.remove_atoms_from_object(delete)
    return
