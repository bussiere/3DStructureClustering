import numpy as np

class Get_Array_Coord:

    

    def add_coord(self, c,r,a,atoms_coord, structure):
        try:
            atoms_coord.append(structure[0][c][r][a].coord)
        except:
            print c,r,a
            print structure[0][c][r][a].coord
            print "the Syntaxe selection has an on more error"
            print "please check again your selection string or your pdb files"
            exit(1)


    def add_res(self, c,res,atoms,atoms_coord, structure):
        if res == "":
            self.get_all_res(c,atoms,atoms_coord, structure)
        else: 
            for r in res:
                if atoms == "":
                    self.get_all_atoms(c,r,atoms_coord, structure)
                else:
                    for a in atoms:
                        self.add_coord(c,r,a,atoms_coord, structure)



    def get_all_chains(self, c, res, atoms,atoms_coord, structure):
        for chain in structure[0].get_list():
            c=chain.id
            self.add_res(c,res,atoms,atoms_coord, structure)



    def get_all_res(self, c,atoms, atoms_coord, structure):
        for res in structure[0][c].get_list():
            r = res.get_id()[1]
            if atoms == "":
                self.get_all_atoms(c,r,atoms_coord, structure)
            else:
                for a in atoms:
                    self.add_coord(c,r,a,atoms_coord, structure)


    def get_all_atoms(self, c,r, atoms_coord, structure):
        for atoms in structure[0][c][r].get_list():
            atoms_coord.append(atoms.coord)


        
    def get_array(self, selection,structure):
        atoms_coord=[]
        for i in xrange(len(selection[0])):
            c=selection[0][i]
            res=selection[1][i]
            atoms=selection[2][i]
            if c == "":
                self.get_all_chains(c,res,atoms, atoms_coord, structure)
            else:
                self.add_res(c,res,atoms,atoms_coord, structure)
        return np.array(atoms_coord)
    
    
if __name__ == "__main__":
    selection=(['A'], [[1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060]], [['CA']])
    try:
        from Bio.PDB.PDBParser import PDBParser
        par = PDBParser()
        structure=par.get_structure("1IHM", "1ihm.pdb")
    except:
        print "please download 1IHM.pdb in your current folder"
        exit(1)
    GAC=Get_Array_Coord()
    atoms_list_array=GAC.get_array(selection,structure)
    print atoms_list_array
    
