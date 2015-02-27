#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DESCRIPTION
This software can be used to find cluster in biological structures (PDB and
Dynamics trajectories like Gromacs and Amber). 

SYNTAXE : 
/!\ NEVER USE SPACES
The Basic Syntaxe is => CHAIN:RES_NUM:ATOMS_NAME
chain, amino acids and atoms are separed with ":"
You can select an amino acid interval with RES_NUM-RES_NUM 
    => 1-30
You can select several atoms type with "-" to (it's apply on the selected AAs)
    => CA-CB
If you wan to select all elements of a group : use "." or nothing.
    => A:.:. (means all atoms of all residues of the chain A)
    => A:: is the same as A:.:.
You can have several selection if you put a "," between them
    => A:.:.,B:1-30:CA (the Entire chain A and only the carbon alpha of the firth
        30 aa of the chain B)
exemples : 
:20-30:CA-CB => Alpha and beta carbon of the amino acids 20 to 30 of all chains
A:1-100:.,B:200-300:CA => alpha carbons of aa 200 to 300 of the chain B AND all
                          atoms of aa 1 to 100 of the chain A.
    
ARGUMENTS : 
    -f : a list of pdb (like "/structure/*.pdb") 
    -l : list of pdb code! /!\ IF YOU USED THAT : -f as to be the PATH of your pdb files
    -s : Selection string
    -a : Alignement type (g = Global, l = Local (on selection string))
    -o : output name for the log file (default is cluster.log)
    
USAGE : 
python Cluster_Analysis.py -f *.pdb -s [SELECTION] -a Alignement_type -o output_log
"""

__author__ = "Thibault TUBIANA"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "20../.."
#==============================================================================
#                     MODULES                            
#==============================================================================
import argparse
from Bio.PDB.PDBParser import PDBParser
import extract_coord as EC
import os,sys
import progressbar as pg
import Bio.PDB.Superimposer as Superimposer
from Bio.SVDSuperimposer import SVDSuperimposer 
from distmat import DistanceMatrix
from nmrclust import NMRClust
#==============================================================================
#                     GLOBAL VARIABLES                            
#==============================================================================
WIDGETS = [pg.Bar('>'), ' ', pg.ETA(), ' ', pg.ReverseBar('<')]
         

#==============================================================================
#                     TOOL FONCTION
#==============================================================================

class NullDevice:
    """
    Silencer to remove stdout for one fonction.
    """
    def write(self, s):
        pass
    def flush(self):
        pass

def silence():
    """
    Returns a decorator for discarding output to the givens stream for the
    execution of the decorated metod.
    """
    def dec(fun):

        def wrapper(*args, **kwargs):
            sys.stdout.flush()
            sys.stderr.flush()
            backup_out = sys.stdout
            backup_err = sys.stderr
            sys.stdout = NullDevice()
            sys.stderr = NullDevice()
            try:
                result = fun(*args, **kwargs)
            except Exception, e:
                raise e
            finally:
                sys.stdout = backup_out
                sys.stderr = backup_err
            return result
        wrapper.__doc__ = fun.__doc__
        wrapper.__name__ = fun.__name__
        return wrapper
    return dec   
    
    
def def_pathfiles(myArgs):
    """
    This function return a list of pdbpath from a list passed in arguments
    or from a text file.
    """    
    listpdb=myArgs['list'] #list of pdb or path to the pdb folder
    
    if(listpdb==None):
        pathFiles=myArgs['files'][0] #warning : list of list.
    #Otherwise if we used a list of pdb code, we take the path to the pdb files,
    #and we add the pdb code
    else:
        pathFiles=[]
        path_cat=myArgs['files'][0][0]
        file_list=open(listpdb,'r')
        #checking the path
        if not (path_cat[-1]=='/'):
            path_cat=path_cat+'/'
            #check if the directory exist
            if not (os.path.exists(path_cat)):
                print("ERROR : the pdb directory that you give is not exist")
                print("    -> please check the directory")
                exit(1)
        for line in file_list:
            #we add th    def create_sphere(self):

#e path for each pdb file
            pathFiles.append(line.lower().split('\n')[0]+'.pdb')
        file_list.close()
    return pathFiles   


    
@silence()
def load_pdb(pdbCode, pdb_file):
    par = PDBParser()
    structure=par.get_structure(pdbCode, pdb_file)
    return(structure)
    

            
            
        
        
#==============================================================================
#                     FONCTIONS
#==============================================================================


def parseArg():
    """
    This fonction will the list of pdb files and the distance
    @return: dictionnary of arguments
    Ex : 
    python Cluster_Analysis.py -f *.pdb -s A:1-30:CA
    """
    arguments=argparse.ArgumentParser(description="\
            Description\
            ")
    arguments.add_argument('-f', "--files", help="list of pdf files or path to the pdb folder (if\
                            you use '-l' argument", required=True, action="append", nargs="*")
    arguments.add_argument('-l','--list', help="OPTIONAL : list of pdb in a file. \
            WARNING! -f argument have to be THE PATH to your pdb files", default=None)
    arguments.add_argument('-s','--select', help="selection syntaxe", required=True)
    arguments.add_argument('-a','--alignement', help="Alignement type : g=Global, l=Local", default="g")
    arguments.add_argument('-o','--output', help="Output File", default="cluster.log")
    args = vars(arguments.parse_args())
    return(args)
    
    
    
    
    
    
def uncode_selection(selection):
    """
    This function will resturn an array of selected atoms
    Syntaxe --> ":" => Substructure; "," => other selection group; "." ==> All;
                "()" => atoms selection of residues (name not numbers); ";" atoms separator
    chain:res-res:atoms
    ex:
    A:1-3:CA-CB ==> atoms CA and CB of res 1,2 and 3 of chain A
    """
    chain=[]
    res=[]
    atoms=[]
    allstruct=False
    for sel in selection.split(','): #for several selection group
        if sel==".":
            chain.append("")
            res.append("")
            atoms.append("")
            allstruct=True
        if not allstruct:
            
            #CHAIN
            chain_sel=sel.split(":")[0]
            if not chain_sel==".":
                chain.append(chain_sel)
            else:
                chain.append("")
            
            #RESIDUES
            try:
                residues_sel=sel.split(":")[1]
                if not residues_sel == "." and not residues_sel=="":
                    if len(residues_sel.split("-")) == 1:
                        res.append([int(residues_sel)])
                    else:
                        res.append(range(int(residues_sel.split("-")[0]),\
                                         int(residues_sel.split("-")[1])+1))
                else:
                    res.append("")
            except:
                res.append("")
            
            #ATOMS
            try: 
                atoms_sel=sel.split(":")[2]
                if not atoms_sel == "." and not atoms_sel=="":
                    atoms.append(sel.split(":")[2].split("-"))
                else:
                    atoms.append("") 
            except:
                atoms.append("")
    return (chain,res,atoms)
                


def alignement(s1,s2):
    """
    This function is used to align two structure (S2 is aligned on S1).
    the new atoms coordinates are injected directly on the structure 2.
    """    
    
    fixed=[x for x in s1.get_atoms()]
    moving=[x for x in s2.get_atoms()]
    sup=Superimposer()
    sup.set_atoms(fixed,moving)
    sup.apply(moving)
    

def create_DM(pdb_atoms_list, alignement_type):
    """
    This function is used to clusterize all structures.
    """
    svd=SVDSuperimposer()
    size=len(pdb_atoms_list)
    FullDM=DistanceMatrix(size)
    pbar = pg.ProgressBar(widgets=WIDGETS, maxval=(size*(size-1)/2)).start()
    counter=0
    for i,frame1 in enumerate(pdb_atoms_list):
        for j, frame2 in enumerate(pdb_atoms_list[i+1:]):
            svd.set(frame1,frame2)
            svd.run()
            rms=svd.get_rms()           #RMS with local alignement
            rms_raw=svd.get_init_rms()  #RMS with the GLOBAL alignement made before
            if alignement_type == "l":
                FullDM.set(i,i+j+1, rms)
            else:
                FullDM.set(i,i+j+1, rms_raw)
            pbar.update(counter)
            counter+=1
    pbar.finish()
    return FullDM
    
    

def wrote_logfile(outpath, cl, PDB_name_list):
    out=open(outpath, "w")
    for i,c in enumerate(cl.clusters):
        out.write("cluster %i\n" %(i+1))
        out.write("    size = %i ; representative frame=%i\n" %(c.size, cl.representative(c)+1))
        members=""
        for m in c.members():
            members=members+str(m+1)+","
        out.write("    ID Members : %s\n" %members)
        out.write("    Spread : %.3f\n" %c.spread)
        out.write("    PDB files:\n")
        for m in c.members():
            out.write("       %s\n" %PDB_name_list[m])
        out.write("\n")
    out.close()
    
###############################################################################
#####                               MAIN                                 ######
###############################################################################

if __name__ == "__main__":
    print "********************************************************"
    print "**********  3D STRUCTURES CLUSTERING 0.9  **************"
    print "********************************************************"
    print ""
    print "Please note :"
    print "  *Gromacs trajectories reader has not been implemented yet"
    print "  *Amber trajectories reader has not been implemented yet"
    print "Run the software with : python Cluster_Analysis.py -f *.pdb -s [SELECTION]\n\n"
    #We get all arguments
    myArgs=parseArg()
    
    pathFiles=def_pathfiles(myArgs)
    selection=uncode_selection(myArgs["select"])
    alignement_type=myArgs["alignement"]    
    outpath=myArgs["output"]    
    
    if alignement_type not in ["g","l"]:
        print "%s is not a recognized parameter for the alignement type"
        print "I will use a global alignement instead!"
        alignement_type="g"
    
    
    
    pdb_atoms_list=range(len(pathFiles))
    PDB_name_list=range(len(pathFiles))
    GAC=EC.Get_Array_Coord()

    structure_list=range(len(pathFiles))

    pbar = pg.ProgressBar(widgets=WIDGETS, maxval=len(pathFiles)).start()
    print "======= PDB loading ======="
    for i,pdb_file in enumerate(pathFiles):
        pdbCode=os.path.split(pdb_file)[1].split('.')[0]
        PDB_name_list[i]=os.path.split(pdb_file)[1]
        structure_list[i]=load_pdb(pdbCode, pdb_file)
        pbar.update(i)
    pbar.finish()
    print "         >done"
    
    
    print "====== GLOBAL PDB Alignement ======"
    print "   >all pdb will be alignement on the first one"
    pbar = pg.ProgressBar(widgets=WIDGETS, maxval=len(structure_list)).start()
    s1=structure_list[0]
    for i in xrange(1,len(structure_list)):
        alignement(s1,structure_list[i])
        pbar.update(i-1)
    pbar.finish()
    print "         >done"    
    
    
    
    pbar = pg.ProgressBar(widgets=WIDGETS, maxval=len(structure_list)).start()
    print "======= COORDINATES EXTRATION AND CONVERSION ======="
    for i,structure in enumerate(structure_list):
        pdb_atoms_list[i]=GAC.get_array(selection,structure)
        pbar.update(i)
    pbar.finish()
    print "         >done"        
        
        
    print "======= DISTANCE MATRIX CALCULATION ========"    
    dm=create_DM(pdb_atoms_list, alignement_type)
    print "         >done"
    
    
    print "====== Clustering ========"
    cl=NMRClust(dm)
    print "         >done"        
    
    print "**** Cluster Results"
    for i,c in enumerate(cl.clusters):
          print "cluster %i" %(i+1)
          print "    size = %i ; representative frame=%i" %(c.size, cl.representative(c)+1)
          print "    Members : ", [x+1 for x in c.members()]
          print "    Spread :", c.spread
    
        
    wrote_logfile(outpath, cl,PDB_name_list)   
        
        
        
        
        
        
        
        
