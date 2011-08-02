#! /usr/bin/env python


###########################################################
###                                                  	###
###              Build incomplete molecules             ###
###          written by T. Pencheva and D.Lagorce 	###
###             http://www.vls3d.com                	###
###                                                  	###
###########################################################

import os
import os.path
import sys
from string import *
import popen2
import time
import shutil

########################################################

def usage():
    return """Usage of BuildMolecules_sp4.py:
    
##########################################
command line = ./BuildMolecules_sp4.py <ammp.param>
##########################################

<ammp.param> is the paramaters file which needs to be edited with your
             configuration. This file is located in your project directory.

   """

########################################################

def params(fhi_params):

    line = fhi_params.readline()
    path_AMMP = line.split()[1]

    line = fhi_params.readline()
    bank = line.split()[1]

    return path_AMMP,bank
    
#######################################################

def preammp_ligand_textfile_maker(AMMP_path):

    command ='pwd'
    fhin,fhout = popen2.popen2(command)
    for line in fhin:
        working_path = line

    atoms_path = AMMP_path + '/progs/preammp/' 
    preammp_fhin = open('preammp_ligand.txt','w')
    preammp_fhin.write(atoms_path + 'atoms.sp4\n')
    preammp_fhin.write(working_path)
    preammp_fhin.write('input_ligand.pdb\n')
    preammp_fhin.write('input_ligand.ammp\n')
    preammp_fhin.close()

#######################################################

def file_organizer(bank_file,energy_file,choice,warnings):

    output_dir = bank_file[:-5] + '_OUTPUT'
    
    if os.path.exists(output_dir):
        command0 = 'rm -rf ' + output_dir
        os.system(command0)
        
    command1 = 'mkdir ' + bank_file[:-5] + '_OUTPUT'
    os.system(command1)

    command3 = 'rm input.mol2'
    os.system(command3)
    command4 = 'rm input_ligand.pdb'
    os.system(command4)
    command5 = 'rm input_ligand.ammp'
    os.system(command5)
    command6 = 'rm output.mol2' 
    os.system(command6)
    command10 = 'rm mol'
    os.system(command10)

########################################################

def Get_Mol_Amount(fhi):

    amount = 0
    line = fhi.readline()
    while line != '':

        if line[:17] == '@<TRIPOS>MOLECULE':
            amount = amount + 1
            line = fhi.readline()
        else:
            line = fhi.readline()
            
    return amount

########################################################

def GetPos(fh_input_LIG):


    index = {}
    mol = 1
    flag = False
    start = fh_input_LIG.tell()
    stop = 0
    head = fh_input_LIG.tell()
    line = fh_input_LIG.readline()
    while line != '':
        line = line[:-1]

        if line == '@<TRIPOS>MOLECULE':
            if flag == False:
                flag = True
                line = fh_input_LIG.readline()
                mol = line[:-1]

            else:
                flag = False
                fh_input_LIG.seek(head)
                stop = head
                index[mol]=(start,stop)
                start = fh_input_LIG.tell()
                line = fh_input_LIG.readline()
                
        else:
            head = fh_input_LIG.tell()
            line = fh_input_LIG.readline()


    stop = head
    index[mol] = (start,stop)
    print '------>\t'+mol
    return index

########################################################
            
def GetLigand_input(dico_pos,key,bank,input):

    fhi = open(bank)
    fho = open(input,'w')
    start,stop = dico_pos[key]
    fhi.seek(start)
    pos = fhi.tell()
    while pos != stop:            
        line = fhi.readline() 
        fho.write(line)   
        pos = fhi.tell()
    fhi.close()
    fho.close()
        
########################################################

def ligand_convertion2pdb(AMMP_path):

    lig_command = AMMP_path + '/bin/mol2_to_templ_sp4'
    os.system(lig_command)

########################################################

def preammp(AMMP_path):

    preammp_command = AMMP_path + '/bin/preammp <preammp_ligand.txt'
    os.system(preammp_command)
  
########################################################
    
def ammp(AMMP_path):

    ammp_command = AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/build_mol2_dgeom.ammp'
    os.system(ammp_command)

    
########################################################

def ammp2mol2(AMMP_path):

    convertion_command= AMMP_path + '/bin/ammp_to_mol2'
    os.system(convertion_command)


########################################################    

def save_energy(AMMP_path):
    save_energy_command = AMMP_path + '/bin/save_energy'
    os.system(save_energy_command)



########################################################    

    
if __name__ == '__main__':

    debut = time.time()
    
    if len(sys.argv)<2 or len(sys.argv)>2:
        print "\n"
        print "raise IOError"
        print usage()

    else:
        file_params = sys.argv[1] ######### fichier parametres 

        fh_params = open(file_params)
        path_of_AMMP,bank = params(fh_params)       #### return dir de VLS_AMMP fichier protein fichier ligand
        fh_bank = open(bank)
        print '-------> COUNTING THE MOLECULES IN THE BANK'
        print '\n'     
        mol_amount = Get_Mol_Amount(fh_bank)
        print '-------> THE BANK CONTAINS '+ str(mol_amount) + ' MOLECULES'
        new_input_protein = 'input_protein.pdb'
        fh_bank.close()        
        fh_bank = open(bank)
        preammp_ligand_textfile_maker(path_of_AMMP)


        built_bank_file = bank[:-5] + '_Built.mol2'
        badenergy_mols = bank[:-5] + '_BadMolecules.mol2'
        fho_badenergy_mols = open(badenergy_mols,'w')
        energy_bank_file = bank[:-5] + '_energy.txt'
        totalwarning = bank[:-5] + '_total_warnings.txt'
        fho_totalwarning = open(totalwarning,'w')

        fho_built = open(built_bank_file,'w')
        fho_energy_file = open(energy_bank_file,'w')
        top = 'LIGAND : NONBOND :TOTAL ENERGY\n'
        fho_energy_file.write(top)



        i = 0
        while i <= mol_amount-1:
            left = mol_amount - i
            print '\n'
            print '------->\t'+str(left)+' MOLECULES LEFT..................'
            print '\n'
            flag = 0
            id = 0
            input = 'input.mol2'
            fho = open(input,'w')
            pos = fh_bank.tell()
            line = fh_bank.readline()
            while line != '':
            
                if line[:17] == '@<TRIPOS>MOLECULE':

                    if flag == 0:
                        flag = 1
                        id = 1
                        fho.write(line)

                    else:
                        fh_bank.seek(pos)
                        break

                elif id == 1:
                    id = 0
                    ligand_name = line[:-1]
                    fho.write(line)
                                    
                else:
                    fho.write(line)
                    
                pos = fh_bank.tell()
                line = fh_bank.readline()

            fho.close()
            energy_file = 'energy.txt'
            print '------->\t CONVERTION TO PDB'
            ligand_convertion2pdb(path_of_AMMP)

            print '------->\t PREAMMP'
            preammp(path_of_AMMP)

            print '------->\t AMMP'
            ammp(path_of_AMMP)
            print '------->\t AMMP2MOL2'
            ammp2mol2(path_of_AMMP)
	    
	    print '------->\t SAVE ENERGY'
            save_energy(path_of_AMMP)
	    

            fhi_energy = open(energy_file,'r')
#            line0 = fhi_energy.readline()
            line0 = fhi_energy.readline()
            line= line0.split()[1]
	    linevdw = line0.split()[0]
            lineangle = float(line) - float(linevdw)
	    print line



            if (float(line) >= 400.0 and lineangle >= 300.0) or line == 'nan':
               	fhi_ammp = open('output.mol2','r')
               	for lines in fhi_ammp:
                  	fho_badenergy_mols.write(lines)
                    	fho_badenergy_mols.flush()
            else:
                fhi_ammp = open('output.mol2','r')
                for lines in fhi_ammp:
                    fho_built.write(lines)
                    fho_built.flush()


            current_line = ('%s : %s')%(ligand_name,line0)

            fho_energy_file.write(current_line)
            fho_energy_file.flush()
            fho_energy_file.write('\n')
	    
	    
            fhi_warning = open('warning','r')

	    for lines in fhi_warning:
        	if lines != '':
	            fho_totalwarning.write(lines)
        	else:
	            break

            fho_totalwarning.flush()



            
            i = i+1

	command3 = 'rm input.mol2'
    	os.system(command3)
    	command4 = 'rm input_ligand.pdb'
    	os.system(command4)
    	command5 = 'rm input_ligand.ammp'
    	os.system(command5)
    	command6 = 'rm output.mol2 output.ammp'
    	os.system(command6)
    	command10 = 'rm mol warning  energy.txt after_opt.txt preammp_ligand.txt '
    	os.system(command10)
	

        fin = time.time()
        temps = (fin - debut)
        screen = ('------------> Done in : %2.2f sec.')%(temps)
        print "\n"
        print screen
        print "\n"

        fho_built.close()
