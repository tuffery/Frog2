#! /usr/bin/env python



########################################################
###                                                  ###
###   Minimization of Protein-ligand Interactions    ###
###              written by D.Lagorce,M.Miteva       ###
###              http://www.vls3d.com                ###
###                                                  ###
########################################################

import os
import os.path
import sys
from string import *
import popen2
import time
import shutil

########################################################

def usage():
    return """Usage of MinProteinLigands:
    
##########################################
command line = ./MinProteinLigands_sp4.py <ammp.param>
##########################################

<ammp.param> is the paramaters file which needs to be edited with your
             configuration. This file is located in your project directory.

             
Possible choice of active/inactive atoms:

   - Case <1>: Make active all atoms of the protein and ligand.
   
   - Case <2>: Make inactive atoms in the protein backbone.
   
   - Case <3>: Make active all atoms of the protein failed in
               the sphere around the ligand.
               
   - Case <4>: Make active all atoms of the protein failed in
               the sphere around the ligand, except ones in the
               protein backbone in the sphere.

    - Case <5>: Make active all atoms of the ligand while the whole
                protein is rigid.

   """

########################################################

def params(fhi_params):

    line = fhi_params.readline()
    ##########
    path_AMMP = line.split()[1]

    line = fhi_params.readline()
    ########
    protein = line.split()[1]

    line = fhi_params.readline()
    ########
    bank = line.split()[1]

    line = fhi_params.readline()
    ########
    choice = atoi(line.split()[1])


    return path_AMMP,protein,bank,choice

    
#######################################################

def preammp_ligand_textfile_maker(AMMP_path):

    command ='pwd'
    fhin,fhout = popen2.popen2(command)
    for line in fhin:
        working_path = line

    atoms_path = AMMP_path + '/progs/preammp/' ############################# 
    preammp_fhin = open('preammp_ligand.txt','w')
    preammp_fhin.write(atoms_path + 'atoms.sp4\n')
    preammp_fhin.write(working_path)
    preammp_fhin.write('input_ligand.pdb\n')
    preammp_fhin.write('input_ligand.ammp\n')
    preammp_fhin.close()


#######################################################
## AMMP_X_threads

def preammp_protein_textfile_maker(AMMP_path,protein_file):

    command ='pwd'
    fhin,fhout = popen2.popen2(command)
    for line in fhin:
        working_path = line
 
    atoms_path = AMMP_path + '/progs/preammp/' ############################# 
    preammp_fhin = open('preammp_protein.txt','w')
    preammp_fhin.write(atoms_path + 'atoms.sp4\n')
    preammp_fhin.write(atoms_path + 'pdb\n')
    preammp_fhin.write(protein_file+'\n')
    preammp_fhin.write('preammpoutput_protein.ammp\n')
    preammp_fhin.close()


########################################################

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
    command7 = 'mv *_minimized* ' + output_dir 
    os.system(command7)
    command8 = 'rm energy.txt'
    os.system(command8)
    command2 = 'mv ' + warnings + ' ' + output_dir
    os.system(command2)
    command9 = 'mv ' + energy_file + ' ' + output_dir
    os.system(command9)
    command10 = 'rm mol'
    os.system(command10)
    command11 = 'rm output.ammp'
    os.system(command11)
    command12 = 'rm before_opt.txt'
    os.system(command12)
    command13 = 'rm after_opt.txt'
    os.system(command13)
   
    if choice >= 3:
        command1 = 'rm active.ammp'
        os.system(command1)
    if choice == 2:
        command2 = 'rm inactive.ammp'
        os.system(command2)


########################################################

def Get_Mol_Amount(fhi):

    amount = 0
    line = fhi.readline()
    while line != '':

        if line[:17] == '@<TRIPOS>MOLECULE':
            amount = amount +1
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
    index[mol]=(start,stop)
    print '------>\t'+mol
    return index


########################################################

def GetList_ligand(dico_pos):

    badlist_ligands = dico.keys()
    list_values = dico.values()
    list_values.sort()
    good_list_ligands = []
    for values in list_values :
        for ligands in badlist_ligands:
            if dico_pos[ligands] == values:
                print ligands
                good_list_ligands.append(ligands)

    return  good_list_ligands


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

def preammp_prot(AMMP_path):

    preammp_command = AMMP_path + '/bin/preammp <preammp_protein.txt'
    os.system(preammp_command)
    
########################################################
    
def ammp(choice,AMMP_path):

    if choice == 1:
        ammp_command = AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/min_case1.ammp'
        
    elif choice == 2:
        active_command = AMMP_path + '/bin/active_case2'
        os.system(active_command)
        ammp_command =  AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/min_case2.ammp'
        
    elif choice == 3:
        active_command = AMMP_path + '/bin/active_case3'
        os.system(active_command)
        ammp_command =  AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/min_case3.ammp'    

    elif choice == 4:
        active_command = AMMP_path + '/bin/active_case4'
        os.system(active_command)
        ammp_command =   AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/min_case4.ammp'  

    elif choice == 5:
        active_command = AMMP_path + '/bin/active_case5'
        os.system(active_command)
        ammp_command =  AMMP_path + '/bin/ammp_nongraph <' + AMMP_path + '/progs/vls_min/min_case5.ammp'    
    else:
        print "raise IOError"
        print usage()

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

def renumb(AMMP_path):

    renumb_cmd = AMMP_path + '/bin/renumb_ligand'
    os.system(renumb_cmd)


########################################################

def save_pdb(AMMP_path,choice):

    if choice == 1:
        save_cmd = AMMP_path + '/bin/ammp_to_pdb_case1'
        flag = 1
    if choice == 2:
        save_cmd = AMMP_path + '/bin/ammp_to_pdb_case2'
        flag = 1
    if choice == 3:
        save_cmd = AMMP_path + '/bin/ammp_to_pdb_case3'
        flag = 1
    if choice == 4:
        save_cmd = AMMP_path + '/bin/ammp_to_pdb_case4'
	flag = 1

    if flag == 1:
	os.system(save_cmd)


########################################################

def warning_file_parser(fhi_warning,fho_totalwarning):

    for lines in fhi_warning:
        if lines != '':
            fho_totalwarning.write(lines)
        else:
            break
    

########################################################

def rename_mol2forAMMP(fhi,fho):

    line = fhi.readline()
    i = 0
    while line != '':
        line = line[:-1]
        if line == '@<TRIPOS>MOLECULE':
            i=i+1
            fho.write(line)
            fho.write('\n')
            line = fhi.readline()
            line = line[:-1]
            namenew = str(i)+'$'+line
            fho.write(namenew)
            fho.write('\n')
	     
        else:
            if line[0:3] != '###':
                fho.write(line)
                fho.write('\n')
                
        line = fhi.readline()

    
if __name__ == '__main__':


    
    if len(sys.argv)<2 or len(sys.argv)>2:
        print "\n"
        print "raise IOError"
        print usage()

    else:
        
        fh_params = open('param.temp')
        path_of_AMMP_temp,protein_temp,bank_part,your_choice_temp = params(fh_params)       
        bank = bank_part[:-4]+'mol2'

        fh_bank = open(bank)
        mol_amount = Get_Mol_Amount(fh_bank)
              
        fh_bank.seek(0,0)
        preammp_ligand_textfile_maker(path_of_AMMP_temp)

      
        minimized_bank_file = bank[:-5] + '_case' + str(your_choice_temp) + '_minimized.mol2'
        pdb_output = bank[:-5] + '_case' + str(your_choice_temp) + '_output_pdb.pdb'
        energy_bank_file = bank[:-5] + '_case' + str(your_choice_temp) + '_energy.txt'
        
        fho_minimized = open(minimized_bank_file,'w')
        fho_energy_file = open(energy_bank_file,'w')
        totalwarning = bank[:-5] + '_case' + str(your_choice_temp) + '_total_warnings.txt'
        fho_totalwarning = open(totalwarning,'w')

        fho_output_pdb = open(pdb_output,'w')


        i = 0
        while i <= mol_amount - 1:
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
            ligand_convertion2pdb(path_of_AMMP_temp)
            print '------->\t PREAMMP'
            preammp(path_of_AMMP_temp)
            print '------->\t RENUMBERING'
            renumb(path_of_AMMP_temp)
            print '------->\t AMMP'
            ammp(your_choice_temp,path_of_AMMP_temp)
            print '------->\t AMMP2MOL2'
            ammp2mol2(path_of_AMMP_temp)
            print '------->\t SAVE ENERGY'
            save_energy(path_of_AMMP_temp)
            fhi_ammp = open('output.mol2','r')
            for lines in fhi_ammp:
                fho_minimized.write(lines)
                fho_minimized.flush()
            fhi_energy = open(energy_file,'r')
            line = fhi_energy.readline()
            current_line = ('%s : %s')%(ligand_name,line)
            fho_energy_file.write(current_line)
            fho_energy_file.flush()
            fho_energy_file.write('\n')
            fho_energy_file.flush()
            fhi_warning = open('warning','r')

	    for lines in fhi_warning:
        	if lines != '':
	            fho_totalwarning.write(lines)
        	else:
	            break

            fho_totalwarning.flush()

            if your_choice_temp <= 4:
                print '------->\t AMMP2PDB'
                save_pdb(path_of_AMMP_temp,your_choice_temp)
                fhi_pdb = open('output_protein.pdb','r')
                fho_output_pdb.write('REMARK ' + ligand_name)
                fho_output_pdb.write('\n')
                for lines in fhi_pdb:
                    fho_output_pdb.write(lines)
                    fho_output_pdb.flush()
                fho_output_pdb.write('\n')
                fhi_pdb.close()
     
            
            i = i+1

        fho_minimized.close()
        fho_energy_file.close()
        fho_output_pdb.close()
	fhi_warning.close()
	fho_totalwarning.close()



