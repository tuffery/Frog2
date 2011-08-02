#! /usr/bin/env python




########################################################
###                                                  ###
###                  AMMOS_VLS3D                     ###
###              written by D.Lagorce, M.Miteva      ###
###              http://www.vls3d.com                ###
###                                                  ###
########################################################


from string import *
import os
import os.path
import sys
import BankDivider
import popen2
import time
import shutil
from subprocess import *

########################################################

def usage():
    return """Usage of AMMOS_ProtLig:
    
##########################################
command line = ./AMMOS_ProtLig_sp4.py <ammp.param>
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

def param_file_creator(path_of_AMMP,bank_part,your_choice,fho):
    
    fho.write('path_of_AMMOS_ProtLig= '+path_of_AMMP+'\n')
    fho.write('protein= input_protein.pdb\n')
    fho.write('bank= '+bank_part+'\n')
    fho.write('case_choice= '+str(your_choice)+'\n')
    fho.close()

    
########################################################

def Params_Manager(fhi_params):

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


def preammp_protein_textfile_maker(AMMP_path,protein_file):

    atoms_path = AMMP_path + '/progs/preammp/' ############################# path /src/ of atoms.sp4 !!!! to check
    preammp_fhin = open('preammp_protein.txt','w')
    preammp_fhin.write(atoms_path + 'atoms.sp4\n')
    preammp_fhin.write(atoms_path + 'pdb\n')
    preammp_fhin.write(protein_file+'\n')
    preammp_fhin.write('preammpoutput_protein.ammp\n')
    preammp_fhin.close()


########################################################

def preammp_prot(AMMP_path):

    preammp_command = AMMP_path + '/bin/preammp <preammp_protein.txt'
    os.system(preammp_command)

    
########################################################

def file_organizer(list_of_bank_part,bank_file,choice,warnings,minimized,energy,pdb):

    output_dir = bank_file[:-5] +'_case_'+ str(choice)+'_OUTPUT'
    
    if os.path.exists(output_dir):
        command0 = 'rm -rf ' + output_dir
        os.system(command0)
    
    command1 = 'mkdir ' + output_dir
    os.system(command1)
    command2 = 'mv '+warnings+' '+output_dir
    command3 = 'mv '+minimized+' '+output_dir
    command4 = 'mv '+energy+' '+output_dir
    command5 = 'mv '+pdb+' '+output_dir
    os.system(command2)
    os.system(command3)
    os.system(command4)
    os.system(command5)
    for bank_part in list_of_bank_part:
	command3 = 'rm -rf '+bank_part
	os.system(command3)
    return output_dir
    
########################################################


def Protein_creator(path_of_AMMP):
    
    ammp_protein_command = path_of_AMMP + '/bin/ammp_nongraph <' + path_of_AMMP + '/progs/vls_min/prepare_protein.ammp' #### AMMP ---> LINK PEPTIDE ---> CREATE PROTEIN.AMMP
    os.system(ammp_protein_command)
    


########################################################
 
def File_Manager(file,list_of_dir):

    for dir in list_of_dir:
        shutil.copy(file,dir)

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



########################################################
        
def log_and_results_joiner(bank_listing,bank,your_choice):

    
    
    total_minimized = bank[:-5]+'_case'+str(your_choice)+'_minimized.mol2'
    total_pdb = bank[:-5]+'_case'+str(your_choice)+'_output.pdb'
    total_energy = bank[:-5]+'_case'+str(your_choice)+'_energy.txt'
    total_warning = bank[:-5]+'_case'+str(your_choice)+'_total_warnings.txt'

    fh_energy = open(total_energy,'w')
    title = 'LIGAND : BEFORE_EXTERNAL : AFTER_EXTERNAL : BEFORE_INTERNAL : AFTER_INTERNAL : BEFORE_TOTAL : AFTER_TOTAL\n'
    fh_energy.write(title)
    fh_energy.close()

    line_minimized = ''
    line_pdb = ''
    line_energy = ''
    line_warning = ''
    
    for bank_part in bank_listing:
        line_minimized +=  bank_part+'/'+bank_part[:-4]+'_case'+str(your_choice)+'_minimized.mol2 '
        line_pdb += bank_part+'/'+bank_part[:-4]+'_case'+str(your_choice)+'_output_pdb.pdb '
        line_energy += bank_part+'/'+bank_part[:-4]+'_case'+str(your_choice)+'_energy.txt '
        line_warning += bank_part+'/'+bank_part[:-4]+'_case'+str(your_choice)+'_total_warnings.txt '
        
        
    command_line_minimized = 'cat ' + line_minimized + '> ' + total_minimized
    command_line_pdb = 'cat ' + line_pdb + '> ' + total_pdb
    command_line_energy = 'cat ' + line_energy + '>> ' + total_energy
    command_line_warning = 'cat ' + line_warning + '> ' + total_warning

    os.system(command_line_minimized)
    os.system(command_line_pdb)
    os.system(command_line_energy)
    os.system(command_line_warning)
    return total_minimized,total_pdb,total_energy,total_warning



########################################################
def Multiconf2Singleconf(fhi_multi,fho_single):

    line_inp = {}
    line_out = {}
    printed_mol = {}

    fho_single.write('LIGAND : BEFORE_EXTERNAL : AFTER_EXTERNAL : BEFORE_INTERNAL : AFTER_INTERNAL : BEFORE_TOTAL : AFTER_TOTAL\n')
  
    numb_inp = 1
    line_inp_curr = fhi_multi.readline()[:-1]
    line_inp [numb_inp] = line_inp_curr
        
    while line_inp_curr != '':
        line_inp_curr = fhi_multi.readline()[:-1]
        numb_inp = numb_inp + 1
        line_inp [numb_inp] = line_inp_curr

    mult_line = line_inp [1].split("$")[1]
    mult_mol = mult_line.split(":")[0]
    printed_mol [1] = mult_mol
    fho_single.write (line_inp [1])	## print first molecule
    fho_single.write ('\n')
    
    i = 2
    numb_printed = 1
	
    while i < numb_inp:
        mult_line = line_inp [i].split(":")[0]
        mult_mol = mult_line.split("$")[1]
	
        j = 1
        while j <= numb_printed:
			
            if mult_mol == printed_mol [j]:
                break
				
            else:	
                j = j + 1
				
        if j > numb_printed:
            numb_printed = numb_printed + 1
            printed_mol [numb_printed] = mult_mol
            fho_single.write (line_inp [i])
            fho_single.write ('\n')

        i = i + 1
        

    
########################################################

    
   
if __name__ == '__main__':

    debut = time.time()
    if len(sys.argv)<2 or len(sys.argv)>2:
        print "\n"
        print "raise IOError"
        print usage()

    else:
        file_params = sys.argv[1] ######### fichier parametres
        shutil.copy(file_params,'param.temp')
        fh_params = open(file_params)
        path_of_AMMP,protein,bank,your_choice = Params_Manager(fh_params)       #### return dir de VLS_AMMP fichier protein fichier ligand
        fh_params.close()


        fh_bank = open(bank)
        bank_renamed = bank[:-5]+'_ID_renamed.mol2'
        fh_bank_renamed = open(bank_renamed,'w')
        rename_mol2forAMMP(fh_bank,fh_bank_renamed)
        fh_bank.close()
        fh_bank_renamed.close()
        
        print '-------> COUNTING THE MOLECULES IN THE BANK'
        print '\n'
        mol_amount = BankDivider.mol_amount(bank_renamed)
        print '-------> THE BANK CONTAINS '+str(mol_amount)+' MOLECULES'

        
        fh_bank = open(bank_renamed)
        print '\n'
        print '\n'
        job = BankDivider.job_amount()
        print job
        print '\n'
    
        res_amount = int(BankDivider.divide_amount(mol_amount,job))
        bankpart_list,basename = BankDivider.divide_bank(bank_renamed,res_amount,job)
        print '\n'
        print '------------> BANK PREPARING FINISHED'
        print "\n"

        list_of_bank_part = BankDivider.bank_moving(bankpart_list)
       
        fh_bank.close()
        
        fh_bank = open(bank_renamed)

        new_input_protein = 'input_protein.pdb'
        shutil.copy(protein,new_input_protein)
        print '------------> COPYING input_protein.pdb'
        
        print '------------> MAKING PROT TEXT_FILE MAKER'
        preammp_protein_textfile_maker(path_of_AMMP,new_input_protein) ### PREPARE PREAMMP PROTEIN INPUT TEXT FILE
        print '------------> PROT PREAMMP'
        preammp_prot(path_of_AMMP)                           ####   PREAMMP --->  PREPARE PROTEIN INPUT FOR AMMP
	print '------------> PROTEIN CREATOR'
        Protein_creator(path_of_AMMP)                        ### CREATE PROTEIN.AMMP

        print '------------> OXT error CORRECTED'
	protein_file = 'protein.ammp'
        File_Manager(protein_file,list_of_bank_part)
        File_Manager(new_input_protein,list_of_bank_part)
        #os.remove(protein_file)
        print '------------> BEGIN LOOP'
        list_of_process = []







        
        command_process = path_of_AMMP+'progs/Python_scripts/MinProteinLigands_sp4.py param.temp'

        working_path = os.getcwd()
        
        for bank_part in list_of_bank_part:

            os.chdir(bank_part)
            fout = open( bank_part[:-4] + ".out" , 'w') ### fh de sortie ###
            ferr =  open (  bank_part[:-4]+ ".err" , 'w') ### fh error ###
            new_fho_params = open('ammp.param','w')
            bank_part_mol2 = bank_part[:-4]+'.mol2'
            param_file_creator(path_of_AMMP,bank_part_mol2,your_choice,new_fho_params)
            shutil.copy('ammp.param','param.temp')
            
            try: 
                list_of_process.append( Popen(command_process,
                                              shell = True,
                                              stdout = fout,
                                              stdin = None,
                                              stderr = ferr,
                                              close_fds = True
                                              ))
                
                print "\n"
                print 'PROCESS --',command_process,'-- LAUNCHED'
            except OSError, err:

                message = "Execution failed for bank :" + bank_part[:-4]+'.mol2'+"\n"
                raise OSError , message
            
            os.chdir(working_path)
            
        for process in list_of_process:

            process.wait()

        
        total_minimized,total_pdb,total_energy,total_warning = log_and_results_joiner(list_of_bank_part,bank,your_choice)

        
        fin = time.time()
        temps = (fin - debut)
        screen = ('------------> Done in : %2.2f sec.')%(temps)
        print "\n"
        print screen
        print "\n"
        
                              
			
        output_dir = file_organizer(list_of_bank_part,bank,your_choice,total_warning,total_minimized,total_energy,total_pdb)
        
################################## Multiconf-> signle conf####################

        os.chdir(output_dir)
        if your_choice == 5:
            os.remove(total_pdb)
        fhi_multi = open(total_energy)
        file_temp = total_energy[:-4]+'.temp'
        fho_multi_temp = open(file_temp,'w')
        line = fhi_multi.readline()
        for lines in fhi_multi:
            fho_multi_temp.write(lines)
        fho_multi_temp.close()
        os.system('sort -n -t":" -k3 '+file_temp+' -o '+file_temp+'2')
        fhi_total = open(file_temp+'2')
        fho_single = open(total_energy[:-4]+'_ranked.txt','w')
        Multiconf2Singleconf(fhi_total,fho_single)
        commandmm = 'rm ' + file_temp + ' '+file_temp+'2'
        os.system(commandmm)
        os.chdir(working_path)
        os.remove('param.temp')
        os.remove('preammp_protein.txt')

## - to verify all open and close files
