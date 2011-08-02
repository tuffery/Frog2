#! /usr/bin/env python


from string import *
import sys
import os
import os.path
import popen2
import operator
import time

def mol_amount(fi):
#### count the molecules in the bank input
    
    amount = 0
    command = 'grep -wc "MOLECULE" '
    command = command + fi
    fhin,fhout = popen2.popen2(command)
    for line in fhin:
        amount = int(line)
    return amount 


def job_amount():
#### how many jobs to run
    
    jobs = raw_input('How many jobs do you want to launch ?  ')
    jobs = int(jobs)
    return jobs
    
    
def divide_amount(res,jobs):
#### calculate the amount of molecules for each part of bank_part
    res_fin =  operator.div(res,jobs)
    line = ('....:::: %s jobs requested ::::....')%(jobs)
    print line
    return res_fin


def divide_bank(file,amount,job):
#### dividing the bank in different bank_part
    i=0
    j=1
    basename = (os.path.basename(file))[:-5]
    newbank = basename + '_part' + str(j) + '.mol2'
    fho = open(newbank,'w')
    fhi = open(file,'r')
    line = fhi.readline()
    bank_part_list = []
    while line != '':
        line = line[:-1]
        if line == '@<TRIPOS>MOLECULE':
            i = i+1
            if j < job: ## test 
                if i <= amount:
                    fho.write(line)
                    fho.write('\n')

                else:
                    if j < job:
                        fhi.seek(head)
                        screen = ('---------> %s written')%(newbank)
                        bank_part_list.append(newbank)
                        print screen
                        fho.close()
                        j = j+1
                        newbank = basename + '_part' + str(j) + '.mol2'
                        screen = ('---------> %s written')%(newbank)
                        fho = open(newbank,'w')
                        i = 0


            else:
                fho.write(line)
                fho.write('\n')                 
        else:
            fho.write(line)
            fho.write('\n')
        
        head = fhi.tell()
        line = fhi.readline()
        
    screen = ('---------> %s written')%(newbank)
    print screen
    bank_part_list.append(newbank)
    return bank_part_list,basename

    fho.close()
    fhi.close()


def bank_moving(bank_part_list):
#### creating different directories from bank_parts
#### and moving each bank_part in the good directory

    list_dir = []
    for bank_part in bank_part_list:
        dirname = bank_part[:-4] + 'dir'
        cmd = 'mkdir '+ dirname
        os.system(cmd)
        cmd2 = 'mv ' + bank_part + ' ' + dirname
        os.system(cmd2)
        list_dir.append(dirname)
    return list_dir
        
         

    

if __name__=='__main__':

    debut = time.time()
    usage = "Usage: You have to type:  ./BankDivider.py [Bank_fileName.mol2]"
    
    if len(sys.argv) != 2:
        print usage
        sys.exit(1)

    if sys.argv[1] == "":
        print "\n"
        print "---- no input file ----"
        print usage
        sys.exit(1)
        
    else:
        bank = sys.argv[1]
        print '\n'
        print 'Please wait, reading and opening the whole bank file'
        print '\n'
    
        fh_bank = open(bank,'r')
    
        res = mol_amount(bank)
    
        print '\n'
        print ('---------> The bank file contains %i molecules')%(res)
        print '\n'
    
        job = job_amount()
        if res <= job:
            job = res

                
        print '\n'
    
        if job != 1 and job != 2 and job != 4:
            print 'choice not supported, restart python command with 1,2 or 4 jobs !'
            print '\n'
        
        if job == 1:
            print '---------> Just one job requested, therefore the bank file will not be modified'
            print '\n'
          
	if job >= 2: 
            res_amount = int(divide_amount(res,job))
            rest = res-job*res_amount 
            print '\n'
            print ('---------> Creating new bank files each containing %i molecules.')%(res_amount)
            job_list,basename = divide_bank(bank,res_amount,job)
            print '\n'
            bank_moving(job_list)
            fin = time.time()
            temps = (fin - debut)/60
            screen = ('------------> Done in : %2.2f min.')%(temps)
            print "\n"
            print screen
            print "\n"
    
