#!/usr/bin/env python3
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import scripts 
import numpy as np

def restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot):

    rst = []
    atm_num = []
    hvy_h = []
    hvy_g = []
    mlines = []
    msk = []
    pdb_file = ('vac.pdb')
    guest_pdb_file = ('vac_guest.pdb')
    host_pdb_file = ('vac_host.pdb')

    # Restraint identifiers
    recep_c = '#Rec_C'
    recep_d = '#Rec_D'
    lign_tr = '#Lig_TR'
    lign_c = '#Lig_C'
    lign_d = '#Lig_D'

    # Change to simulation directory
    if stage != 'fe':
      os.chdir(guest)

    # Find anchors
    with open('anchors-%s.txt' % guest.lower(), 'r') as f:
      data = f.readline().split()    
      
      l1_atom = data[0].strip()   
      l2_atom = data[1].strip()   
      l3_atom = data[2].strip()   
      lig_res = str((int(final_host_num) + 2))
      if comp == 'e' or comp == 'v':
        lig_res = str((int(final_host_num) + 3))
      L1 = ':'+lig_res+'@'+l1_atom
      L2 = ':'+lig_res+'@'+l2_atom
      L3 = ':'+lig_res+'@'+l3_atom

    if (comp != 'c' and comp != 'r'):

      # Adjust host anchors to the new residue numbering

      h1_resid = H1.split('@')[0][1:]
      h2_resid = H2.split('@')[0][1:]
      h3_resid = H3.split('@')[0][1:]

      h1_atom = H1.split('@')[1]
      h2_atom = H2.split('@')[1]
      h3_atom = H3.split('@')[1]

      p1_resid = str(int(h1_resid) + 1)
      p2_resid = str(int(h2_resid) + 1)
      p3_resid = str(int(h3_resid) + 1)

      H1 = ":"+p1_resid+"@"+h1_atom
      H2 = ":"+p2_resid+"@"+h2_atom
      H3 = ":"+p3_resid+"@"+h3_atom


      # Get host heavy atoms
      with open('./full.pdb') as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        for i in range(0, len(lines)):
          if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
            if int(lines[i][22:26].strip()) >= 2 and int(lines[i][22:26].strip()) <= final_host_num+1:
              data = lines[i][12:16].strip()
              if data[0] != 'H':
                hvy_h.append(lines[i][6:11].strip()) 


      if (comp == 'e' or comp == 'v'):

        hvy_h = []
        hvy_g = []

        # Adjust anchors

        p1_resid = str(int(h1_resid) + 2)
        p2_resid = str(int(h2_resid) + 2)
        p3_resid = str(int(h3_resid) + 2)

        H1 = ":"+p1_resid+"@"+h1_atom
        H2 = ":"+p2_resid+"@"+h2_atom
        H3 = ":"+p3_resid+"@"+h3_atom

        # Get host heavy atoms
        with open('./full.pdb') as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line) # Non-blank lines in a list   
          for i in range(0, len(lines)):
            if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
              if int(lines[i][22:26].strip()) >= 3 and int(lines[i][22:26].strip()) <= final_host_num+2:
                data = lines[i][12:16].strip()
                if data[0] != 'H':
                  hvy_h.append(lines[i][6:11].strip()) 

        # Get bulk guest heavy atoms
        with open('./full.pdb') as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line) # Non-blank lines in a list   
          if comp == 'e':
            for i in range(0, len(lines)):
              if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                if lines[i][22:26].strip() == str(int(lig_res) + 2):
                  data = lines[i][12:16].strip()
                  if data[0] != 'H':
                    hvy_g.append(lines[i][6:11].strip()) 
          if comp == 'v':
            for i in range(0, len(lines)):
              if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                if lines[i][22:26].strip() == str(int(lig_res) + 1):
                  data = lines[i][12:16].strip()
                  if data[0] != 'H':
                    hvy_g.append(lines[i][6:11].strip()) 

    # Adjust anchors for ligand only
    if (comp == 'c'):
      L1 = L1.replace(':'+lig_res, ':1') 
      L2 = L2.replace(':'+lig_res, ':1') 
      L3 = L3.replace(':'+lig_res, ':1') 

    
    # Get a relation between atom number and masks
    atm_num = scripts.num_to_mask(pdb_file)
    guest_atm_num = scripts.num_to_mask(guest_pdb_file)

    # Get number of ligand atoms
    with open('./vac_guest.pdb') as myfile:
        data = myfile.readlines()
        vac_atoms = int(data[-3][6:11].strip())

    # Define conformatonal restraints on the host

    if host_rest_type == 'distances' and comp != 'c':
      rst.append(''+H1+' '+H2+'') 
      rst.append(''+H2+' '+H3+'') 
      rst.append(''+H3+' '+H1+'') 
      nd = 3
    elif host_rest_type == 'dihedrals' and comp != 'c':
      host_atm_num = scripts.num_to_mask(host_pdb_file)
      # Get host dihedral restraints from ligand parameter/pdb file
      spool = 0
      with open('./vac_host.prmtop') as fin:
          lines = (line.rstrip() for line in fin)
          lines = list(line for line in lines if line) # Non-blank lines in a list   
          for line in lines:
            if 'FLAG DIHEDRALS_WITHOUT_HYDROGEN' in line:
              spool=1
            elif 'FLAG EXCLUDED_ATOMS_LIST' in line:
              spool=0
            if spool != 0 and (len(line.split()) > 3):
              mlines.append(line)


      for i in range(0, len(mlines)):
        data = mlines[i].split()
        if int(data[3]) > 0:
          anum = []
          for j in range (0, len(data)): 
            anum.append(abs(int(data[j])//3)+1)      
          msk.append('%s %s %s %s' %(host_atm_num[anum[0]], host_atm_num[anum[1]], host_atm_num[anum[2]], host_atm_num[anum[3]]))   
             
      for i in range(0, len(mlines)):
        data = mlines[i].split()
        if len(data) > 7:
          if int(data[8]) > 0:
            anum = []
            for j in range (0, len(data)): 
              anum.append(abs(int(data[j])//3)+1)      
            msk.append('%s %s %s %s' %(host_atm_num[anum[5]], host_atm_num[anum[6]], host_atm_num[anum[7]], host_atm_num[anum[8]]))   
             
      excl = msk[:]
      ind = 0
      mat = []
      for i in range(0, len(excl)):
         data = excl[i].split()
         for j in range(0, len(excl)):   
           if j == i:
             break 
           data2 = excl[j].split()
           if (data[1] == data2[1] and data[2] == data2[2]) or (data[1] == data2[2] and data[2] == data2[1]):
             ind = 0
             for k in range(0, len(mat)):
               if mat[k] == j:
                 ind = 1
             if ind == 0:
               mat.append(j) 

      for i in range(0, len(mat)):
        msk[mat[i]]= ''

      if (comp == 'a' or comp == 'l' or comp == 't' or comp == 'q'):
        msk = list(filter(None, msk))
#        for m in msk:
#          print (m)
        mskt = list(msk)
        for i in range(0, len(mskt)):
          data = mskt[i].split()
          for j in range(0, len(data)):
            for k in range(0, final_host_num):
              num = int(data[j].split()[0][1])
              if num == k+1:
                p = k+1
                n = k+2
                data[j] = data[j].replace(':'+str(p)+'',':'+str(n)+'')
                break
          mskt[i] = '%s %s %s %s' %(data[0], data[1], data[2], data[3])
        msk = list(mskt)
      elif (comp == 'e' or comp == 'v'):
        msk = list(filter(None, msk))
        mskt = list(msk)
        for i in range(0, len(mskt)):
          data = mskt[i].split()
          for j in range(0, len(data)):
            for k in range(0, final_host_num):
              num = int(data[j].split()[0][1])
              if num == k+1:
                p = k+1
                n = k+3
                data[j] = data[j].replace(':'+str(p)+'',':'+str(n)+'')
                break
            mskt[i] = '%s %s %s %s' %(data[0], data[1], data[2], data[3])
        msk = list(mskt)
      else: # c or r
        msk = list(filter(None, msk))

      for i in range(0, len(msk)):
        rst.append(msk[i])
      
      rsta = rst[:]

      nd = len(rst) 

    # Define translational/rotational restraints on the guest


    rst.append(''+H1+' '+L1+'')
    rst.append(''+H2+' '+H1+' '+L1+'')
    rst.append(''+H3+' '+H2+' '+H1+' '+L1+'')
    rst.append(''+H1+' '+L1+' '+L2+'')
    rst.append(''+H2+' '+H1+' '+L1+' '+L2+'')
    rst.append(''+H1+' '+L1+' '+L2+' '+L3+'')

    # New restraints for ligand only
    if (comp == 'c'):
      rst = []
       
    # Get guest dihedral restraints from ligand parameter/pdb file

    mlines = []
    msk = []
    spool = 0
    with open('./vac_guest.prmtop') as fin:
        lines = (line.rstrip() for line in fin)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        for line in lines:
          if 'FLAG DIHEDRALS_WITHOUT_HYDROGEN' in line:
            spool=1
          elif 'FLAG EXCLUDED_ATOMS_LIST' in line:
            spool=0
          if spool != 0 and (len(line.split()) > 3):
            mlines.append(line)


    for i in range(0, len(mlines)):
      data = mlines[i].split()
      if int(data[3]) > 0:
        anum = []
        for j in range (0, len(data)): 
          anum.append(abs(int(data[j])//3)+1)      
        msk.append('%s %s %s %s' %(guest_atm_num[anum[0]], guest_atm_num[anum[1]], guest_atm_num[anum[2]], guest_atm_num[anum[3]]))   
           
    for i in range(0, len(mlines)):
      data = mlines[i].split()
      if len(data) > 7:
        if int(data[8]) > 0:
          anum = []
          for j in range (0, len(data)): 
            anum.append(abs(int(data[j])//3)+1)      
          msk.append('%s %s %s %s' %(guest_atm_num[anum[5]], guest_atm_num[anum[6]], guest_atm_num[anum[7]], guest_atm_num[anum[8]]))   
           
    excl = msk[:]
    ind = 0
    mat = []
    for i in range(0, len(excl)):
       data = excl[i].split()
       for j in range(0, len(excl)):   
         if j == i:
           break 
         data2 = excl[j].split()
         if (data[1] == data2[1] and data[2] == data2[2]) or (data[1] == data2[2] and data[2] == data2[1]):
           ind = 0
           for k in range(0, len(mat)):
             if mat[k] == j:
               ind = 1
           if ind == 0:
             mat.append(j) 

    for i in range(0, len(mat)):
      msk[mat[i]]= ''

    if (comp != 'c'):
      msk = filter(None, msk) 
      msk = [m.replace(':1',':'+lig_res) for m in msk]


    for i in range(0, len(msk)):
      rst.append(msk[i])

    # New restraints for protein only
    if (comp == 'r'):
      if host_rest_type == 'distances':
        rst = []
        rst.append(''+H1+' '+H2+'') 
        rst.append(''+H2+' '+H3+'') 
        rst.append(''+H3+' '+H1+'') 
      elif host_rest_type == 'dihedrals':
        rst = []
        for i in range(0, len(rsta)):
          rst.append(rsta[i])   

    # Get initial restraint values for references

    assign_file = open('assign.in', 'w')
    assign_file.write('%s  %s  %s  %s  %s  %s  %s\n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3))
    assign_file.write('parm full.hmr.prmtop\n')
    assign_file.write('trajin full.inpcrd\n')
    for i in range(0, len(rst)):
        arr = rst[i].split()
        if len(arr) == 2:
          assign_file.write('%s %s %s'%('distance r'+str(i), rst[i], 'noimage out assign.dat\n'))
        if len(arr) == 3:
          assign_file.write('%s %s %s'%('angle r'+str(i), rst[i], 'out assign.dat\n'))
        if len(arr) == 4:
          assign_file.write('%s %s %s'%('dihedral r'+str(i), rst[i], 'out assign.dat\n'))

    assign_file.close()
    sp.call('cpptraj -i assign.in > assign.log', shell=True)


    # Assign reference values for restraints
    with open('./assign.dat') as fin:
        lines = (line.rstrip() for line in fin)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        vals = lines[1].split()
        vals.append(vals.pop(0))
        del vals[-1]

    if stage == 'equil':
      rdhf = weight*rest[0]/100
      rdsf = weight*rest[1]/100
      ldf =  weight*rest[2]/100
      laf =  weight*rest[3]/100
      ldhf = weight*rest[4]/100
      rcom = rest[5]
    elif comp == 'l' or comp == 'c':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = 0
      laf = 0
      ldhf = weight*rest[4]/100
      rcom = rest[5]
    elif comp == 'a' or comp == 'r':
      rdhf = weight*rest[0]/100
      rdsf = weight*rest[1]/100
      ldf = 0
      laf = 0
      ldhf = 0
      rcom = rest[5]
    elif comp == 't':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = weight*rest[2]/100
      laf = weight*rest[3]/100
      ldhf = rest[4]
      rcom = rest[5]
    elif comp == 'v' or comp == 'e':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = rest[2]
      laf = rest[3]
      ldhf = rest[4]
      rcom = rest[5]
      lcom = rest[6]

    if guest_rot == 'yes':
      laf1 = 0
    else:
      laf1 = laf

    # Write AMBER restraint file for the full system 
    if (comp != 'c' and comp != 'r'):
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  %s  %s  %s \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
        data = rst[i].split()
        # Host conformational restraints
        if i < nd:
          if len(data) == 2:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), rdsf, rdsf, recep_c))
          elif len(data) == 4:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, rdhf, rdhf, recep_d))
        # Guest translational/rotational restraints
        if i >= nd and i < 5+nd and comp != 'a':
          if len(data) == 2:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldf, ldf, lign_tr))
          elif len(data) == 3:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','  
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
          elif len(data) == 4:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))
          if comp == 'e':
            if i == (nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldf, ldf, lign_tr))
            if i == (1+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
            if i == (2+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))
            if i == (3+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
            if i == (4+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))
        if i == (5+nd) and comp != 'a':
          if len(data) == 4:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf1, laf1, lign_tr))
          if comp == 'e':
            if i == (5+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf1, laf1, lign_tr))
        # Guest conformational restraints
        elif i >= 6+nd and comp != 'a':
          nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
          if comp == 'v':
            nums2 = str(atm_num.index(data[0])+vac_atoms)+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','
            disang_file.write('%s %-23s '%('&rst iat=', nums2))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
          if comp == 'e':
            nums2 = str(atm_num.index(data[0])+vac_atoms)+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','
            disang_file.write('%s %-23s '%('&rst iat=', nums2))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
            nums3 = str(atm_num.index(data[0])+2*vac_atoms)+','+str(atm_num.index(data[1])+2*vac_atoms)+','+str(atm_num.index(data[2])+2*vac_atoms)+','+str(atm_num.index(data[3])+2*vac_atoms)+','
            disang_file.write('%s %-23s '%('&rst iat=', nums3))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
            nums4 = str(atm_num.index(data[0])+3*vac_atoms)+','+str(atm_num.index(data[1])+3*vac_atoms)+','+str(atm_num.index(data[2])+3*vac_atoms)+','+str(atm_num.index(data[3])+3*vac_atoms)+','
            disang_file.write('%s %-23s '%('&rst iat=', nums4))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))


      # COM restraints
      if comp != 'r' and comp != 'c':
        cv_file = open('cv.in', 'w')
        cv_file.write('cv_file \n')
        cv_file.write('&colvar \n')
        cv_file.write(' cv_type = \'COM_DISTANCE\' \n')
        cv_file.write(' cv_ni = %s, cv_i = 1,0,' % str(len(hvy_h)+2))
        for i in range(0, len(hvy_h)):
          cv_file.write(hvy_h[i])
          cv_file.write(',')
        cv_file.write('\n')
        cv_file.write(' anchor_position = %10.4f, %10.4f, %10.4f, %10.4f \n' % (float(0.0), float(0.0), float(0.0), float(999.0)))
        cv_file.write(' anchor_strength = %10.4f, %10.4f, \n' % (rcom, rcom))
        cv_file.write('/ \n')
        if comp == 'e' or comp == 'v':
          cv_file.write('&colvar \n')
          cv_file.write(' cv_type = \'COM_DISTANCE\' \n')
          cv_file.write(' cv_ni = %s, cv_i = 2,0,' % str(len(hvy_g)+2))
          for i in range(0, len(hvy_g)):
            cv_file.write(hvy_g[i])
            cv_file.write(',')
          cv_file.write('\n')
          cv_file.write(' anchor_position = %10.4f, %10.4f, %10.4f, %10.4f \n' % (float(0.0), float(0.0), float(0.0), float(999.0)))
          cv_file.write(' anchor_strength = %10.4f, %10.4f, \n' % (lcom, lcom))
          cv_file.write('/ \n')
        cv_file.close()



      # Analysis of simulations

      if (comp != 'l' and comp != 'a'):
        restraints_file = open('restraints.in', 'w')
        restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage))
        restraints_file.write('noexitonerror\n')
        restraints_file.write('parm vac.prmtop\n')
        for i in range(2,11):
          restraints_file.write('trajin md%02.0f.nc\n' % i)
        for i in range(nd, 6+nd):
          arr = rst[i].split()
          if len(arr) == 2:
            restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
          if len(arr) == 3:
            restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
          if len(arr) == 4:
            restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
      elif (comp == 'a'):
        restraints_file = open('restraints.in', 'w')
        restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage))
        restraints_file.write('noexitonerror\n')
        restraints_file.write('parm vac.prmtop\n')
        for i in range(2,11):
          restraints_file.write('trajin md%02.0f.nc\n' % i)
        for i in range(0, nd):
          arr = rst[i].split()
          if len(arr) == 2:
            restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
          if len(arr) == 3:
            restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
          if len(arr) == 4:
            restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
      elif (comp == 'l'):
        restraints_file = open('restraints.in', 'w')
        restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage))
        restraints_file.write('noexitonerror\n')
        restraints_file.write('parm vac.prmtop\n')
        for i in range(2,11):
          restraints_file.write('trajin md%02.0f.nc\n' % i)
        for i in range(6+nd, len(rst)):
          arr = rst[i].split()
          if len(arr) == 2:
            restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
          if len(arr) == 3:
            restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
          if len(arr) == 4:
            restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))

    elif comp == 'c':
      while '' in rst:
        rst.remove('')
    # Write restraint file for ligand system 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  %s  %s  %s \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
        data = rst[i].split()
        # Ligand conformational restraints
        nums = str(guest_atm_num.index(data[0]))+','+str(guest_atm_num.index(data[1]))+','+str(guest_atm_num.index(data[2]))+','+str(guest_atm_num.index(data[3]))+','
        disang_file.write('%s %-23s '%('&rst iat=', nums))
        disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
      # Analysis of simulations
      restraints_file = open('restraints.in', 'w')
      restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage))
      restraints_file.write('noexitonerror\n')
      restraints_file.write('parm vac.prmtop\n')
      for i in range(2,11):
        restraints_file.write('trajin md%02.0f.nc\n' % i)
      for i in range(0, len(rst)):
        arr = rst[i].split()
        if len(arr) == 2:
          restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
        if len(arr) == 3:
          restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
        if len(arr) == 4:
          restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
    elif comp == 'r':
      # Write restraint file for host system 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
        data = rst[i].split()
        # Host conformational restraints
        if len(data) == 2:
          nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), rdsf, rdsf, recep_c))
        if len(data) == 4:
          nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180 , rdhf, rdhf, recep_d))
      # Analysis of simulations
      restraints_file = open('restraints.in', 'w')
      restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', H1, H2, H3, L1, L2, L3, 'stage = '+stage))
      restraints_file.write('noexitonerror\n')
      restraints_file.write('parm vac.prmtop\n')
      for i in range(2,11):
        restraints_file.write('trajin md%02.0f.nc\n' % i)
      for i in range(0, len(rst)):
        arr = rst[i].split()
        if len(arr) == 2:
          restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
        if len(arr) == 3:
          restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
        if len(arr) == 4:
          restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))

    disang_file.write('\n')
    disang_file.close()

    if stage != 'fe':
      os.chdir('../')

def sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, rng):

    if stage != 'fe':
      # Copy folder for running simulations
      if not os.path.exists('run_files'):
        try:
          shutil.copytree('../run_files', './run_files')
        # Directories are the same
        except shutil.Error as e:
          print('Directory not copied. Error: %s' % e)
        # Any error saying that the directory doesn't exist
        except OSError as e:
          print('Directory not copied. Error: %s' % e)
      if hmr == 'no': 
        replacement = 'full.prmtop'
        for dname, dirs, files in os.walk('./run_files'):
          for fname in files:
            fpath = os.path.join(dname, fname)
            with open(fpath) as f:
              s = f.read()
              s = s.replace('full.hmr.prmtop', replacement)
            with open(fpath, "w") as f:
              f.write(s)
      elif hmr == 'yes': 
        replacement = 'full.hmr.prmtop'
        for dname, dirs, files in os.walk('./run_files'):
          for fname in files:
            fpath = os.path.join(dname, fname)
            with open(fpath) as f:
              s = f.read()
              s = s.replace('full.prmtop', replacement)
            with open(fpath, "w") as f:
              f.write(s)
      # Change to simulation directory
      os.chdir(guest)

    # Find anchors
    with open('disang.rest', 'r') as f:
        data = f.readline().split()    
        L1 = data[6].strip()   
        L2 = data[7].strip()   
        L3 = data[8].strip()   


    # Get number of atoms in vacuum
    with open('./vac.pdb') as myfile:
        data = myfile.readlines()
        vac_atoms = data[-3][6:11].strip()

    # Create minimization and NPT equilibration files for big box and small ligand box 
    if comp != 'c' and comp != 'r':
      with open("../amber_files/mini.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open("../amber_files/therm1.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open("../amber_files/therm2.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)))
    else:
      with open("../amber_files/mini-lig.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line)
      with open("../amber_files/therm1-lig.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line)
      with open("../amber_files/therm2-lig.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line.replace('_temperature_', str(temperature)))


    # Create gradual release files for equilibrium
    if (stage == 'equil'):
      for i in range(0, num_sim):
        with open('../amber_files/mdin-equil', "rt") as fin:
          with open("./mdin-%02d" %int(i), "wt") as fout:
            if i == (num_sim-1):
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang%02d' %int(i)))
            else:
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang%02d' %int(i)))

    # Create free energy files
    if (stage == 'fe'):
      if (comp != 'c' and comp != 'r'):
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-rest', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang'))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))
      else:
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-lig', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang'))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))

    # Create running scripts for local and server
    if (stage == 'fe'):
      if (comp != 'c' and comp != 'r'): 
        with open('../run_files/local-'+stage+'.bash', "rt") as fin:
          with open("./run-local.bash", "wt") as fout:
            for line in fin:
              fout.write(line)
        with open('../run_files/PBS-'+stage, "rt") as fin:
          with open("./PBS-run", "wt") as fout:
            for line in fin:
              fout.write(line.replace('STAGE', guest).replace('POSE', '%s%02d' %(comp, int(win))))
      else:
        with open('../run_files/local-lig.bash', "rt") as fin:
          with open("./run-local.bash", "wt") as fout:
            for line in fin:
              fout.write(line)
        with open('../run_files/PBS-lig', "rt") as fin:
          with open("./PBS-run", "wt") as fout:
            for line in fin:
              fout.write(line.replace('STAGE', guest).replace('POSE', '%s%02d' %(comp, int(win))))
    else:
      with open('../run_files/local-'+stage+'.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line.replace('RANGE', str(rng)))
      with open('../run_files/PBS-'+stage, "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', stage).replace('POSE', guest).replace('RANGE', str(rng)))

    os.chdir('../')


def dec_files(temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, weight, lambdas, ntwx):

    # Find anchors
    with open('disang.rest', 'r') as f:
        data = f.readline().split()    
        L1 = data[6].strip()   
        L2 = data[7].strip()   
        L3 = data[8].strip()   

    # Get number of atoms in vacuum
    with open('./vac.pdb') as myfile:
        data = myfile.readlines()
        vac_atoms = data[-3][6:11].strip()


    if (comp == 'v'):
      for i in range(0, num_sim+1):
        with open('./vac.pdb') as myfile:
          data = myfile.readlines()
          mk2 = int(data[-3][22:26].strip())
          mk1 = int(mk2 - 1)
        with open('../amber_files/mdin-lj', "rt") as fin:
          with open("./mdin-%02d" %int(i), "wt") as fout:
            if i == 1 or i == 0:
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
            else:
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
        mdin = open("./mdin-%02d" %int(i), 'a')
        mdin.write('  mbar_states = %02d\n' %len(lambdas))
        mdin.write('  mbar_lambda = ')
        for i in range(0, len(lambdas)):
          mdin.write(' %6.5f,' %(lambdas[i]))
        mdin.write('\n')
        mdin.write('  infe = 1,\n')
        mdin.write(' /\n')
        mdin.write(' &pmd \n')
        mdin.write(' output_file = \'cmass.txt\' \n')
        mdin.write(' output_freq = %02d \n' % int(ntwx))
        mdin.write(' cv_file = \'cv.in\' \n')
        mdin.write(' /\n')
        mdin.write(' &wt type = \'END\' , /\n')
        mdin.write('DISANG=disang.rest\n')
        mdin.write('LISTOUT=POUT\n')

    if (comp == 'e'):
      for i in range(0, num_sim+1):
        with open('./vac.pdb') as myfile:
          data = myfile.readlines()
          mk4 = int(data[-3][22:26].strip())
          mk3 = int(mk4 - 1)
          mk2 = int(mk4 - 2)
          mk1 = int(mk4 - 3)
        with open('../amber_files/mdin-ch', "rt") as fin:
          with open("./mdin-%02d" %int(i), "wt") as fout:
            if i == 1 or i == 0:
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
            else:
              for line in fin:
                fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
        mdin = open("./mdin-%02d" %int(i), 'a')
        mdin.write('  mbar_states = %02d\n' %len(lambdas))
        mdin.write('  mbar_lambda = ')
        for i in range(0, len(lambdas)):
          mdin.write(' %6.5f,' %(lambdas[i]))
        mdin.write('\n')
        mdin.write('  infe = 1,\n')
        mdin.write(' /\n')
        mdin.write(' &pmd \n')
        mdin.write(' output_file = \'cmass.txt\' \n')
        mdin.write(' output_freq = %02d \n' % int(ntwx))
        mdin.write(' cv_file = \'cv.in\' \n')
        mdin.write(' /\n')
        mdin.write(' &wt type = \'END\' , /\n')
        mdin.write('DISANG=disang.rest\n')
        mdin.write('LISTOUT=POUT\n')


    if (comp == 'v'):
      # Create heating and NPT equilibration files for vdw decoupling
      with open("../amber_files/eqnpt-lj.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      with open("../amber_files/heat-lj.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      # Create running scripts for local and server
      with open('../run_files/local-sdr.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-sdr', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', guest).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'e'):
      # Create heating and NPT equilibration files for charge decoupling
      with open("../amber_files/eqnpt-ch.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
      with open("../amber_files/heat-ch.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
      # Create running scripts for local and server
      with open('../run_files/local-sdr.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-sdr', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', guest).replace('POSE', '%s%02d' %(comp, int(win))))

    os.chdir('../')
    

