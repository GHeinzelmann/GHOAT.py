#!/usr/bin/env python3
import datetime as dt
import glob as glob
import os as os
import re as re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import scripts 

def build_equil(guest, host, mol, H1, H2, H3, min_adis, max_adis, l1_range, amber_ff, final_host_num, guest_charge, sdr_dist):


    # Create equilibrium directory
    if not os.path.exists('equil'):
      os.makedirs('equil')
    os.chdir('equil')
    if os.path.exists('./build_files'):
      shutil.rmtree('./build_files')
    try:
      shutil.copytree('../build_files', './build_files')
    # Directories are the same
    except shutil.Error as e:
      print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
      print('Directory not copied. Error: %s' % e)
    os.chdir('build_files')

    # Copy structure and separate guest and host
    shutil.copy('../../structures/%s-%s.pdb' %(host.lower(), guest.lower()), './')
    with open("separate-ini.tcl", "rt") as fin:
      with open("separate.tcl", "wt") as fout:
        for line in fin:
          fout.write(line.replace('MMM', mol).replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('FIRST','1').replace('LAST', str(final_host_num)))
    sp.call('vmd -dispdev text -e separate.tcl', shell=True)

 
    h1_resid = H1.split('@')[0][1:]
    h2_resid = H2.split('@')[0][1:]
    h3_resid = H3.split('@')[0][1:]

    h1_atom = H1.split('@')[1]
    h2_atom = H2.split('@')[1]
    h3_atom = H3.split('@')[1]
   
    print('Receptor anchors:')
    print(H1)
    print(H2)
    print(H3)

    # Replace names in initial files and VMD scripts
    with open("prep-ini.tcl", "rt") as fin:
      with open("prep.tcl", "wt") as fout:
        for line in fin:
          fout.write(line.replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('MMM', mol).replace('NN1', h1_atom).replace('H1A', h1_resid).replace('NN2', h2_atom).replace('H2A', h2_resid).replace('NN3', h3_atom).replace('H3A', h3_resid).replace('FIRST','1').replace('LAST', str(final_host_num)).replace('STAGE','equil').replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('RANG','%4.2f' %l1_range).replace('SDRD', str(sdr_dist)))

    # Get parameters for host
    shutil.copy('../../parameters/%s.mol2' %(host.lower()), './')  # AM1-BCC charges problematic for cyclic hosts, using provided charges only
    if not os.path.exists('../../parameters/%s.frcmod' %host.lower()):
      if amber_ff == 'gaff':
        sp.call('parmchk2 -i '+host.lower()+'.mol2 -f mol2 -o '+host.lower()+'.frcmod -s 1', shell=True)
      elif amber_ff == 'gaff2':
        sp.call('parmchk2 -i '+host.lower()+'.mol2 -f mol2 -o '+host.lower()+'.frcmod -s 2', shell=True)
    else:
      shutil.copy('../../parameters/%s.frcmod' %(host.lower()), './')

    # Get parameters for guest
    if not os.path.exists('../../parameters/%s.mol2' %guest.lower()):
      print('Antechamber parameters command: antechamber -i '+guest.lower()+'.pdb -fi pdb -o '+guest.lower()+'.mol2 -fo mol2 -c bcc -s 2 -at '+amber_ff.lower()+' -nc %s' % guest_charge)
      sp.call('antechamber -i '+guest.lower()+'.pdb -fi pdb -o '+guest.lower()+'.mol2 -fo mol2 -c bcc -s 2 -at '+amber_ff.lower()+' -nc %s' % guest_charge, shell=True)
    else:
      shutil.copy('../../parameters/%s.mol2' %(guest.lower()), './') # Provided charges for guest
    if not os.path.exists('../../parameters/%s.frcmod' %guest.lower()):
      if amber_ff == 'gaff':
        sp.call('parmchk2 -i '+guest.lower()+'.mol2 -f mol2 -o '+guest.lower()+'.frcmod -s 1', shell=True)
      elif amber_ff == 'gaff2':
        sp.call('parmchk2 -i '+guest.lower()+'.mol2 -f mol2 -o '+guest.lower()+'.frcmod -s 2', shell=True)
    else:
      shutil.copy('../../parameters/%s.frcmod' %(guest.lower()), './')

    # Put complex in AMBER format and find ligand anchor atoms
    sp.call('vmd -dispdev text -e prep.tcl', shell=True)

    # Save parameters in ff folder
    if not os.path.exists('../ff/'):
      os.makedirs('../ff/')
    shutil.copy('./%s.mol2' %(host.lower()), '../ff/')
    shutil.copy('./%s.frcmod' %(host.lower()), '../ff/')
    shutil.copy('./%s.mol2' %(guest.lower()), '../ff/')
    shutil.copy('./%s.frcmod' %(guest.lower()), '../ff/')

    # Check size of anchor file 
    anchor_file = 'anchors.txt'
    if os.stat(anchor_file).st_size == 0:
      os.chdir('../')
      return 'anch1'
    f = open(anchor_file, 'r')
    for line in f:
      splitdata = line.split()
      if len(splitdata) < 3:
        os.rename('./anchors.txt', 'anchors-'+guest.lower()+'.txt')
        os.chdir('../')
        return 'anch2'

    os.rename('./anchors.txt', 'anchors-'+guest.lower()+'.txt')
    os.chdir('../')

    # Create simulation directory
    if not os.path.exists(guest):
      os.makedirs(guest)
    os.chdir(guest)
    shutil.copy('../build_files/%s-%s-aligned.pdb' %(host.lower(), guest.lower()), './build-ini.pdb')
    shutil.copy('../build_files/%s.pdb' %(guest.lower()), './')
    shutil.copy('../build_files/%s.pdb' %(host.lower()), './')
    shutil.copy('../build_files/anchors-%s.txt' %(guest.lower()), './')
    shutil.copy('../build_files/dum1.pdb', './')
    shutil.copy('../build_files/dum.mol2', './')
    shutil.copy('../build_files/dum.frcmod', './')

    dum_coords = []
    recep_coords = []
    lig_coords = []
    dum_atomlist = []
    lig_atomlist = []
    recep_atomlist = []
    dum_atom = 0
    lig_atom = 0
    recep_atom = 0
    total_atom = 0
    resid_lig = 0
    resname_lig = mol
    resname_list = []
    resid_list = []

    # Read coordinates for dummy atoms
    for i in range(1, 2):
      shutil.copy('../build_files/dum'+str(i)+'.pdb', './')
      with open('dum'+str(i)+'.pdb') as dum_in:
        lines = (line.rstrip() for line in dum_in)
        lines = list(line for line in lines if line)
        dum_coords.append((float(lines[1][30:38].strip()), float(lines[1][38:46].strip()), float(lines[1][46:54].strip())))
        dum_atomlist.append(lines[1][12:16].strip())
        resname_list.append(lines[1][17:20].strip())
        resid_list.append(float(lines[1][22:26].strip()))
        dum_atom += 1
        total_atom += 1


    # Read coordinates from aligned system
    with open('build-ini.pdb') as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line) # Non-blank lines in a list   

    # Count atoms of receptor and ligand
    for i in range(0, len(lines)):
      if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
              if (lines[i][17:20].strip() != mol) and (lines[i][17:20].strip() != 'DUM'):
                 recep_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                 recep_atomlist.append(lines[i][12:16].strip())
                 resname_list.append(lines[i][17:20].strip())
                 resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                 recep_last = int(lines[i][22:26].strip())
                 recep_atom += 1
                 total_atom += 1
              elif lines[i][17:20].strip() == mol:
                 lig_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                 lig_atomlist.append(lines[i][12:16].strip())
                 resname_list.append(lines[i][17:20].strip())
                 resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                 lig_atom += 1
                 total_atom += 1

    coords = dum_coords + recep_coords + lig_coords
    atom_namelist = dum_atomlist + recep_atomlist + lig_atomlist
    lig_resid = recep_last + dum_atom + 1

    # Write the new pdb file

    build_file = open('build.pdb', 'w')

    # Positions for the dummy atoms 
    for i in range(0, 1):
       build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
       build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))
       build_file.write('%6.2f%6.2f\n'%(0, 0))
       build_file.write('TER\n')

    # Positions of the receptor atoms
    for i in range(dum_atom , dum_atom + recep_atom):
        build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
        build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))

        build_file.write('%6.2f%6.2f\n'%(0, 0))
    build_file.write('TER\n')

    # Positions of the ligand atoms
    for i in range(dum_atom + recep_atom, total_atom):
        build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],mol, float(lig_resid)))
        build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]),float(coords[i][2])))

        build_file.write('%6.2f%6.2f\n'%(0, 0))

    build_file.write('TER\n')
    build_file.write('END\n')
    build_file.close()

    os.chdir('../')

    return 'all'

    
def build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist):


    # Get files or finding new anchors and building some systems

    if not os.path.exists('../build_files'):
      try:
        shutil.copytree('../../../build_files', '../build_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      os.chdir('../build_files') 

      # Replace names in initial files and VMD scripts
      h1_resid = H1.split('@')[0][1:]
      h2_resid = H2.split('@')[0][1:]
      h3_resid = H3.split('@')[0][1:]

      h1_atom = H1.split('@')[1]
      h2_atom = H2.split('@')[1]
      h3_atom = H3.split('@')[1]
     
      print('Receptor anchors:')
      print(H1)
      print(H2)
      print(H3)
    
      with open("prep-ini.tcl", "rt") as fin:
        with open("prep.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('MMM', mol).replace('NN1', h1_atom).replace('H1A', str(int(h1_resid)+1)).replace('NN2', h2_atom).replace('H2A', str(int(h2_resid)+1)).replace('NN3', h3_atom).replace('H3A', str(int(h3_resid)+1)).replace('FIRST','2').replace('LAST', str(int(final_host_num)+1)).replace('STAGE','equil').replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('RANG','%4.2f' %l1_range).replace('SDRD', str(sdr_dist)))
      with open("separate-ini.tcl", "rt") as fin:
        with open("separate.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('MMM', mol).replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('FIRST','2').replace('LAST', str(int(final_host_num)+1)))

      # Get parameters from equilibrium
      if not os.path.exists('../ff'):
        os.makedirs('../ff')
      shutil.copy('../../../equil/ff/%s.mol2' %(guest.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(guest.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.mol2' %(host.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(host.lower()), '../ff/')

      # Get parameter file and final state from equilibrium
      for file in glob.glob('../../../equil/%s/full*.prmtop' %guest.lower()):
        shutil.copy(file, './')
      for file in glob.glob('../../../equil/%s/vac*' %guest.lower()):
        shutil.copy(file, './')
      shutil.copy('../../../equil/%s/md%02d.rst7' %(guest.lower(), fwin), './')
      sp.call('cpptraj -p full.prmtop -y md%02d.rst7 -x %s-%s.pdb > cpptraj.log' %(fwin, host.lower(), guest.lower()), shell=True)
      sp.call('vmd -dispdev text -e separate.tcl', shell=True)
      sp.call('vmd -dispdev text -e prep.tcl', shell=True)
      os.rename('./anchors.txt', 'anchors-'+guest.lower()+'.txt')

      os.chdir('../rest/')

      # Check size of anchor file 
      anchor_file = '../build_files/anchors-'+guest.lower()+'.txt'
      if os.stat(anchor_file).st_size == 0:
        return 'anch1'
      f = open(anchor_file, 'r')
      for line in f:
        splitdata = line.split()
        if len(splitdata) < 3:
          return 'anch2'
      
    # Copy and replace simulation files for the first window
    if int(win) == 0:
      if os.path.exists('amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../../../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_amber_ff_', amber_ff)
          with open(fpath, "w") as f:
            f.write(s)

    if not os.path.exists('run_files'):
      try:
        shutil.copytree('../../../run_files', './run_files')
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

     
    if (comp == 'a' or comp == 'l' or comp == 't'):
      # Create window directory
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      # Copy a few files and define new reference state
      if int(win) == 0:
        for file in glob.glob('../../build_files/full*.prmtop'):
          shutil.copy(file, './')
        for file in glob.glob('../../build_files/vac*'):
          shutil.copy(file, './')
        shutil.copy('../../build_files/md%02d.rst7' %fwin, './md00.rst7')
        shutil.copy('../../build_files/anchors-'+guest+'.txt', './')
        sp.call('cpptraj -p full.prmtop -y md00.rst7 -x full.rst7 > cpptraj1.log', shell=True)
        shutil.copy('./full.rst7', './full.inpcrd')
        sp.call('cpptraj -p full.prmtop -y md00.rst7 -x full.pdb > cpptraj2.log', shell=True)
      else:
        for file in glob.glob('../%s00/*' %comp):
          shutil.copy(file, './')
    elif comp == 'c':
      # Copy files to c00 to create new box for ligand and copy to the different windows 
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      if int(win) == 0:
        shutil.copy('../../build_files/'+guest+'.pdb', './')
        shutil.copy('../../build_files/'+host+'-'+guest+'.pdb', './')
        shutil.copy('../../build_files/anchors-'+guest+'.txt', './')
      else:
        for file in glob.glob('../c00/*'):
          shutil.copy(file, './')
    elif comp == 'r':
      # Copy files to r00 to create new box for host and copy to the different windows 
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      if int(win) == 0:
        # Get files and parameters for building
        shutil.copy('../../build_files/'+guest+'.pdb', './')
        shutil.copy('../../build_files/'+host+'.pdb', './')
        shutil.copy('../../build_files/'+host+'.pdb', './build.pdb')
        shutil.copy('../../build_files/'+host+'-'+guest+'.pdb', './')
        shutil.copy('../../build_files/anchors-'+guest+'.txt', './')
        shutil.copy('../../build_files/dum.mol2', './')
        shutil.copy('../../build_files/dum.frcmod', './')
        for file in glob.glob('../../ff/*.frcmod'):
          shutil.copy(file, './')
        for file in glob.glob('../../ff/*.mol2'):
          shutil.copy(file, './')
      else:
        for file in glob.glob('../r00/*'):
          shutil.copy(file, './')

def build_dec(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist):



    # Get files or finding new anchors and building some systems

    if not os.path.exists('../build_files'):
      try:
        shutil.copytree('../../../build_files', '../build_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      os.chdir('../build_files')

      # Replace names in initial files and VMD scripts
      h1_resid = H1.split('@')[0][1:]
      h2_resid = H2.split('@')[0][1:]
      h3_resid = H3.split('@')[0][1:]

      h1_atom = H1.split('@')[1]
      h2_atom = H2.split('@')[1]
      h3_atom = H3.split('@')[1]

      print('Receptor anchors:')
      print(H1)
      print(H2)
      print(H3)

      with open("prep-ini.tcl", "rt") as fin:
        with open("prep.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('MMM', mol).replace('NN1', h1_atom).replace('H1A', str(int(h1_resid)+1)).replace('NN2', h2_atom).replace('H2A', str(int(h2_resid)+1)).replace('NN3', h3_atom).replace('H3A', str(int(h3_resid)+1)).replace('FIRST','2').replace('LAST', str(int(final_host_num)+1)).replace('STAGE','equil').replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('RANG','%4.2f' %l1_range).replace('SDRD', str(sdr_dist)))
      with open("separate-ini.tcl", "rt") as fin:
        with open("separate.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('MMM', mol).replace('hhhh', host.lower()).replace('gggg', guest.lower()).replace('FIRST','2').replace('LAST', str(int(final_host_num)+1)))

      # Get parameters from equilibrium
      if not os.path.exists('../ff'):
        os.makedirs('../ff')
      shutil.copy('../../../equil/ff/%s.mol2' %(guest.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(guest.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.mol2' %(host.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(host.lower()), '../ff/')


      # Get parameter file and final state from equilibrium
      for file in glob.glob('../../../equil/%s/full*.prmtop' %guest.lower()):
        shutil.copy(file, './')
      for file in glob.glob('../../../equil/%s/vac*' %guest.lower()):
        shutil.copy(file, './')
      shutil.copy('../../../equil/%s/md%02d.rst7' %(guest.lower(), fwin), './')
      sp.call('cpptraj -p full.prmtop -y md%02d.rst7 -x %s-%s.pdb > cpptraj.log' %(fwin, host.lower(), guest.lower()), shell=True)
      sp.call('vmd -dispdev text -e separate.tcl', shell=True)
      sp.call('vmd -dispdev text -e prep.tcl', shell=True)
      os.rename('./anchors.txt', 'anchors-'+guest.lower()+'.txt')

      os.chdir('../sdr/')

    # Copy and replace simulation files for the first window
    if int(win) == 0:
      if os.path.exists('./amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../../../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_amber_ff_', amber_ff)
          with open(fpath, "w") as f:
            f.write(s)

    if not os.path.exists('./run_files'):
      try:
        shutil.copytree('../../../run_files', './run_files')
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

    if not os.path.exists('./%s%02d' %(comp, int(win))):
      os.makedirs('./%s%02d' %(comp, int(win)))
    os.chdir('./%s%02d' %(comp, int(win)))
    if int(win) == 0:
      for file in glob.glob('../../build_files/vac*'):
        shutil.copy(file, './')
      shutil.copy('../../build_files/'+guest+'.pdb', './')
      shutil.copy('../../build_files/'+host+'-'+guest+'.pdb', './')
      shutil.copy('../../build_files/'+host+'-'+guest+'-aligned.pdb', './build-ini.pdb')
      shutil.copy('../../build_files/anchors-'+guest+'.txt', './')
      shutil.copy('../../build_files/dum.frcmod', './')
      shutil.copy('../../build_files/dum.mol2', './')
      for file in glob.glob('../../ff/*.frcmod'):
        shutil.copy(file, './')
      for file in glob.glob('../../ff/*.mol2'):
        shutil.copy(file, './')

      dum_coords = []
      recep_coords = []
      lig_coords = []
      dum_atomlist = []
      lig_atomlist = []
      recep_atomlist = []
      dum_atom = 0
      lig_atom = 0
      recep_atom = 0
      total_atom = 0
      resid_lig = 0
      resname_lig = mol
      resname_list = []
      resid_list = []

      # Read coordinates for dummy atoms
      for i in range(1, 3):
        shutil.copy('../../build_files/dum'+str(i)+'.pdb', './')
        with open('dum'+str(i)+'.pdb') as dum_in:
          lines = (line.rstrip() for line in dum_in)
          lines = list(line for line in lines if line)
          dum_coords.append((float(lines[1][30:38].strip()), float(lines[1][38:46].strip()), float(lines[1][46:54].strip())))
          dum_atomlist.append(lines[1][12:16].strip())
          resname_list.append(lines[1][17:20].strip())
          resid_list.append(float(lines[1][22:26].strip()))
          dum_atom += 1
          total_atom += 1


      # Read coordinates from aligned system
      with open('build-ini.pdb') as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list   

      # Count atoms of receptor and ligand
      for i in range(0, len(lines)):
        if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                if (lines[i][17:20].strip() != mol) and (lines[i][17:20].strip() != 'DUM'):
                   recep_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                   recep_atomlist.append(lines[i][12:16].strip())
                   resname_list.append(lines[i][17:20].strip())
                   resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                   recep_last = int(lines[i][22:26].strip())
                   recep_atom += 1
                   total_atom += 1
                elif lines[i][17:20].strip() == mol:
                   lig_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                   lig_atomlist.append(lines[i][12:16].strip())
                   resname_list.append(lines[i][17:20].strip())
                   resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                   lig_atom += 1
                   total_atom += 1

      coords = dum_coords + recep_coords + lig_coords
      atom_namelist = dum_atomlist + recep_atomlist + lig_atomlist
      lig_resid = recep_last + dum_atom + 1

      # Write the new pdb file

      build_file = open('build.pdb', 'w')

      # Positions for the dummy atoms 
      for i in range(0, 2):
         build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
         build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))
         build_file.write('%6.2f%6.2f\n'%(0, 0))
         build_file.write('TER\n')

      # Positions of the receptor atoms
      for i in range(dum_atom , dum_atom + recep_atom):
          build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
          build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))

          build_file.write('%6.2f%6.2f\n'%(0, 0))
      build_file.write('TER\n')

      # Positions of the ligand atoms
      for i in range(dum_atom + recep_atom, total_atom):
          build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],mol, float(lig_resid)))
          build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]),float(coords[i][2])))

          build_file.write('%6.2f%6.2f\n'%(0, 0))

      build_file.write('TER\n')

      # Extra guests for decoupling


      build_file = open('build.pdb', 'a')
      if (comp == 'e'):
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+1)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2])))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        build_file.write('TER\n')
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+2)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        build_file.write('TER\n')
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+3)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        print('Creating new system for decharging...')
      else:
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid + 1)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        print('Creating new system for vdw decoupling...')
      build_file.write('TER\n')
      build_file.write('END\n')
      build_file.close()
    else:
      for file in glob.glob('../'+comp+'00/*'):
        shutil.copy(file, './')


def create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt):
    
    # Copy and replace simulation files
    if stage != 'fe':
      if os.path.exists('amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_amber_ff_', amber_ff)
          with open(fpath, "w") as f:
            f.write(s)
      os.chdir(guest)

    # Copy tleap files that are used for restraint generation and analysis
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap_vac.in')
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap_vac_guest.in')
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap_vac_host.in')
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap.in')

    # Copy host and guest parameter files
    for file in glob.glob('../ff/*'):
       shutil.copy(file, './')

    # Append tleap file for vacuum
    tleap_vac = open('tleap_vac.in', 'a')
    tleap_vac.write('# Load the guest parameters\n')        
    tleap_vac.write('loadamberparams %s.frcmod\n'%(guest.lower()))
    tleap_vac.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), guest.lower()))
    tleap_vac.write('# Load the host parameters\n')        
    tleap_vac.write('loadamberparams %s.frcmod\n'%(host.lower()))
    tleap_vac.write('%s = loadmol2 %s.mol2\n\n'%(hmol.upper(), host.lower()))
    tleap_vac.write('model = loadpdb build.pdb\n\n')
    tleap_vac.write('check model\n')
    tleap_vac.write('savepdb model vac.pdb\n')
    tleap_vac.write('saveamberparm model vac.prmtop vac.inpcrd\n')
    tleap_vac.write('quit\n')
    tleap_vac.close()

    # Append tleap file for guest only
    tleap_vac_ligand = open('tleap_vac_guest.in', 'a')
    tleap_vac_ligand.write('# Load the guest parameters\n')        
    tleap_vac_ligand.write('loadamberparams %s.frcmod\n'%(guest.lower()))
    tleap_vac_ligand.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), guest.lower()))
    tleap_vac_ligand.write('model = loadpdb %s.pdb\n\n' %(guest.lower()))
    tleap_vac_ligand.write('check model\n')
    tleap_vac_ligand.write('savepdb model vac_guest.pdb\n')
    tleap_vac_ligand.write('saveamberparm model vac_guest.prmtop vac_guest.inpcrd\n')
    tleap_vac_ligand.write('quit\n')
    tleap_vac_ligand.close()

    # Append tleap file for host only
    tleap_vac_ligand = open('tleap_vac_host.in', 'a')
    tleap_vac_ligand.write('# Load the host parameters\n')        
    tleap_vac_ligand.write('loadamberparams %s.frcmod\n'%(host.lower()))
    tleap_vac_ligand.write('%s = loadmol2 %s.mol2\n\n'%(hmol.upper(), host.lower()))
    tleap_vac_ligand.write('model = loadpdb %s.pdb\n\n' %(host.lower()))
    tleap_vac_ligand.write('check model\n')
    tleap_vac_ligand.write('savepdb model vac_host.pdb\n')
    tleap_vac_ligand.write('saveamberparm model vac_host.prmtop vac_host.inpcrd\n')
    tleap_vac_ligand.write('quit\n')
    tleap_vac_ligand.close()


    # Generate complex in vacuum
    p = sp.call('tleap -s -f tleap_vac.in > tleap_vac.log', shell=True)

    # Generate individual structures in vacuum
    p = sp.call('tleap -s -f tleap_vac_guest.in > tleap_vac_guest.log', shell=True)
    p = sp.call('tleap -s -f tleap_vac_host.in > tleap_vac_host.log', shell=True)

    # Find out how many cations/anions are needed for neutralization
    neu_cat = 0
    neu_ani = 0
    f = open('tleap_vac.log', 'r')
    for line in f:
        if "The unperturbed charge of the unit" in line:
            splitline = line.split()
            if float(splitline[6].strip('\'\",.:;#()][')) < 0:
                neu_cat = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
            elif float(splitline[6].strip('\'\",.:;#()][')) > 0:
                neu_ani = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
    f.close()
     
    # Get ligand removed charge when doing LJ calculations
    lig_cat = 0
    lig_ani = 0
    f = open('tleap_vac_guest.log', 'r')
    for line in f:
        if "The unperturbed charge of the unit" in line:
            splitline = line.split()
            if float(splitline[6].strip('\'\",.:;#()][')) < 0:
                lig_cat = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
            elif float(splitline[6].strip('\'\",.:;#()][')) > 0:
                lig_ani = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
    f.close()
     
    # Adjust ions for LJ and electrostatic Calculations (avoid neutralizing plasma)
    if comp == 'v':
      charge_neut = neu_cat - neu_ani - 2*lig_cat + 2*lig_ani
      neu_cat = 0
      neu_ani = 0
      if charge_neut > 0:
        neu_cat = abs(charge_neut)
      if charge_neut < 0:
        neu_ani = abs(charge_neut)
    if comp == 'e':
      charge_neut = neu_cat - neu_ani - 3*lig_cat + 3*lig_ani
      neu_cat = 0
      neu_ani = 0
#      print (charge_neut)
      if charge_neut > 0:
        neu_cat = abs(charge_neut)
      if charge_neut < 0:
        neu_ani = abs(charge_neut)
 
    # Define volume density for different water models
    if water_model == 'TIP3P':
       water_box = water_model.upper()+'BOX'
       ratio = 0.0576
    elif water_model == 'SPCE': 
       water_box = 'SPCBOX'
       ratio = 0.0576
    elif water_model == 'TIP4PEW': 
       water_box = water_model.upper()+'BOX'
       ratio = 0.0573

    # Update target number of residues according to the ion definitions 
    if (neut == 'no'):
      target_num = int(num_waters - neu_cat + neu_ani + 2*int(ion_def[2])) 
    elif (neut == 'yes'):
      target_num = int(num_waters + neu_cat + neu_ani)
    
    
    # Number of cations and anions   
    num_cat = ion_def[2]
    num_ani = ion_def[2] - neu_cat + neu_ani
    

    # Create the first box guess to get the initial number of waters and cross sectional area
    buff = 50.0  
    scripts.write_tleap(mol, hmol, guest, host, water_model, water_box, buff, buffer_x, buffer_y)
    num_added = scripts.check_tleap()
    cross_area = scripts.cross_sectional_area()

    # Define a few parameters for solvation iteration
    count = 0
    max_count = 10
    rem_limit = 16
    factor = 1
    ind = 0.90   
    buff_diff = 1.0  

    # Iterate to get the correct number of waters
    while num_added != target_num:
        count += 1
        if count > max_count:
        # Try different parameters
             rem_limit += 4
             if ind > 0.5:
               ind = ind - 0.02
             else:
               ind = 0.90
             factor = 1
             max_count = max_count + 10
        tleap_remove = None
        # Manually remove waters if inside removal limit
        if num_added > target_num and (num_added - target_num) < rem_limit:
            difference = num_added - target_num
            tleap_remove = [target_num + 1 + i for i in range(difference)]
            scripts.write_tleap(mol, hmol, guest, host, water_model, water_box, buff, buffer_x, buffer_y, tleap_remove)
            scripts.check_tleap()
            break
        # Set new buffer size based on chosen water density
        res_diff = num_added - target_num - (rem_limit/2)
        buff_diff = res_diff/(ratio*cross_area)
        buff -= (buff_diff * factor)
        if buff < 0:
           print ('Not enough water molecules to fill the system in the z direction, please increase the number of water molecules')
           sys.exit(1)
        # Set relaxation factor  
        factor = ind * factor
        # Get number of waters
        scripts.write_tleap(mol, hmol, guest, host, water_model, water_box, buff, buffer_x, buffer_y)
        num_added = scripts.check_tleap()       
 
    # Write the final tleap file with the correct system size and removed water molecules
    shutil.copy('tleap.in', 'tleap_solvate.in')
    tleap_solvate = open('tleap_solvate.in', 'a')
    tleap_solvate.write('# Load the guest parameters\n')        
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(guest.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), guest.lower()))
    tleap_solvate.write('# Load the host parameters\n')        
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(host.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(hmol.upper(), host.lower()))
    tleap_solvate.write('model = loadpdb build.pdb\n\n')
    tleap_solvate.write('# Load the water and jc ion parameters\n')        
    tleap_solvate.write('source leaprc.water.%s\n'%(water_model.lower()))
    tleap_solvate.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tleap_solvate.write('# Create water box with chosen model\n')
    tleap_solvate.write('solvatebox model ' + water_box + ' {'+ str(buffer_x) +' '+ str(buffer_y) +' '+ str(buff) +'}\n\n')
    if tleap_remove is not None:
        tleap_solvate.write('# Remove a few waters manually\n')
        for water in tleap_remove:
            tleap_solvate.write('remove model model.%s\n' % water)
        tleap_solvate.write('\n')
    # Ionize/neutralize system
    if (neut == 'no'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[0], num_cat))
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[1], num_ani))
    elif (neut == 'yes'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        if neu_cat != 0:
          tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[0], neu_cat))
        if neu_ani != 0:
          tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[1], neu_ani))
    tleap_solvate.write('\n')
    tleap_solvate.write('desc model\n')
    tleap_solvate.write('savepdb model full.pdb\n')
    tleap_solvate.write('saveamberparm model full.prmtop full.inpcrd\n')
    tleap_solvate.write('quit')
    tleap_solvate.close()
    p = sp.call('tleap -s -f tleap_solvate.in > tleap_solvate.log', shell=True)
    
    f = open('tleap_solvate.log', 'r')
    for line in f:
        if "Could not open file" in line:
           print ('WARNING!!!')
           print (line)
           sys.exit(1)
        if "WARNING: The unperturbed charge of the unit:" in line:
           print (line)
           print ('The system is not neutralized properly after solvation')
        if "addIonsRand: Argument #2 is type String must be of type: [unit]" in line:
           print('Aborted.The ion types specified in the input file could be wrong.')              
           print('Please check the tleap_solvate.log file, and the ion types specified in the input file.\n')
           sys.exit(1)
    f.close()


    # Apply hydrogen mass repartitioning
    print('Applying mass repartitioning...')
    shutil.copy('../amber_files/parmed-hmr.in', './')
    sp.call('parmed -O -n -i parmed-hmr.in > parmed-hmr.log', shell=True)

    if stage != 'fe':
      os.chdir('../')


def guest_box(guest, mol, lig_buffer, water_model, neut, ion_lig, comp, amber_ff):

    # Define volume density for different water models
    if water_model == 'TIP3P':
       water_box = water_model.upper()+'BOX'
    elif water_model == 'SPCE':
       water_box = 'SPCBOX'
    elif water_model == 'TIP4PEW':
       water_box = water_model.upper()+'BOX'

    # Copy ligand parameter files
    for file in glob.glob('../../ff/%s.*' %guest.lower()):
        shutil.copy(file, './')

    # Write and run tleap file
    tleap_solvate = open('tleap_solvate.in', 'a')
    tleap_solvate.write('source leaprc.'+amber_ff+'\n\n')
    tleap_solvate.write('# Load the ligand parameters\n')
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(guest.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), guest.lower()))
    tleap_solvate.write('model = loadpdb %s.pdb\n\n' %(guest.lower()))
    tleap_solvate.write('# Load the water and jc ion parameters\n')
    tleap_solvate.write('source leaprc.water.%s\n'%(water_model.lower()))
    tleap_solvate.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tleap_solvate.write('check model\n')
    tleap_solvate.write('savepdb model vac.pdb\n')
    tleap_solvate.write('saveamberparm model vac.prmtop vac.inpcrd\n\n')
    tleap_solvate.write('# Create water box with chosen model\n')
    tleap_solvate.write('solvatebox model ' + water_box + ' '+str(lig_buffer)+'\n\n')
    if (neut == 'no'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_lig[0], ion_lig[2]))
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_lig[1]))
    elif (neut == 'yes'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_lig[0]))
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_lig[1]))
    tleap_solvate.write('\n')
    tleap_solvate.write('desc model\n')
    tleap_solvate.write('savepdb model full.pdb\n')
    tleap_solvate.write('saveamberparm model full.prmtop full.inpcrd\n')
    tleap_solvate.write('quit\n')
    tleap_solvate.close()
    p = sp.call('tleap -s -f tleap_solvate.in > tleap_solvate.log', shell=True)

    # Apply hydrogen mass repartitioning
    print('Applying mass repartitioning...')
    shutil.copy('../amber_files/parmed-hmr.in', './')
    sp.call('parmed -O -n -i parmed-hmr.in > parmed-hmr.log', shell=True)

    # Copy a few files for consistency
    shutil.copy('./vac.pdb','./vac_guest.pdb')
    shutil.copy('./vac.prmtop','./vac_guest.prmtop')

