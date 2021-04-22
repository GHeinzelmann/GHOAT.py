#!/usr/bin/env python3
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import build 
from lib import scripts 
from lib import setup
from lib import analysis

ion_def = []
guest_list = []
guest_list_code = []
guest_list_charge = []
guest_def = []
code_def = []
charge_def = []
release_eq = []
attach_rest = []
lambdas = []  
weights = []  
components = []  
aa1_guests = []  
aa2_guests = []  
rec_dihcf_force = 0
rec_discf_force = 0
guest_charge = 0
guest_rot = 'no'
final_host_num = 1

# Read arguments that define input file and stage
if len(sys.argv) < 5:
  scripts.help_message()
  sys.exit(0)
for i in [1, 3]:
  if '-i' == sys.argv[i].lower():
    input_file = sys.argv[i + 1]
  elif '-s' == sys.argv[i].lower():
    stage = sys.argv[i + 1]
  else:
    scripts.help_message()
    sys.exit(1)

# Open input file
with open(input_file) as f_in:       
    # Remove spaces and tabs
    lines = (line.strip(' \t\n\r') for line in f_in)
    lines = list(line for line in lines if line)  # Non-blank lines in a list

for i in range(0, len(lines)):
    # split line using the equal sign, and remove text after #
    if not lines[i][0] == '#':
        lines[i] = lines[i].split('#')[0].split('=')

# Read parameters from input file 
for i in range(0, len(lines)):
    if not lines[i][0] == '#':
        lines[i][0] = lines[i][0].strip().lower()
        lines[i][1] = lines[i][1].strip()
        if lines[i][0] == 'temperature':
            temperature = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'eq_steps1':
            eq_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'eq_steps2':
            eq_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'a_steps1':
            a_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'a_steps2':
            a_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l_steps1':
            l_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l_steps2':
            l_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 't_steps1':
            t_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 't_steps2':
            t_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'c_steps1':
            c_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'c_steps2':
            c_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'r_steps1':
            r_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'r_steps2':
            r_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps1':
            e_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps2':
            e_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'v_steps1':
            v_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'v_steps2':
            v_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'guest_list':
            newline = lines[i][1].strip('\'\"-,.:;#()][').split(',')
            for j in range(0, len(newline)):
               guest_list.append(newline[j].lower())
        elif lines[i][0] == 'guest_list_code':
            newline = lines[i][1].strip('\'\"-,.:;#()][').split(',')
            for j in range(0, len(newline)):
               guest_list_code.append(newline[j])
        elif lines[i][0] == 'guest_list_charge':
            newline = lines[i][1].strip('\'\"-,.:;#()][').split(',')
            for j in range(0, len(newline)):
               guest_list_charge.append(newline[j])
        elif lines[i][0] == 'calc_type':
            calc_type = lines[i][1].lower()
        elif lines[i][0] == 'host':
            host = lines[i][1].lower()
        elif lines[i][0] == 'final_host_num':
            final_host_num = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'h1':
            H1 = lines[i][1]
        elif lines[i][0] == 'h2':
            H2 = lines[i][1]
        elif lines[i][0] == 'h3':
            H3 = lines[i][1]
        elif lines[i][0] == 'host_code':
            hmol = lines[i][1]
        elif lines[i][0] == 'fe_type':
            if lines[i][1].lower() == 'rest':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'sdr':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'all':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'custom':
                fe_type = lines[i][1].lower()
            else:
                print('Free energy type not recognized, please choose all, rest (only restraints) or sdr (simultaneous decoupling/recoupling), or custom')
                sys.exit(1)
        elif lines[i][0] == 'dec_int':
            if lines[i][1].lower() == 'mbar':
                dec_int = lines[i][1].lower()
            elif lines[i][1].lower() == 'ti':
                dec_int = lines[i][1].lower()
            else:
                print('Double decoupling type not recognized, please choose ti or mbar')
                sys.exit(1)
        elif lines[i][0] == 'blocks':
            blocks = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'hmr':
            if lines[i][1].lower() == 'yes':
                hmr = 'yes'
            elif lines[i][1].lower() == 'no':
                hmr = 'no'
            else:
                print('Wrong input! Please use yes or no to indicate whether hydrogen mass repartitioning '
                      'will be used.')
                sys.exit(1)
        elif lines[i][0] == 'water_model':
            if lines[i][1].lower() == 'tip3p':
                water_model = lines[i][1].upper()
            elif lines[i][1].lower() == 'tip4pew':
                water_model = lines[i][1].upper()
            elif lines[i][1].lower() == 'spce':
                water_model = lines[i][1].upper()
            else:
                print('Water model not supported. Please choose TIP3P, TIP4PEW or SPCE')
                sys.exit(1)
        elif lines[i][0] == 'num_waters':
            num_waters = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'neutralize_only':
            if lines[i][1].lower() == 'yes':
                neut = 'yes'
            elif lines[i][1].lower() == 'no':
                neut = 'no'
            else:
                print('Wrong input! Please choose neutralization only or add extra ions')
                sys.exit(1)
        elif lines[i][0] == 'guest_rot':
            if lines[i][1].lower() == 'yes':
                guest_rot = 'yes'
            elif lines[i][1].lower() == 'no':
                guest_rot = 'no'
            else:
                print('Wrong input! Please choose yes or no to allow guest rotation')
                sys.exit(1)
        elif lines[i][0] == 'cation':
            cation = lines[i][1]
        elif lines[i][0] == 'anion':
            anion = lines[i][1]
        elif lines[i][0] == 'num_cations':
            num_cations = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'num_cat_ligbox':
            num_cat_ligbox = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'buffer_x':
            buffer_x = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'buffer_y':
            buffer_y = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_buffer':
            lig_buffer = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_dihcf_force':
            rec_dihcf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_discf_force':
            rec_discf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_distance_force':
            lig_distance_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_angle_force':
            lig_angle_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_dihcf_force':
            lig_dihcf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_com_force':
            rec_com_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_com_force':
            lig_com_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_range':
            l1_range = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'min_adis':
            min_adis = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'max_adis':
            max_adis = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'host_rest_type':
            if lines[i][1].lower() == 'dihedrals':
                host_rest_type = 'dihedrals'
            elif lines[i][1].lower() == 'distances':
                host_rest_type = 'distances'
            else:
                print('Wrong input! Please use distances or dihedrals to define host restraints')
                sys.exit(1)
        elif lines[i][0] == 'release_eq':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                release_eq.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'attach_rest':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                attach_rest.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'lambdas':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                lambdas.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'weights':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                weights.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'components':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                components.append(strip_line[j])
        elif lines[i][0] == 'sdr_dist':
            sdr_dist = scripts.check_input('float', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'ntpr':
            ntpr = lines[i][1]
        elif lines[i][0] == 'ntwr':
            ntwr = lines[i][1]
        elif lines[i][0] == 'ntwe':
            ntwe = lines[i][1]
        elif lines[i][0] == 'ntwx':
            ntwx = lines[i][1]
        elif lines[i][0] == 'cut':
            cut = lines[i][1]
        elif lines[i][0] == 'gamma_ln':
            gamma_ln = lines[i][1]
        elif lines[i][0] == 'barostat':
            barostat = lines[i][1]
        elif lines[i][0] == 'amber_ff':
            amber_ff = lines[i][1]
        elif lines[i][0] == 'dt':
            dt = lines[i][1]

# Number of simulations, 1 equilibrium and 1 production
num_sim = 2
rng = 0

# Define free energy components
if fe_type == 'rest':
  components = ['c', 'a', 'l', 't', 'r'] 
elif fe_type == 'sdr':
  components = ['e', 'v'] 
elif fe_type == 'all':
  components = ['a', 'c', 'r', 'l', 't', 'e', 'v'] 

# Create restraint definitions
rest = [rec_dihcf_force, rec_discf_force, lig_distance_force, lig_angle_force, lig_dihcf_force, rec_com_force, lig_com_force]

# Create ion definitions
ion_def = [cation, anion, num_cations]
ion_lig = [cation, anion, num_cat_ligbox]

# Define number of steps for all stages
dic_steps1 = {}
dic_steps2 = {}
dic_steps1['a'] = a_steps1
dic_steps2['a'] = a_steps2
dic_steps1['l'] = l_steps1
dic_steps2['l'] = l_steps2
dic_steps1['t'] = t_steps1
dic_steps2['t'] = t_steps2
dic_steps1['c'] = c_steps1
dic_steps2['c'] = c_steps2
dic_steps1['r'] = r_steps1
dic_steps2['r'] = r_steps2
dic_steps1['v'] = v_steps1
dic_steps2['v'] = v_steps2
dic_steps1['e'] = e_steps1
dic_steps2['e'] = e_steps2

# Organize guests
for i in range(0, len(guest_list)):
  guest_def.append(guest_list[i])

for i in range(0, len(guest_list_code)):
  code_def.append(guest_list_code[i])

if stage == 'equil':
  comp = 'q'
  win = 0

  # Read charges for charge parameterization (if needed)
  if guest_list_charge:
    for i in range(0, len(guest_list_charge)):
      charge_def.append(guest_list_charge[i])

  # Create equilibrium systems for all guests listed in the input file
  for i in range(0, len(guest_def)):
    rng = len(release_eq) - 1
    guest = guest_def[i]
    mol = code_def[i]
    if guest_list_charge:
      guest_charge = charge_def[i]
    if not os.path.exists('./structures/'+host+'-'+guest+'.pdb'):
      continue
    print('Setting up '+str(guest_def[i]))
    # Get number of simulations
    num_sim = len(release_eq)
    # Create aligned initial complex
    anch = build.build_equil(guest, host, mol, H1, H2, H3, min_adis, max_adis, l1_range, amber_ff, final_host_num, guest_charge, sdr_dist)
    if anch == 'anch1':
      aa1_guests.append(guest)
      os.chdir('../')
      continue
    if anch == 'anch2':
      aa2_guests.append(guest)
      os.chdir('../')
      continue
    # Solvate system with ions
    print('Creating box...')
    build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
    # Apply restraints and prepare simulation files
    print('Equil release weights:')
    for i in range(0, len(release_eq)):
      weight = release_eq[i]
      print('%s' %str(weight))
      setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
      shutil.copy('./'+guest+'/disang.rest', './'+guest+'/disang%02d.rest' %int(i))
    shutil.copy('./'+guest+'/disang%02d.rest' %int(0), './'+guest+'/disang.rest')
    setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, eq_steps1, eq_steps2, rng)
    os.chdir('../')
  if len(aa1_guests) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand first anchor L1 for', aa1_guests)
    print ('The ligand is most likely not in the defined binding site in these systems.')
  if len(aa2_guests) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand L2 or L3 anchors for', aa2_guests)
    print ('Try reducing the min_adis parameter in the input file.')
elif stage == 'fe':
  # Create systems for all guests after equilibration
  if not os.path.exists('fe'):
    os.makedirs('fe')
  os.chdir('fe')
  for i in range(0, len(guest_def)):
    guest = guest_def[i]
    mol = code_def[i]
    fwin = len(release_eq) - 1
    if not os.path.exists('../equil/'+guest):
      continue
    print('Setting up '+str(guest_def[i]))
    # Create and move to guest directory
    if not os.path.exists(guest):
      os.makedirs(guest)
    os.chdir(guest)
    # Generate folder and restraints for all components and windows
    for j in range(0, len(components)):
      comp = components[j]
      # Guest conformational release in a small box
      if (comp == 'c'):
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            anch = build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            if anch == 'anch1':
              aa1_guests.append(guest)
              break
            if anch == 'anch2':
              aa2_guests.append(guest)
              break
            print('Creating box for guest only...')
            build.guest_box(guest, mol, lig_buffer, water_model, neut, ion_lig, comp, amber_ff)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, c_steps1, c_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, c_steps1, c_steps2, rng)
        os.chdir('../')
        if len(aa1_guests) != 0:
          break
        if len(aa2_guests) != 0:
          break
      # Host conformational release in a separate box
      elif (comp == 'r'):
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            anch = build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            if anch == 'anch1':
              aa1_guests.append(guest)
              break
            if anch == 'anch2':
              aa2_guests.append(guest)
              break
            print('Creating box for apo host...')
            build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, r_steps1, r_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, r_steps1, r_steps2, rng)
        os.chdir('../')
        if len(aa1_guests) != 0:
          break
        if len(aa2_guests) != 0:
          break
      # Bound state restraints
      elif (comp == 'a' or comp == 'l' or comp == 't'):
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
          anch = build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
          if anch == 'anch1':
            aa1_guests.append(guest)
            break
          if anch == 'anch2':
            aa2_guests.append(guest)
            break
          setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
          steps1 = dic_steps1[comp]
          steps2 = dic_steps2[comp]
          setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, rng)
        os.chdir('../')
        if len(aa1_guests) != 0:
          break
        if len(aa2_guests) != 0:
          break
      elif (comp == 'e' or comp == 'v'):
        steps1 = dic_steps1[comp]
        steps2 = dic_steps2[comp]
        if not os.path.exists('sdr'):
          os.makedirs('sdr')
        os.chdir('sdr')
        for k in range(0, len(lambdas)):
          weight = lambdas[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
            anch = build.build_dec(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.dec_files(temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, weight, lambdas, ntwx)
          else:
            print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
            anch = build.build_dec(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.dec_files(temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, weight, lambdas, ntwx)
        os.chdir('../')
    os.chdir('../')
  if len(aa1_guests) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand first anchor L1 for', aa1_guests)
    print ('The ligand most likely left the binding site during equilibration of these systems.')
  if len(aa2_guests) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand L2 or L3 anchors for', aa2_guests)
    print ('Try reducing the min_adis parameter in the input file.')
elif stage == 'analysis':
  # Free energies MBAR/TI and analytical calculations
  for i in range(0, len(guest_def)):
    guest = guest_def[i]
    analysis.fe_values(blocks, components, temperature, guest, attach_rest, lambdas, weights, dec_int, rest, guest_rot)
    os.chdir('../../')
