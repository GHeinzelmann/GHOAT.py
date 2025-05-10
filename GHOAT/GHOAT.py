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
import numpy as np

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


# Defaults

a_steps1 = 0
a_steps2 = 0
l_steps1 = 0
l_steps2 = 0
t_steps1 = 0
t_steps2 = 0
c_steps1 = 0
c_steps2 = 0
r_steps1 = 0
r_steps2 = 0
e_steps1 = 0
e_steps2 = 0
v_steps1 = 0
v_steps2 = 0
rec_dihcf_force = 0
rec_discf_force = 0
guest_charge = 0
guest_rot = 'no'
final_host_num = 1
buffer_z = 0
num_waters = 0
ion_conc = 0.0
software = 'amber'
dlambda = 0.001

a_itera1 = 0
a_itera2 = 0
l_itera1 = 0
l_itera2 = 0
t_itera1 = 0
t_itera2 = 0
c_itera1 = 0
c_itera2 = 0
r_itera1 = 0
r_itera2 = 0
e_itera1 = 0
e_itera2 = 0
v_itera1 = 0
v_itera2 = 0

ntpr = '1000'
ntwr = '10000'
ntwe = '0'
ntwx = '2500'
cut = '9.0'
barostat = '2'
ti_points = 0

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
        # OpenMM only
        elif lines[i][0] == 'a_itera1':
            a_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'a_itera2':
            a_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'l_itera1':
            l_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'l_itera2':
            l_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 't_itera1':
            t_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 't_itera2':
            t_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'c_itera1':
            c_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'c_itera2':
            c_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'r_itera1':
            r_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'r_itera2':
            r_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'e_itera1':
            e_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'e_itera2':
            e_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'v_itera1':
            v_itera1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'v_itera2':
            v_itera2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'itera_steps':
            itera_steps = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'itcheck':
            itcheck = lines[i][1]
        ####
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
        elif lines[i][0] == 'ti_points':
            ti_points = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
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
        elif lines[i][0] == 'ion_conc':
            ion_conc = scripts.check_input('float', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'buffer_x':
            buffer_x = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'buffer_y':
            buffer_y = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'buffer_z':
            buffer_z = scripts.check_input('float', lines[i][1], input_file, lines[i][0])
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
            if lines[i][1].lower() == 'gaff':
                amber_ff = 'gaff'
            elif lines[i][1].lower() == 'gaff2':
                amber_ff = 'gaff2'
            else:
                print('Wrong input! Available options for force-field are gaff and gaff2')
                sys.exit(1)
        elif lines[i][0] == 'dt':
            dt = lines[i][1]
        elif lines[i][0] == 'software':
            if lines[i][1].lower() == 'openmm':
                software = lines[i][1].lower()
            elif lines[i][1].lower() == 'amber':
                software = lines[i][1].lower()
            else:
                print('Simulation software not recognized, please choose openmm or amber')
                sys.exit(1)


if num_waters == 0 and buffer_z == 0:
  print('Wrong input! Please choose either a number of water molecules or a z buffer value.')
  sys.exit(1)

if num_waters != 0 and buffer_z != 0:
  print('Wrong input! Please choose either a number of water molecules or a z buffer value.')
  sys.exit(1)

# Number of simulations, 1 equilibrium and 1 production
num_sim = 2
rng = 0

# Define free energy components
if fe_type == 'rest':
  components = ['c', 'a', 'l', 't', 'r'] 
elif fe_type == 'sdr':
  components = ['e', 'v'] 
elif fe_type == 'all':
  components = ['a', 'l', 't', 'e', 'v', 'c', 'r'] 

# Create restraint definitions
rest = [rec_dihcf_force, rec_discf_force, lig_distance_force, lig_angle_force, lig_dihcf_force, rec_com_force, lig_com_force]

# Create ion definitions
ion_def = [cation, anion, ion_conc]

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

# Define number of steps for all stages
dic_itera1 = {}
dic_itera2 = {}
dic_itera1['a'] = a_itera1
dic_itera2['a'] = a_itera2
dic_itera1['l'] = l_itera1
dic_itera2['l'] = l_itera2
dic_itera1['t'] = t_itera1
dic_itera2['t'] = t_itera2
dic_itera1['c'] = c_itera1
dic_itera2['c'] = c_itera2
dic_itera1['r'] = r_itera1
dic_itera2['r'] = r_itera2
dic_itera1['v'] = v_itera1
dic_itera2['v'] = v_itera2
dic_itera1['e'] = e_itera1
dic_itera2['e'] = e_itera2

# Organize guests
for i in range(0, len(guest_list)):
  guest_def.append(guest_list[i])

for i in range(0, len(guest_list_code)):
  code_def.append(guest_list_code[i])


# Obtain Gaussian Quadrature lambdas and weights

if dec_int == 'ti':
  if ti_points != 0:
    lambdas = []
    weights = []
    x,y = np.polynomial.legendre.leggauss(ti_points)
    # Adjust Gaussian lambdas
    for i in range(0, len(x)):
      lambdas.append(float((x[i]+1)/2))
    # Adjust Gaussian weights
    for i in range(0, len(y)):
      weights.append(float(y[i]/2))
  else:
    print('Wrong input! Please choose a positive integer for the ti_points variable when using the TI-GQ method')
    sys.exit(1)
  print('lambda values:', lambdas)
  print('Gaussian weights:', weights)
elif dec_int == 'mbar':
  if lambdas == []:
    print('Wrong input! Please choose a set of lambda values when using the MBAR method')
    sys.exit(1)
  if ti_points != 0:
    print('Wrong input! Do not define the ti_points variable when applying the MBAR method, instead choose a set of lambda values')
    sys.exit(1)
  print('lambda values:', lambdas)


# Adjust components and windows for OpenMM

if software == 'openmm' and stage == 'fe':
  components_inp = list(components)
  print(components_inp)
  components = ['t', 'c', 'v']
  attach_rest_inp = list(attach_rest)
  print(attach_rest_inp)
  attach_rest = [ 100.0 ]
  lambdas_inp = list(lambdas)
  print(lambdas_inp)
  lambdas = [ 0.0 ]
  dt = str(float(dt)*1000)
  print(dt)
  cut = str(float(cut)/10)
  print(cut)

  # Convert equil output file
  os.chdir('equil')
  for i in range(0, len(guest_def)):
    guest = guest_def[i]
    rng = len(release_eq) - 1
    if os.path.exists(guest):
      os.chdir(guest)
      convert_file = open('convert.in', 'w')
      convert_file.write('parm full.prmtop\n')
      convert_file.write('trajin md%02d.dcd\n' %rng)
      convert_file.write('trajout md%02d.rst7 onlyframes 10\n' %rng)
      convert_file.close()
      sp.call('cpptraj -i convert.in > convert.log', shell=True)
      os.chdir('../')
  os.chdir('../')

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
    build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, buffer_z, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
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
            build.guest_box(guest, mol, lig_buffer, water_model, neut, ion_def, comp, amber_ff)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, c_steps1, c_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, c_steps1, c_steps2, rng)
        os.chdir('../')
        if anch != 'all':
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
            build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, buffer_z, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, r_steps1, r_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_rest(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.sim_files(hmr, temperature, mol, num_sim, guest, comp, win, stage, r_steps1, r_steps2, rng)
        os.chdir('../')
        if anch != 'all':
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
        if anch != 'all':
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
            if anch == 'anch1':
              aa1_guests.append(guest)
              break
            if anch == 'anch2':
              aa2_guests.append(guest)
              break
            build.create_box(comp, hmr, guest, host, mol, hmol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, buffer_z, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, final_host_num)
            setup.restraints(guest, host, host_rest_type, final_host_num, H1, H2, H3, rest, weight, stage, mol, comp, sdr_dist, guest_rot)
            setup.dec_files(temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, weight, lambdas, ntwx)
          else:
            print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
            build.build_dec(fwin, min_adis, max_adis, l1_range, H1, H2, H3, hmr, hmol, mol, host, guest, final_host_num, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, amber_ff, dt, sdr_dist)
            setup.dec_files(temperature, mol, num_sim, guest, comp, win, stage, steps1, steps2, weight, lambdas, ntwx)
        os.chdir('../')
        if anch != 'all':
          break
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
  # Free energy analysis for OpenMM
  if software == 'openmm':
    for i in range(0, len(guest_def)):
      guest = guest_def[i]
      analysis.fe_openmm(components, temperature, guest, rest, attach_rest, lambdas, guest_rot, dic_itera1, dic_itera2, itera_steps, dt, dlambda, dec_int, weights, blocks, ti_points)
      os.chdir('../../')
  else:
  # Free energies MBAR/TI and analytical calculations
    for i in range(0, len(guest_def)):
      guest = guest_def[i]
      analysis.fe_values(blocks, components, temperature, guest, attach_rest, lambdas, weights, dec_int, rest, guest_rot, dic_steps1, dic_steps2, dt)
      os.chdir('../../')

# Convert equilibration folders to openmm

if software == 'openmm' and stage == 'equil':

  # Adjust a few variables
  cut = str(float(cut)/10)
  dt = str(float(dt)*1000)

  os.chdir('equil')
  for i in range(0, len(guest_def)):
    guest = guest_def[i]
    rng = len(release_eq) - 1
    if os.path.exists(guest):
      print(guest)
      os.rename(guest, guest+'-amber')
      os.mkdir(guest)
      os.chdir(guest)
#      shutil.copy('../'+guest+'-amber/equil-%s.pdb' % mol.lower(), './')
      shutil.copy('../'+guest+'-amber/cv.in', './')
      shutil.copy('../'+guest+'-amber/assign.dat', './')
      for file in glob.glob('../'+guest+'-amber/vac*'):
        shutil.copy(file, './')
      for file in glob.glob('../'+guest+'-amber/full*'):
        shutil.copy(file, './')
      for file in glob.glob('../'+guest+'-amber/disang*'):
        shutil.copy(file, './')
      for file in glob.glob('../'+guest+'-amber/build*'):
        shutil.copy(file, './')
      for file in glob.glob('../'+guest+'-amber/tleap_solvate*'):
        shutil.copy(file, './')
      fin = open('../../run_files/local-equil-op.bash', "rt")
      data = fin.read()
      data = data.replace('RANGE', '%02d' %rng)
      fin.close()
      fin = open('run-local.bash', "wt")
      fin.write(data)
      fin.close()
      fin = open('../../run_files/PBS-Op', "rt")
      data = fin.read()
      data = data.replace('STAGE', stage).replace('POSE', guest)
      fin.close()
      fin = open('PBS-run', "wt")
      fin.write(data)
      fin.close()
      fin = open('../../run_files/SLURMM-Op', "rt")
      data = fin.read()
      data = data.replace('STAGE', stage).replace('POSE', guest)
      fin.close()
      fin = open('SLURMM-run', "wt")
      fin.write(data)
      fin.close()
      for j in range(0, len(release_eq)):
        fin = open('../../lib/equil.py', "rt")
        data = fin.read()
        data = data.replace('LIG', mol.upper()).replace('TMPRT', str(temperature)).replace('TSTP', str(dt)).replace('GAMMA_LN', str(gamma_ln)).replace('STG','%02d' %j).replace('CTF', cut)
        if hmr == 'yes':
          data = data.replace('PRMFL', 'full.hmr.prmtop')
        else:
          data = data.replace('PRMFL', 'full.prmtop')
        if j == rng:
          data = data.replace('TOTST', str(eq_steps2))
        else:
          data = data.replace('TOTST', str(eq_steps1))
        fin.close()
        fin = open('equil-%02d.py' %j, "wt")
        fin.write(data)
        fin.close()
      os.chdir('../')
      shutil.rmtree('./'+guest+'-amber')
  print(os.getcwd())


if software == 'openmm' and stage == 'fe':

  # Redefine input arrays

  components = list(components_inp)
  attach_rest = list(attach_rest_inp)
  lambdas = list(lambdas_inp)
  lambdas_rest = []
  for i in attach_rest:
    lbd_rst=float(i)/float(100)
    lambdas_rest.append(lbd_rst)
  Input = lambdas_rest
  lambdas_rest = ['{:.5f}'.format(elem) for elem in Input]

  # Start script

  print('')
  print('#############################')
  print('## OpenMM patch for BAT.py ##')
  print('#############################')
  print('')
  print('Components: ', components)
  print('')
  print('Decoupling lambdas: ', lambdas)
  print('')
  print('Restraint lambdas: ', lambdas_rest)
  print('')

  # Generate folder and restraints for all components and windows
  for i in range(0, len(guest_def)):
    if not os.path.exists(guest_def[i]):
      continue
    os.chdir(guest_def[i])
    for j in range(0, len(components)):
      comp = components[j]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'r' or comp == 'c':
          if not os.path.exists('rest'):
            os.makedirs('rest')
          os.chdir('rest')
          if not os.path.exists(comp+'-comp'):
            os.makedirs(comp+'-comp')
          os.chdir(comp+'-comp')
          itera1 = dic_itera1[comp]
          itera2 = dic_itera2[comp]
          shutil.copy('../../../../run_files/local-rest-op.bash', './run-local.bash')
          fin = open('../../../../run_files/PBS-Op', "rt")
          data = fin.read()
          data = data.replace('POSE', comp).replace('STAGE', guest_def[i])
          fin.close()
          fin = open('PBS-run', "wt")
          fin.write(data)
          fin.close()
          fin = open('../../../../run_files/SLURMM-Op', "rt")
          data = fin.read()
          data = data.replace('POSE', comp).replace('STAGE', guest_def[i])
          fin.close()
          fin = open('SLURMM-run', "wt")
          fin.write(data)
          fin.close()
          fin = open('../../../../lib/rest.py', "rt")
          data = fin.read()
          data = data.replace('LAMBDAS', '[%s]' % ' , '.join(map(str, lambdas_rest))).replace('LIG', mol.upper()).replace('TMPRT', str(temperature)).replace('TSTP', str(dt)).replace('SPITR', str(itera_steps)).replace('PRIT', str(itera2)).replace('EQIT', str(itera1)).replace('ITCH', str(itcheck)).replace('GAMMA_LN', str(gamma_ln)).replace('CMPN', str(comp)).replace('CTF', cut).replace('BLCKS', str(blocks))
          if hmr == 'yes':
            data = data.replace('PRMFL', 'full.hmr.prmtop')
          else:
            data = data.replace('PRMFL', 'full.prmtop')
          fin.close()
          fin = open('rest.py', "wt")
          fin.write(data)
          fin.close()
          if comp == 'c':
            shutil.copy('../../../../'+stage+'/'+guest_def[i]+'/rest/c00/disang.rest', './')
            for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/rest/c00/full*'):
              shutil.copy(file, './')
          else:
            shutil.copy('../../../../'+stage+'/'+guest_def[i]+'/rest/t00/disang.rest', './')
            shutil.copy('../../../../'+stage+'/'+guest_def[i]+'/rest/t00/cv.in', './')
            for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/rest/t00/full*'):
              shutil.copy(file, './')
      elif comp == 'e' or comp == 'v':
        if not os.path.exists('sdr'):
          os.makedirs('sdr')
        os.chdir('sdr')
        if dec_int == 'mbar':
          if not os.path.exists(comp+'-comp'):
            os.makedirs(comp+'-comp')
          os.chdir(comp+'-comp')
          itera1 = dic_itera1[comp]
          itera2 = dic_itera2[comp]
          shutil.copy('../../../../run_files/local-sdr-op.bash', './run-local.bash')
          fin = open('../../../../run_files/PBS-Op', "rt")
          data = fin.read()
          data = data.replace('POSE', comp).replace('STAGE', guest_def[i])
          fin.close()
          fin = open('PBS-run', "wt")
          fin.write(data)
          fin.close()
          fin = open('../../../../run_files/SLURMM-Op', "rt")
          data = fin.read()
          data = data.replace('POSE', comp).replace('STAGE', guest_def[i])
          fin.close()
          fin = open('SLURMM-run', "wt")
          fin.write(data)
          fin.close()
          fin = open('../../../../lib/sdr.py', "rt")
          data = fin.read()
          data = data.replace('LAMBDAS', '[%s]' % ' , '.join(map(str, lambdas))).replace('LIG', mol.upper()).replace('TMPRT', str(temperature)).replace('TSTP', str(dt)).replace('SPITR', str(itera_steps)).replace('PRIT', str(itera2)).replace('EQIT', str(itera1)).replace('ITCH', str(itcheck)).replace('GAMMA_LN', str(gamma_ln)).replace('CMPN', str(comp)).replace('CTF', cut).replace('BLCKS', str(blocks))
          if hmr == 'yes':
            data = data.replace('PRMFL', 'full.hmr.prmtop')
          else:
            data = data.replace('PRMFL', 'full.prmtop')
          fin.close()
          fin = open('sdr.py', "wt")
          fin.write(data)
          fin.close()
          shutil.copy('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/disang.rest', './')
          shutil.copy('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/cv.in', './')
          for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/full*'):
            shutil.copy(file, './')
          for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/vac*'):
            shutil.copy(file, './')
          for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/tleap_solvate*'):
            shutil.copy(file, './')
          for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/build*'):
            shutil.copy(file, './')
          for file in glob.glob('../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/full*'):
            shutil.copy(file, './')
        elif dec_int == 'ti':
          if not os.path.exists(comp+'-comp'):
            os.makedirs(comp+'-comp')
          os.chdir(comp+'-comp')
          itera1 = int(dic_itera1[comp]*itera_steps)
          itera2 = int(dic_itera2[comp]/2)
          for k in range(0, len(lambdas)):
            if not os.path.exists('%s%02d' %(comp, int(k))):
              os.makedirs('%s%02d' %(comp, int(k)))
            os.chdir('%s%02d' %(comp, int(k)))
            shutil.copy('../../../../../run_files/local-sdr-op-ti.bash', './run-local.bash')
            fin = open('../../../../../run_files/SLURMM-Op', "rt")
            data = fin.read()
            data = data.replace('STAGE', guest_def[i]).replace('POSE', '%s%02d' %(comp, int(k)))
            fin.close()
            fin = open("SLURMM-run", "wt")
            fin.write(data)
            fin.close()
            fin = open('../../../../../run_files/PBS-Op', "rt")
            data = fin.read()
            data = data.replace('STAGE', guest_def[i]).replace('POSE', '%s%02d' %(comp, int(k)))
            fin.close()
            fin = open("PBS-run", "wt")
            fin.write(data)
            fin.close()
            fin = open('../../../../../lib/equil-sdr.py', "rt")
            data = fin.read()
            data = data.replace('LBD0', '%8.6f' % lambdas[k]).replace('LIG', mol.upper()).replace('TMPRT', str(temperature)).replace('TSTP', str(dt)).replace('SPITR', str(itera_steps)).replace('PRIT', str(itera2)).replace('EQIT', str(itera1)).replace('ITCH', str(itcheck)).replace('GAMMA_LN', str(gamma_ln)).replace('CMPN', str(comp)).replace('CTF', cut)
            if hmr == 'yes':
              data = data.replace('PRMFL', 'full.hmr.prmtop')
            else:
              data = data.replace('PRMFL', 'full.prmtop')
            fin.close()
            fin = open('equil-sdr.py', "wt")
            fin.write(data)
            fin.close()
            fin = open('../../../../../lib/sdr-ti.py', "rt")
            data = fin.read()
            # "Split" initial lambda into two close windows 
            lambda1 = float(lambdas[k] - dlambda/2)
            lambda2 = float(lambdas[k] + dlambda/2)
            data = data.replace('LBD1', '%8.6f' % lambda1).replace('LBD2', '%8.6f' % lambda2).replace('LIG', mol.upper()).replace('TMPRT', str(temperature)).replace('TSTP', str(dt)).replace('SPITR', str(itera_steps)).replace('PRIT', str(itera2)).replace('EQIT', str(itera1)).replace('ITCH', str(itcheck)).replace('GAMMA_LN', str(gamma_ln)).replace('CMPN', str(comp)).replace('CTF', cut).replace('BLCKS', str(blocks))
            if hmr == 'yes':
              data = data.replace('PRMFL', 'full.hmr.prmtop')
            else:
              data = data.replace('PRMFL', 'full.prmtop')
            fin.close()
            fin = open('sdr-ti.py', "wt")
            fin.write(data)
            fin.close()
            shutil.copy('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/disang.rest', './')
            shutil.copy('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/cv.in', './')
            for file in glob.glob('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/full*'):
              shutil.copy(file, './')
            for file in glob.glob('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/vac*'):
              shutil.copy(file, './')
            for file in glob.glob('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/tleap_solvate*'):
              shutil.copy(file, './')
            for file in glob.glob('../../../../../'+stage+'/'+guest_def[i]+'/sdr/v00/build*'):
              shutil.copy(file, './')
            os.chdir('../')
      os.chdir('../../')
    # Clean up amber windows
    dirpath = os.path.join('rest', 't00')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
      shutil.rmtree(dirpath)
    dirpath = os.path.join('rest', 'amber_files')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
      shutil.rmtree(dirpath)
    dirpath = os.path.join('rest', 'c00')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
      shutil.rmtree(dirpath)
    dirpath = os.path.join('sdr', 'v00')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
      shutil.rmtree(dirpath)
    dirpath = os.path.join('sdr', 'amber_files')
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
      shutil.rmtree(dirpath)
    os.chdir('../')



