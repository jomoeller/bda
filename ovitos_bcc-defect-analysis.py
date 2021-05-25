##########################################################################################
#                                ______  ______    ___                                   #
#                                | ___ \ |  _  \  / _ \                                  #
#                                | |_/ / | | | | / /_\ \                                 #
#                                | ___ \ | | | | |  _  |                                 #
#                                | |_/ / | |/ /  | | | |                                 #
#                                \____/  |___/   \_| |_/                                 #
#                                                                                        #
#                                  BCC Defect Analysis                                   #
#      - A novel method for identifying defects in body-centered cubic crystals -        #
#                                                                                        #
# Developed and written by:                                                              #
#     Johannes J. Möller (johannes.moeller@fau.de),                                      #
#     Department of Materials Science and Engineering, Institute I,                      #
#     Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU), Germany.                  #
#                                                                                        #
# If you used the BDA method to analyze your simulation results,                         #
# please cite the BDA in your publications as follows:                                   #
#     J.J. Möller and E. Bitzek                                                          #
#     BDA: A novel method for identifying defects in body-centered cubic crystals        #
#     MethodsX 3 (2016), 279-288                                                         #
#     http://dx.doi.org/10.1016/j.mex.2016.03.013                                        #
# If possible, please also include a link to the website:                                #
#     http://jomoeller.github.io/bda/                                                    #
##########################################################################################

import os, sys, subprocess, argparse, platform
import time
from numpy import *
import ovito
from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
# import threading
import multiprocessing
from multiprocessing import Process, Value, Array, Pool
#from array import array

##########################################################################################
# TESTS FOR SURFACE ATOMS
##########################################################################################

def is_surface(i):
  global atom_coord,atom_csp,atom_cna
  coord=atom_coord[i]
  csp=atom_csp[i]
  cna=atom_cna[i]
  if cna != 3 and (coord <= 11 or is_neighbor2surface(i)):
    atom_defect[i]=srf
    return True
  else: 
    return False

def is_neighbor2surface(i):
  if atom_cna[i] != 3 and atom_coord[i] < 14:
    count=0
    neighbors=get_neighbors(i)
    for n in neighbors:
      if atom_cna[n] != 3 and (atom_defect[n]==srf or atom_coord[n]<=11):
        count+=1
        atom_defect[n]=srf # might be redundant
    if count>=4:
      atom_defect[i]=srf 
      return True
    else:
      return False
  else:
    return False  

##########################################################################################
# TEST FOR NON-SCREW DISLOCATION
##########################################################################################

def is_dislo(i):
  global atom_coord
  coord=atom_coord[i]
  if coord >=12 and coord != 14: # 12-, 13-, and 15-coordinated atoms
    nr_14 = 0
    nr_non14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect
    for n in neighbors:
      if atom_coord[n]!=14:
        nr_non14+=1
      if atom_coord[n]==14:
        nr_14+=1
    if nr_non14 > nr_14:
      atom_defect[i]=dsl
      return True
    else:
      return False
  elif coord == 14: 
    nr_14 = 0
    nr_non14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect    
    for n in neighbors:
      if atom_coord[n] >= 12 and atom_coord[n]!=14:
        nr_non14+=1
      if atom_coord[n]==14:
        nr_14+=1
    if nr_non14 >=4 and nr_14 <= 6 and nr_perfect <= 4:
      atom_defect[i]=dsl
      return True
    else:
      return False    
  else:
    return False

##########################################################################################
# TEST FOR VACANCY
##########################################################################################

def is_vac(i):
  global atom_coord,atom_csp,atom_cna
  coord=atom_coord[i]
  csp=atom_csp[i]
  ############################
  # Mono-vacancy
  ############################
  if coord == 13 and csp < 1:
    nr_perfect = 0
    nr_12_4 = 0 # for vacancy row
    nr_13 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    # 
    for n in neighbors:
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] > 4: 
        nr_13+=1
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] > 4:  
        nr_12_4+=1
    if (nr_13==4 and nr_perfect == 9) or (nr_13==2 and nr_12_4 == 2 and nr_perfect == 9):
      atom_defect[i]=vcn 
      return True
  elif coord == 13 and csp > 4:
    nr_12_1 = 0 # for vacancy row
    nr_12_4 = 0 # for vacancy row
    nr_13_1 = 0
    nr_13_4 = 0
    nr_13 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] < 1:  
        nr_12_1+=1
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] > 4:  
        nr_12_4+=1
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] < 1:  
        nr_13_1+=1
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] > 4:  
        nr_13_4+=1
      if atom_cna[n] != 3 and atom_coord[n]==13:  
        nr_13+=1
    if (nr_13_1 == 3 and nr_13_4 == 3 and nr_perfect == 7) or (nr_12_1 == 2 and nr_12_4 == 2 and nr_13_1 == 1 and nr_13_4 == 1 and nr_perfect == 7) or (nr_13 == 6 and nr_perfect == 7) or (nr_13_4 == 4 and nr_perfect > 7):
      atom_defect[i]=vcn 
      return True
    else: 
      return False   
  ############################
  # Di-vacancy (= vacancy row)
  ############################
  elif coord == 12 and csp > 4:       
    nr_12_1 = 0
    nr_12_4 = 0
    nr_13_1 = 0
    nr_13_4 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] < 1:  
        nr_12_1+=1     
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] > 4:  
        nr_12_4+=1     
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] < 1:  
        nr_13_1+=1
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] > 4:  
        nr_13_4+=1
    if nr_12_1 == 2 and nr_12_4 == 1 and nr_13_1 == 2 and nr_13_4 == 4 and nr_perfect == 3:
      atom_defect[i]=vcn 
      return True
    else: 
      return False   
  elif coord == 12 and csp < 1:  
    nr_12_4 = 0     
    nr_13_1 = 0
    nr_13_4 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n] != 3 and atom_coord[n]==12 and atom_csp[n] > 4:  
        nr_12_4+=1
      if atom_cna[n] != 3 and atom_coord[n]==13 and atom_csp[n] > 4:  
        nr_13_4+=1
    if nr_12_4 == 2 and nr_13_4 == 4 and nr_perfect == 6:
      atom_defect[i]=vcn 
      return True
    else: 
      return False   
    
  
  else: 
    return False

##########################################################################################
# TEST FOR TWIN BOUNDARY AND SCREW DISLOCATION
##########################################################################################

def is_twin(i):
  global atom_coord,atom_cna,atom_csp
  coord=atom_coord[i]
  cna=atom_cna[i]
  csp=atom_csp[i]
  if coord == 13 and csp > 4.5:
    nr_13 = 0
    nr_14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n]!=3 and atom_coord[n]==13: 
        nr_13+=1
      if atom_cna[n]!=3 and atom_coord[n]==14: 
        nr_14+=1
    if nr_13 == 5 and nr_14 == 2 and nr_perfect == 6: 
      atom_defect[i]=twn
      return True
  if coord == 14 and csp > 8:
    nr_14 = 0
    nr_non14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    if nr_perfect <= 8: 
      atom_defect[i]=twn
      return True
  if coord == 14:
    nr_13 = 0
    nr_14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n]!=3 and atom_coord[n]==13: # ideally: 0
        nr_13+=1
      if atom_cna[n]!=3 and atom_coord[n]==14: # ideally: 6
        nr_14+=1
    if nr_perfect >= 6 and nr_perfect <= 9 and (nr_14 >= 4 or (nr_13 == 4 and nr_14 == 2)): 
      atom_defect[i]=twn
      return True
    else: 
      return False
  elif coord == 13 and csp < 1:
    nr_14 = 0
    nr_non14 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_cna[n]!=3 and atom_coord[n]==14 and atom_csp[n] > 8: # ideally: 6
        nr_14+=1
    if nr_14 == 4: 
      atom_defect[i]=twn
      return True
    else: 
      return False
  else: 
    return False

##########################################################################################
# TEST FOR PLANAR FAULT
##########################################################################################

def is_planarfault(i):
  global atom_cna,atom_coord
  coord=atom_coord[i]
  cna=atom_cna[i]

# OLD DEFINITION:
#     elif nr_13 >=6 and nr_13 <= 7 and nr_perfect >= 6 and nr_perfect <= 7:
#       # a multi-layer planar fault at the surface:
#       atom_defect[i]=plf 
#       return True

  if cna != 3 and coord == 12:
    nr_12 = 0
    nr_13 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_coord[n]==12:
        nr_12+=1
      if atom_coord[n]==13:
        nr_13+=1
    if nr_12>=9 and nr_perfect == 0:
      # the interior of a multi-layer planar fault:
      atom_defect[i]=plf 
      return True
    elif (nr_12 >= 3 and nr_12 <= 6) and (nr_13 >= 7 and nr_13 <= 9)  and nr_perfect == 0:
      # the interior of a multi-layer planar fault at the partial tip
      atom_defect[i]=plf 
      return True
    elif nr_12 >= 6 and nr_13 >= 3 and nr_perfect == 0:
      # the interior of a multi-layer planar fault but adjacent to top layer:
      atom_defect[i]=plf 
      return True
    else:
      return False
  elif coord == 13:
    nr_12 = 0
    nr_13 = 0
    neighbors=get_neighbors(i)
    nr_nonperfect = len(neighbors)
    nr_perfect = coord - nr_nonperfect  
    for n in neighbors:
      if atom_coord[n]==12:
        nr_12+=1
      if atom_coord[n]==13:
        nr_13+=1
    if nr_12 + nr_13 == 9 and nr_13 >= 7:
      # a two-layer planar fault:
      atom_defect[i]=plf 
      return True
    elif nr_13 == 6 and nr_12 == 3 and nr_perfect == 4:
      # the top layer of a multi-layer planar fault
      atom_defect[i]=plf 
      return True
    elif nr_13 == 6 and nr_12 <= 1 and nr_perfect >= 6:
      # the top layer of a multi-layer planar fault: partial and edge of partials
      atom_defect[i]=plf 
      return True
    elif nr_13 >= 7 and nr_12 <= 4 and nr_perfect <= 3:
      # the top layer of a multi-layer planar fault: partial and edge of partials
      atom_defect[i]=plf 
      return True
      #     elif nr_13 == 6 and nr_perfect == 7:
#       # the top layer of a multi-layer planar fault:
#       atom_defect[i]=plf 
#       return True
    else:
      return False
  else: 
    return False

##########################################################################################
# TEST FOR MOST COMMON NEIGHBOR DEFECT
##########################################################################################

def common_neighbor_defect(i):
  global atom_defect
  defect_count=[0,0,0,0,0,0,0]
  neighbors=get_neighbors(i)
  for n in neighbors:
    for d in range(1,6):
      if atom_defect[n] == d:
        defect_count[d]+=1
  max_count=0
  most_common=[]
  for d in range(1,6):
    if defect_count[d] > max_count and defect_count[d] >= 3:
      most_common=[d]
      max_count=defect_count[d]
    elif defect_count[d] == max_count and defect_count[d] >= 3:
      most_common.append(d)
#  print(i,most_common,max_count)
  if len(most_common)==1:
    return most_common[0]
  else:
#    return atom_defect[i]
    return els

##########################################################################################
# OUTPUT ATOMS
##########################################################################################

def write_atom(i):
  global atom_nrs,atom_types,atom_masses,atom_pos,atom_cna,atom_coord,atom_csp,atom_defect

  # write line to outfile:
  outstr="%10d %3d %12.10f %12.6f %12.6f %12.6f %2d %2d %10.6f %d\n" % \
         (atom_nrs[i],atom_types[i],atom_masses[i],\
          atom_pos[i][0],atom_pos[i][1],atom_pos[i][2],\
          atom_cna[i],atom_coord[i],atom_csp[i],atom_defect[i])
  f.write(outstr)

##########################################################################################
# GET AN ATOMS NEIGHBORS
##########################################################################################

def get_neighbors(i):
  global neighbor_finder
  global atom_coord
  return [neigh.index for neigh in neighbor_finder.find(i) if atom_cna[neigh.index] != 3 or atom_coord[neigh.index] != 14]
  
##########################################################################################
# Function to control option parsing in Python
##########################################################################################

def controller():
  global VERBOSE,bc,br,alats,filenames,include_perfect,keep_unidentified


  p = argparse.ArgumentParser(description='BDA (BCC Defect Analysis) - A novel method for identifying defects in body-centered cubic crystals. Developed and written by Johannes J. Moeller, johannes.moeller@fau.de. Please visit http://jomoeller.github.io/bda/ for further information.',
                                     prog='ovitos_bcc-defect-analysis_v2.py',
                                    usage= '%(prog)s [options]')
  p.add_argument('-c','--config',nargs='+',help='Atomistic configuration(s) in IMD format',required=True)
  p.add_argument('-b','--boundary-conditions',nargs=3,help='Boundary conditions (0:free|1:periodic)',type=int,default=[0,0,0],metavar=('X','Y','Z'))
  p.add_argument('-a','--lattice-parameter',nargs=1,help='BCC Lattice parameter',type=float)
  p.add_argument('-p','--potential',nargs=1,help='Potential',type=str)
  p.add_argument('-r','--boundary-region',nargs=1,help='Regions to cut away from non-periodic boundaries (default: 5)',type=float,default=[5])
  p.add_argument('-i','--include-perfect',help='Include perfect lattice atoms in exported files',action='store_true')
  p.add_argument('-k','--keep-unidentified',help='Keep unidentified and do no optimization loops',action='store_true')

  args=p.parse_args()
#  print(args.config,args.boundary_conditions,args.lattice_parameter,args.potential)

  if args.config:
    filenames = args.config

  if args.boundary_conditions:
    bc = args.boundary_conditions

  if args.boundary_region:
    br = args.boundary_region

  if args.include_perfect:
    include_perfect = True
  else:
    include_perfect = False

  if args.keep_unidentified:
    keep_unidentified = True
  else:
    keep_unidentified = False

  known_potentials = [['Chiesa','DD_CS3-33','Men-II','Chamati','Gordon','MPG20','Marinica11','Rosato'],[2.8665,2.8665,2.8553,2.8661,2.85516,2.85516,2.814767,2.86650]] 
  
  alats=[]
  if args.lattice_parameter:
    for file in filenames:
      alats.append(args.lattice_parameter[0])
    #print("Lattice parameter for all configurations: ",alats[0])
  elif args.potential:
    for pot in known_potentials[0]:
      if pot in args.potential:
        alats.append(known_potentials[1][known_potentials[0].index(pot)])
#        print("Lattice parameter of ",pot," potential for all configurations: ",alats[0])
  else:
    for pot in known_potentials[0]:
      for file in filenames:
        if pot in file:
          alats.append(known_potentials[1][known_potentials[0].index(pot)])
  #        print("Lattice parameter of ",pot," potential for",file," : ",alats[-1])
    if alats == []:
      errstr=[str(pot) for pot in known_potentials[0]]
      p.error('Either --lattice-parameter or --potential is required or the filename must contain one of the recognizable potential names: '+str(errstr))


##########################################################################################
# MAIN PART
##########################################################################################

def main():
  global atom_nrs,atom_types,atom_masses,atom_pos,atom_coord,atom_csp,atom_cna,atom_defect,atom_neighbors,neighbor_finder
  global lasti,max_neighbors,blk,srf,vcn,dsl,twn,plf,els,include_perfect,keep_unidentified
  global f,filenames,alats,bc,br
  # Handle arguments passed to the script:
  controller()

  # checking for current Ovito version:
  print("This is the BCC Defect Analysis working with OVITO", ovito.version_string)

  # Handle non-periodic boundary conditions:
  if bc[0] == 0: xtrafo=1.1
  else: xtrafo=1
  if bc[1] == 0: ytrafo=1.1
  else: ytrafo=1
  if bc[2] == 0: ztrafo=1.1
  else: ztrafo=1

  # Define numbers for defects:
  blk=0
  srf=1
  vcn=2
  dsl=3
  twn=4
  plf=5
  els=6

  for file in filenames:
    stime=0
    node = None
    print("Working on file: ", file)
#    print(include_perfect)
    # Get the corresponding lattice parameter and cutoff radius:
    alat=alats[filenames.index(file)]
    nearest_neighbors=8
    nn_cutoff=(sqrt(3)/2+1)/2*alat
    nn2_cutoff=(sqrt(2)+1)/2*alat
    print("Using lattice parameter %.4f Angstroms (cutoff for coordination analyis: %.4f Angstroms)" % (alat,nn2_cutoff))

    # Import the file to OVITO and immediately remove it from the viewport (to possibly save memory)  
    print("Importing file...", end="",flush=True)
    time1=time.time()
    node = import_file(file)
    data = node.compute()
    box = asarray(data.cell)
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    # Get the min and max values of the imported configuration:
    pos_min=amin(data.particles.positions, axis=0)
    pos_max=amax(data.particles.positions, axis=0)
    dist=[0,0,0]
    slice=[0,0,0]
    for i in range(3):
      dist[i]=(pos_max[i]+pos_min[i])/2
      slice[i]=(pos_max[i]-pos_min[i])-2*br[0]

    # define the modifiers to be applied:
    trafo=AffineTransformationModifier(operate_on = {'cell'},transformation=[[xtrafo,0,0,0],[0,ytrafo,0,0],[0,0,ztrafo,0]])
    trafo2=AffineTransformationModifier(operate_on = {'cell'},transformation=[[1/xtrafo,0,0,0],[0,1/ytrafo,0,0],[0,0,1/ztrafo,0]])
    csp=CentroSymmetryModifier(num_neighbors = nearest_neighbors)
    cna=CommonNeighborAnalysisModifier(mode = CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
    coord=CoordinationNumberModifier(cutoff = nn2_cutoff)

    # append the modifiers to the node:
    node.modifiers.append(trafo)
    node.modifiers.append(cna)
    print("Computing adaptive common neighbor analysis...", end="",flush=True)
    time1=time.time()
    node.compute() 
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    node.modifiers.append(coord)
    print("Computing coordination analysis...", end="",flush=True)
    time1=time.time()
    node.compute() 
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    print("Computing centrosymmetry parameter...", end="",flush=True)
    time1=time.time()
    node.modifiers.append(csp)
    node.compute() 
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    node.modifiers.append(trafo2)

    # cut away the non-periodic boundary regions if desired:
    if br != 0:
      if bc[0] == 0 and br[0] != 0: 
        slice1=SliceModifier(normal=(1,0,0),slice_width=slice[0]) 
        node.modifiers.append(slice1)
        slice1.distance=dist[0]
      if bc[1] == 0 and br[0] != 0:     
        slice2=SliceModifier(normal=(0,1,0),slice_width=slice[1]) 
        node.modifiers.append(slice2)
        slice2.distance=dist[1]
      if bc[2] == 0 and br[0] != 0: 
        slice3=SliceModifier(normal=(0,0,1),slice_width=slice[2]) 
        node.modifiers.append(slice3)
        slice3.distance=dist[2]
    
    print("Cutting away atoms at non-periodic boundaries...", end="",flush=True)
    time1=time.time()
    data = node.compute()
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    # We start here with the output in case also perfect atoms should be included in the output:
    nr_atoms=data.particles.count	# will be overwritten later on
    filename=file + ".bda"
    f = open(filename, 'w')
    f.write('#F A 1 1 1 3 0 4 \n')
    f.write('#C number type mass x y z cna coord csp defect\n')
    f.write('#X           '+str.format("{0:" ">12.6f}",box[0][0])+' '+str.format("{0:" ">12.6f}",box[0][1])+' '+str.format("{0:" ">12.6f}",box[0][2])+'\n')
    f.write('#Y           '+str.format("{0:" ">12.6f}",box[1][0])+' '+str.format("{0:" ">12.6f}",box[1][1])+' '+str.format("{0:" ">12.6f}",box[1][2])+'\n')
    f.write('#Z           '+str.format("{0:" ">12.6f}",box[2][0])+' '+str.format("{0:" ">12.6f}",box[2][1])+' '+str.format("{0:" ">12.6f}",box[2][2])+'\n')
    f.write('##\n')
    f.write('##\n')
    f.write('#E\n')

    if include_perfect:
    
      print("Writing %d atoms in perfect bcc environment..." % nr_atoms, end="", flush=True)
      atom_nrs=data.particles.identifiers
      atom_types=data.particles.particle_types
      atom_masses=data.particles.masses
      atom_pos=data.particles.positions
      atom_cna=Array('i',data.particles.structure_types)
      atom_coord=Array('i',data.particles['Coordination'])
      atom_csp=Array('f',data.particles['Centrosymmetry'])
      atom_defect=Array('i',[-1]*len(atom_nrs))

      for i in range(nr_atoms):
        if atom_cna[i]==3 and atom_coord[i]==14:
          atom_defect[i]=0
          write_atom(i)
        
    select_perfect=SelectExpressionModifier(expression = 'StructureType==3&&Coordination==14')
    delete_selected=DeleteSelectedModifier()
    node.modifiers.append(select_perfect)
    node.modifiers.append(delete_selected)
    print("Deleting atoms in perfect bcc environment...", end="",flush=True)
    time1=time.time()
    data = node.compute()    
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    print("Preparing non-bcc neighbor finder...", end="",flush=True) 
    time1=time.time()
    neighbor_finder = CutoffNeighborFinder(nn2_cutoff, data)
    time2=time.time()
    ntime=(time2-time1)
    stime+=ntime
    print(" done in %.1f seconds!" % ntime)

    # exporting the node to the file:
    print("Exporting values of ACNA, CN, and CSP to file: ", file + ".ccc")
    export_file(node, file + ".ccc", "imd")

    # We continue to work on the remaining atoms:

    # Oh my god, this is sooo advanced:
    max_neighbors=max(data.particles['Coordination'])
    nr_atoms=data.particles.count
    print("Numer of remaining atoms:", nr_atoms)

    print("Identifying defects...") 

    atom_nrs=data.particles.identifiers
    atom_types=data.particles.particle_types
    atom_masses=data.particles.masses
    atom_pos=data.particles.positions
    atom_cna=Array('i',data.particles.structure_types)
    atom_coord=Array('i',data.particles['Coordination'])
    atom_csp=Array('f',data.particles['Centrosymmetry'])
    atom_defect=Array('i',[-1]*len(atom_nrs))
#    atom_neighbors=Array('i',[-1]*len(atom_nrs))
    atom_neighbors=full((len(atom_nrs),max_neighbors), -1,dtype=int)
#    print(atom_neighbors)
    
    identified=[]
    unidentified=[]
    time1=time.time()
    defect_atoms=0
  
    for i in range(nr_atoms):

      # first check if this atom has already been tested:
      if atom_defect[i] == -1:
        # check only those atoms that are non-bcc according to CNA and CN:
        if atom_cna[i]!=3 or atom_coord[i]!=14:
          defect_atoms+=1
          if is_surface(i):
            write_atom(i)
            pass
          elif is_vac(i):
            identified.append(i)
            pass
          elif is_twin(i):
            identified.append(i)
            pass
          elif is_planarfault(i):
            identified.append(i)
            pass
          elif is_dislo(i):
            identified.append(i)
            pass
          else:
            # not yet identified defects
            # can be twins or dislocations
            atom_defect[i]=els
            unidentified.append(i)
        else:
          # write bulk atoms:
          atom_defect[i]=blk
          write_atom(i)
      else:
        write_atom(i)
 
    print("Number of non-surface defect atoms: ", defect_atoms,"(",defect_atoms/nr_atoms*100,"% of all atoms)")
    print("Identified defect atoms after initial run: ", len(identified),"(",len(identified)/defect_atoms*100,"% )")
    print("Unidentified defect atoms after initial run: ", len(unidentified),"(",len(unidentified)/defect_atoms*100,"% )")

    # Check if an atom's defect is the most common one of its neighbors and occurs >= 3 times
    # else throw it into the list of unidentified atoms.
  
    for i in identified:
      cd = common_neighbor_defect(i)
      if atom_defect[i] == cd:
        write_atom(i)
      else:      
        unidentified.append(i)
        # defect type is not changed here, but later when all atoms have been checked. 
        # Otherwise, the changing behavior could be cascade like.	 
  
    print("Unidentified defect atoms after re-checking already identified atoms: ", len(unidentified),"(",len(unidentified)/defect_atoms*100,"% )")

    loop_count=0
    unidentified_orig=unidentified
    llen=0
 
    for i in unidentified_orig:
      atom_defect[i]=els
 
    # Comment the following while loop for debugging purposes:
    if not keep_unidentified:
      while len(unidentified)/defect_atoms > 0.005 and len(unidentified) != llen:
        loop_count+=1
        print("Entering loop nr.",loop_count)
        list=unidentified
        llen=len(unidentified)
        unidentified=[]
        for i in list:
          cd = common_neighbor_defect(i)
          if atom_defect[i] != cd:
            atom_defect[i] = cd 
          else:
            atom_defect[i]=els
            unidentified.append(i)
        print("Unidentified atoms after loop nr.",loop_count,": ",len(unidentified),"(",len(unidentified)/defect_atoms*100,"% )")

    for i in unidentified_orig:
      write_atom(i)
   
    time2=time.time()      
    f.close()
    print("All bulk and (un)identified atoms written into file: ", filename) 
    print("Took %0.1f seconds" % (time2-time1))

#This idiom means the below code only runs when executed from command line
if __name__ == '__main__':
    main()
