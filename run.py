import os
import sys

cord={}

# print(sys.argv)
chainlen=int(sys.argv[1])
n=int(sys.argv[2])
repeat=int(sys.argv[3])
run_name=sys.argv[4]
clustersize=int(sys.argv[5])


with open('%s.data' % run_name, 'r') as datafile:
    lines = datafile.readlines()
    
    for line in lines:
        content = line.split()
        if content and content[-1] == 'xhi':
            size = 2*float(content[1])
        elif len(content) == 9:
            # print(content)
            if int(content[0]) % chainlen == 0 or int(content[0]) % chainlen == 1:
                cord[int(content[0])] = [float(content[i])+size*float(content[i+3]) for i in range(3,6)]



for i in range(repeat):
    os.system('mkdir %d' % i)
    os.system('cp %s.data cluster.py %d' % (run_name,i))
    os.system('cd %d; python3 cluster.py %d %s.data %d' % (i,chainlen,run_name,clustersize))
    
    # for k in range(30):
    cordm=[]
    element=[]
    with open('%d/result.data' % i, 'r') as datafile:
        lines = datafile.readlines()

        for line in lines:
            content = line.split()
            
            # print(content)
            
            cordm.append([float(content[i]) for i in range(1,4)])
            element.append([int(content[i]) for i in range(4,8)])

    # print(cordm)
    # print(element)

    
    # print('%d/loop_prod_%d' % (i,m))
    with open('%d/%s' % (i,run_name), 'w') as lp:
        
        lp.write('units lj\n')
        lp.write('dimension 3\n')
        lp.write('boundary p p p\n')	
        lp.write('atom_style bond\n\n')	

        lp.write('neighbor 0.50 bin\n')                               
        lp.write('neigh_modify every 1 delay 1 check yes\n\n')         

        lp.write('read_data %s.data extra/bond/per/atom 4 extra/special/per/atom 4 extra/bond/types 1 extra/atom/types 2\n\n' % run_name)

        lp.write('mass 3 1.0\n')
        lp.write('mass 4 1.0\n')

        for m in range(len(cordm)):
            c = element[m]
            lp.write('create_atoms 4 single %f %f %f\n' % (cordm[m][0],cordm[m][1],cordm[m][2]))
            for k in range(4):
                lp.write('set atom %d type 3\n' % c[k])
                lp.write('create_bonds single/bond 3 %d %d\n' % (c[k], m + n*chainlen+1))
            # lp.write('set atom %d type 3\n' % c[1])
            # lp.write('create_bonds 3 %d %d\n' % (c[1], m + n*chainlen))
            # lp.write('set atom %d type 3\n' % c[2])
            # lp.write('create_bonds 3 %d %d\n' % (c[2], m + n*chainlen))
            # lp.write('set atom %d type 3\n' % c[3])
            # lp.write('create_bonds 3 %d %d\n' % (c[3], m + n*chainlen))
            
        

        # lp.write('group end type 3 4')

        # c = element[m]
            # for j in range(4):
            #     lp.write('group %s%d id %d\n' % ('a',c[j],c[j]))
            # lp.write('\n')
            # lp.write('group %s%d id %d %d %d %d\n' % ('c',m,c[0],c[1],c[2],c[3]))

        lp.write('write_data %s_before.data nocoeff' % run_name)

        lp.write('\npair_style hybrid lj/cut 2.5 soft 1.0\n')        
        lp.write('pair_coeff *3 *3 lj/cut 1.0 1.0\n\n')
        # lp.write('pair_coeff 3 4 morse 30.0 2.0 1.0\n')
        lp.write('pair_coeff * 4 soft 50.0\n')
        # lp.write('pair_coeff 4 4 lj/cut 1.0 1.0\n')

        

        lp.write('bond_style hybrid fene harmonic\n')	
        lp.write('bond_coeff 1 fene 30.0 1.5 1.0 1.0\n')	
        lp.write('bond_coeff 2 fene 30.0 1.5 1.0 1.0\n')
        
        lp.write('special_bonds fene\n\n')				


        # lp.write('dump            1 all xtc 20 %s.xtc\n' % run_name)
        # lp.write('dump_modify     1 unwrap yes\n\n')


        lp.write('fix 1 all npt temp 1.0 1.0 1.0 iso 0.0 0.0 1.0\n')
        # for m in range(len(cordm)):
        # lp.write('fix bond end bond/create 1 3 4 1.3 1 iparam 2 3 jparam 4 4\n\n')

        # for m, c in enumerate(element):
        # for m in range(len(cordm)):
        #     c = element[m]
        #     for j in range(4):
        #         x = cordm[m][0] - size*round((cordm[m][0]-cord[c[j]][0])/size)
        #         y = cordm[m][1] - size*round((cordm[m][1]-cord[c[j]][1])/size)
        #         z = cordm[m][2] - size*round((cordm[m][2]-cord[c[j]][2])/size)

        #         # if j == 0:
        #         lp.write('fix f%d a%d spring tether 30 %f %f %f 0.5\n' % (c[j], c[j], x, y, z))
                # else:
                    # lp.write('fix f%d a%d spring tether 40 %f %f %f 0.96\n' % (c[j], c[j], x, y, z))


        lp.write('\ntimestep 0.01\n')
        lp.write('thermo_style custom step temp bonds press pxx pyy pzz xlo xhi ylo yhi zlo zhi lx ly lz\n')
        lp.write('thermo 200\n\n')		

        elastic=2.0
        for ki in range(9):
            lp.write('bond_coeff 3 harmonic %f 0.96\n' % elastic)
            lp.write('run 1000\n\n')
            elastic*=2
        
        # lp.write('bond_coeff 3 fene 30.0 1.5 1.0 1.0\n')
        # lp.write('run 1000\n\n')

        lp.write('write_data %s_gel.data nocoeff\n\n' % run_name)
    
    with open('%d/%s.pbs' % (i,run_name), 'w') as pbsfile:
        pbsfile.write('#!/bin/bash\n')
        pbsfile.write('#SBATCH -p common,scavenger\n')
        pbsfile.write('#SBATCH -J %s\n' % run_name)
        pbsfile.write('#SBATCH -n 1\n')
        pbsfile.write('#SBATCH -N 1\n')
        pbsfile.write('#SBATCH --mem 4G\n')
        pbsfile.write('#SBATCH -o %s.log\n' % run_name)
        pbsfile.write('#SBATCH -e %s.err\n' % run_name)
        pbsfile.write('#SBATCH --mail-type=FAIL\n')
        pbsfile.write('#SBATCH --mail-user=dc314@duke.edu\n\n')

        pbsfile.write('module load OpenMPI/2.0.3\n')
        pbsfile.write('module load LAMMPS/30April19\n\n')
        # for m in range(len(cordm)):
        pbsfile.write('mpirun -np 1 lmp_mpi -in %s\n' % run_name)
    os.system('cd %d; sbatch %s.pbs' % (i,run_name))
