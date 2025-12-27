#!/bin/bash
#SBATCH --job-name=Co2FeSn_opt
#SBATCH --error=error-ecutrho8-%j.err
#SBATCH --output=salida-ecutrho8-%j.out
#SBATCH --partition=compute
#SBATCH --time=2-10
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1

module load quantum-espresso/7.1

kpoints_values=(6 8 10 12 14 16 18 20)

for k in "${kpoints_values[@]}"; do

    cat > "Co2FeSn_${k}_kpoints.in" << EOF
&control
    calculation = 'scf'
    prefix = 'Co2FeSn'
    tstress = .true.
    tprnfor = .true.
    outdir = './outdir'
    pseudo_dir = './pseudos'
    nstep = 100
    restart_mode = 'from_scratch'
/
&system
    ibrav = 14, 
    celldm(1) = 7.80362,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.5,
    celldm(5) = 0.5,
    celldm(6) = 0.5,
    nat = 4, ntyp = 3,
    ecutwfc = 100
    ecutrho = 800
    occupations = 'smearing',
    smearing = 'gaussian',
    degauss = 0.02,          
    nspin = 2,
    starting_magnetization(1) = 0.50,
    starting_magnetization(2) = 0.50,
/
&electrons
    electron_maxstep=200
    conv_thr=1.0d-8
    mixing_beta=0.2
/
&ions
/
ATOMIC_SPECIES
Co  58.933   Co.pbe-spn-kjpaw_psl.0.3.1.UPF
Fe  55.845  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
Sn 118.710  Sn.pbe-dn-kjpaw_psl.1.0.0.UPF

ATOMIC_POSITIONS (crystal)
Fe  0.500000000  0.500000000  0.500000000
Sn  0.000000000  0.000000000  0.000000000
Co  0.250000000  0.250000000  0.250000000
Co  0.750000000  0.750000000  0.750000000
K_POINTS automatic
${k} ${k} ${k} 0 0 0
EOF

    mpiexec.hydra -np 20 pw.x -in Co2FeSn_${k}_kpoints.in > Co2FeSn_${k}_kpoints.out

done

