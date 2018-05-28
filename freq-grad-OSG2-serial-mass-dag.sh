#!/bin/sh

dos2unix -q -k $1

#LIST OF ATOMIC MASS

cat >> atomic_mass.txt << EOF
H 1.008
He 4.003
Li 6.941
Be 9.012
B 10.811
C 12.011
N 14.007
O 15.999
F 18.998
Ne 20.180
Na 22.990
Mg 24.305
Al 26.982
Si 28.086
P 30.974
S 32.065
Cl 35.453
Ar 39.948
K 39.098
Ca 40.078
Sc 44.956
Ti 47.867
V 50.942
Cr 51.996
Mn 54.938
Fe 55.845
Co 58.933
Ni 58.693
Cu 63.546
Zn 65.390
Ga 69.723
Ge 72.640
As 74.922
Se 78.960
Br 79.904
Kr 83.800
Rb 85.468
Sr 87.620
Y 88.906
Zr 91.224
Nb 92.906
Mo 95.940
Tc 98.000
Ru 101.070
Rh 102.906
Pd 106.420
Ag 107.868
Cd 112.411
In 114.818
Sn 118.710
Sb 121.760
Te 127.600
I 126.905
Xe 131.293
Cs 132.906
Ba 137.327
La 138.906
Ce 140.116
Pr 140.908
Nd 144.240
Pm 145.000
Sm 150.360
Eu 151.964
Gd 157.250
Tb 158.925
Dy 162.500
Ho 164.930
Er 167.259
Tm 168.934
Yb 173.040
Lu 174.967
Hf 178.490
Ta 180.948
W 183.840
Re 186.207
Os 190.230
Ir 192.217
Pt 195.078
Au 196.967
Hg 200.590
Tl 204.383
Pb 207.200
Bi 208.980
Po 209.000
At 210.000
Rn 222.000
Fr 223.000
Ra 226.000
Ac 227.000
Th 232.038
Pa 231.036
U 238.029
Np 237.000
Pu 244.000
Am 243.000
Cm 247.000
Bk 247.000
Cf 251.000
Es 252.000
Fm 257.000
Md 258.000
No 259.000
Lr 262.000
Rf 261.000
Db 262.000
Sg 266.000
Bh 264.000
Hs 277.000
Mt 268.000
EOF

#DEFINING VARIABLES

FileDir=$(pwd)

CYCLES=2

INPUTFILE_OLD=$1

if [ `grep -i -c "base" $INPUTFILE_OLD` -gt 0 ]
then
  PREFIX_OLD=`grep -v "#" $INPUTFILE_OLD | grep -m 1 -i "base" | sed -e 's/.*base "\(.*\)"/\1/' | awk '{ print $1}'`
else
  PREFIX_OLD=${INPUTFILE_OLD%.???}
fi

XYZFILE=$PREFIX_OLD.xyz
GBWFILE=$PREFIX_OLD.gbw

MOLECULE=$(echo $INPUTFILE_OLD | cut -d_ -f1)

DETAIL=$(echo $INPUTFILE_OLD | cut -d. -f1 | sed 's/.*_\(.*\)/\1/')

CHARGE=`tac $INPUTFILE_OLD | grep -m1 '*xyz*' | awk '{ print $2 }'`
SPIN=`tac $INPUTFILE_OLD | grep -m1 '*xyz*' | awk '{ print $3 }'`

if [ $# -gt 1 ] && [ $(expr "$2" : ".*brokensym.*") -gt 0 ]
then
  BROKENSYM=$(echo $2 | sed 's/_/ /')
  SPIN=$(echo "$(echo $BROKENSYM | awk '{print $2}' | sed 's/\(.\),\(.\)/\1/') - $(echo $BROKENSYM | awk '{print $2}' | sed 's/\(.\),\(.\)/\2/') + 1" | bc)
  BS=bs
  echo "Broken symmetry calculation is requested with $BROKENSYM"
  sleep 3
  SPIN=${SPIN}
fi

PROJECT=${MOLECULE}_${CHARGE}_${SPIN}${BS}_${DETAIL}_numfreq-grad

CURRENT_TIME=$(date +"%Y_%m_%d-%H_%M_%S")
WorkDir=$HOME/outputs/$PROJECT/$CURRENT_TIME

echo "Work Directory is $WorkDir"

mkdir -p $WorkDir

cp $XYZFILE $WorkDir
cp atomic_mass.txt $WorkDir

rm atomic_mass.txt

if [ -f $PREFIX_OLD.gbw ]
then
  cp $GBWFILE $WorkDir
else
  GBWFILE=""
  GUESS="guess pmodel"
fi

cd $WorkDir

#CREATING PYTHON SCRIPT

cat >> eig.py << EOF
#!/usr/bin/env python

import sys
import numpy as np

def diag(variable):
        input = np.loadtxt('I_matrix.txt',dtype='i',delimiter=',')
        w,v = np.linalg.eig(input)
        I = np.prod(w)
        print(I)

data = sys.stdin.read()
diag(data)
EOF

chmod +x "eig.py"

#CREATING EXECUTABLE

cat >> $PROJECT.sh << EOF
#!/bin/sh

tar xzf orca-3.0.3-grad-only.tar.gz
rm orca-3.0.3-grad-only.tar.gz

export RSH_COMMAND=ssh
export PATH=\$PWD/orca_3_0_3_linux_x86-64_grad_only:\$PATH

\$PWD/orca_3_0_3_linux_x86-64_grad_only/orca \$1 > \${1%.???}.log
sleep 3

rm -r orca*
rm *numfreq*tar.gz
rm *gbw1
rm *uco
rm *qro
rm *unso
rm *unoloc
rm *uno

EOF

chmod +x "$PROJECT.sh"

##################################
#JOB1+2: CREATING GRAD INPUTFILES#
##################################

#DEFINING METHODS

FUNCTIONAL_FREQ="b3lyp d3"
BS_S_FREQ="def2-svp zora"
BS_AUX_J_S_FREQ=def2-svp/j
BS_AUX_C_S_FREQ=def2-svp/c
BS_L_FREQ="def2-tzvpp zora"
BS_AUX_J_L_FREQ=def2-tzvpp/j
BS_AUX_C_L_FREQ=def2-tzvpp/c

TM1="Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn"
TM2="Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd"
Ln="La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu"
TM3="Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg"
An="Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr"
Element_list=`echo "$TM1,$TM2,$Ln,$TM3,$An" | sed 's/,/\n/g'`

for elements in $Element_list
do

ELEMENT=`grep -m1 -w "$elements " $XYZFILE | awk '{ print $1 }'`

if [ -z "$ELEMENT" ]
then
  :
else
cat >> basis_freq.txt << EOF
  newgto $ELEMENT "$BS_L_FREQ" end
  newauxgto $ELEMENT "$BS_AUX_J_L_FREQ" end
  newauxgto $ELEMENT "$BS_AUX_C_L_FREQ" end
EOF
fi

done

BASIS_FREQ=`cat basis_freq.txt`

rm basis_freq.txt

###########################################################
#JOB2: CREATING DISTURBED XYZFILES AND INPUTFILES FOR THEM#
###########################################################

if [[ $(sed '1q;d' $XYZFILE | awk '{print $1}') =~ ''[0-9]+$'' ]]; then
  TOTAL_ATOMS=$(echo "$(wc -l $XYZFILE | awk '{ print $1 }') - 2" | bc)
else
  TOTAL_ATOMS=$(wc -l $XYZFILE | awk '{ print $1 }')
fi

SIX_N=$(echo "${TOTAL_ATOMS} * 6" | bc)

tail -n ${TOTAL_ATOMS} $XYZFILE > $XYZFILE.tmp

XYZ=$(cat $XYZFILE.tmp)

#CREATING INPUTS

BASE_DISTURB=${MOLECULE}_${CHARGE}_${SPIN}${BS}_RIJCOSX-B3LYP-D3_BS_${DETAIL}.grad
BASE_OLD=${MOLECULE}_${CHARGE}_${SPIN}${BS}_RIJCOSX-B3LYP-D3_BS_${DETAIL}.sp

NPROCS_DISTURB=1

SIX_N_FROM_0=$(echo "$TOTAL_ATOMS * 6 - 1" | bc)

for j in $(seq 1 1 $CYCLES); do

BASE_FREQ_OLD=$BASE_OLD.$j
INPUTFILE_FREQ_OLD=$BASE_FREQ_OLD.inp
XYZFILE_FREQ_OLD=$BASE_FREQ_OLD.xyz

TARGET=$INPUTFILE_FREQ_OLD

MOINP=""

if [ -f $PREFIX_OLD.gbw ]; then
  MOINP=$(echo -e "!moread\n\n%moinp \"$GBWFILE\"\n\n")
fi

XYZ=$(cat $XYZFILE.tmp)

cat >> $TARGET << EOF
!uks $FUNCTIONAL_FREQ $BS_S_FREQ grid4 tightscf
!rijcosx $BS_AUX_J_S_FREQ $BS_AUX_C_S_FREQ gridx4
!sp xyzfile
$MOINP
%pal nprocs $NPROCS_DISTURB end

%base "${INPUTFILE_FREQ_OLD%.???}"

%basis
$BASIS_FREQ
  end

%scf
  maxiter 500
  $BROKENSYM
  end

*xyz $CHARGE $SPIN
$XYZ
*

EOF

#ADDING JOBS TO DAG

cat >> $PROJECT.dag << EOF
JOB freq_sp-$j $PROJECT.old.sub
VARS freq_sp-$j PREFIX="$BASE_FREQ_OLD"
SCRIPT POST freq_sp-$j $PROJECT.old-check.sh $BASE_OLD.$j
RETRY freq_sp-$j 5

EOF
for i in $(seq 0 1 $SIX_N_FROM_0); do

BASE_FREQ_DISTURB=$BASE_DISTURB.$i.$j
INPUTFILE_FREQ_DISTURB=$BASE_FREQ_DISTURB.inp
XYZFILE_FREQ_DISTURB=$BASE_FREQ_DISTURB.xyz

#CREATING DISTURBED XYZFILES

ROW_i=$(echo "$(expr $i / 6) + 1" | bc)
COLUMN_i=$(echo "$(expr $(expr $i / 2) % 3) + 2" | bc)
NUM_OLD_i=$(sed "${ROW_i}q;d" $XYZFILE.tmp | awk "{print \$${COLUMN_i}}")
SIGN_i=$(expr $i % 2)

if [ $SIGN_i -eq 0 ]; then
  NUM_NEW_i=$(echo "${NUM_OLD_i}+0.005*0.529177249" | bc -l)
else
  NUM_NEW_i=$(echo "${NUM_OLD_i}-0.005*0.529177249" | bc -l)
fi

XYZ=$(awk "NR==${ROW_i}{\$${COLUMN_i}=${NUM_NEW_i}}1" ${XYZFILE}.tmp)

#CREATING INPUTFILES

TARGET=$INPUTFILE_FREQ_DISTURB

MOINP=""

if [ -f $PREFIX_OLD.gbw ]; then
  MOINP=$(echo -e "!moread\n\n%moinp \"$GBWFILE\"\n\n")
fi

cat >> $TARGET << EOF
!uks $FUNCTIONAL_FREQ $BS_S_FREQ grid4 tightscf
!rijcosx $BS_AUX_J_S_FREQ $BS_AUX_C_S_FREQ gridx4
!opt xyzfile
$MOINP
%pal nprocs $NPROCS_DISTURB end

%base "${INPUTFILE_FREQ_DISTURB%.???}"

%basis
$BASIS_FREQ
  end

%scf
  maxiter 500
  $BROKENSYM
  end

%geom
  maxiter 1
  end

*xyz $CHARGE $SPIN
$XYZ
*

EOF

#ADDING JOBS TO DAG

cat >> $PROJECT.dag << EOF
JOB freq_disturb-$i-$j $PROJECT.disturb.sub
VARS freq_disturb-$i-$j PREFIX="$BASE_FREQ_DISTURB"
SCRIPT POST freq_disturb-$i-$j $PROJECT.disturb-check.sh $BASE_DISTURB.$i
RETRY freq_disturb-$i-$j 5

EOF

done
done

#PREPARING FOR JOB1 SUBMISSION

cat >> $PROJECT.old.sub << EOF
universe = vanilla
log = $PROJECT.disturb-\$(Cluster).\$(Process).logfile
error = $PROJECT.disturb-\$(Cluster).\$(Process).errfile
transfer_input_files = http://proxy.chtc.wisc.edu/SQUID/tyang29/orca-3.0.3-grad-only.tar.gz,\$(PREFIX).inp, $GBWFILE

executable = $PROJECT.sh
arguments = \$(PREFIX).inp

output = $PROJECT.disturb-\$(Cluster).\$(Process).outfile

should_transfer_files = YES
when_to_transfer_output = ON_EXIT

notification = Always
notify_user = tyang29@wisc.edu

request_cpus = $NPROCS_DISTURB
request_memory = $((2*$NPROCS_DISTURB))GB
request_disk = $((30*$NPROCS_DISTURB))GB


+WantGlideIn = true
periodic_release = (JobStatus == 5) && ((CurrentTime - EnteredCurrentStatus) > 180) && (HoldReasonCode == 12 || HoldReasonCode == 13) && (NumJobStarts < 5)

queue

EOF

#PREPARING FOR JOB1 SUBMISSION

cat >> $PROJECT.disturb.sub << EOF
universe = vanilla
log = $PROJECT.disturb-\$(Cluster).\$(Process).logfile
error = $PROJECT.disturb-\$(Cluster).\$(Process).errfile
transfer_input_files = http://proxy.chtc.wisc.edu/SQUID/tyang29/orca-3.0.3-grad-only.tar.gz,\$(PREFIX).inp, $GBWFILE

executable = $PROJECT.sh
arguments = \$(PREFIX).inp

output = $PROJECT.disturb-\$(Cluster).\$(Process).outfile

should_transfer_files = YES
when_to_transfer_output = ON_EXIT

notification = Always
notify_user = tyang29@wisc.edu

request_cpus = $NPROCS_DISTURB
request_memory = $((2*$NPROCS_DISTURB))GB
request_disk = $((30*$NPROCS_DISTURB))GB


+WantGlideIn = true
periodic_release = (JobStatus == 5) && ((CurrentTime - EnteredCurrentStatus) > 180) && (HoldReasonCode == 12 || HoldReasonCode == 13) && (NumJobStarts < 5)

queue

EOF

#########################################################
#CREATING POSTSCRIPT TO CHECK DISTURB FREQUENCY HESSFILES
#########################################################

OLD_CHECK=$PROJECT.old-check.sh

cat >> $OLD_CHECK << EOF
#!/bin/sh

rm *tmp*
EXIT_CODE=1
for j in \$(find ./ -name "\$1.log" -exec basename {} \;); do
while [ \$(grep -c 'SCF CONVERGED' \$j) -gt 0 ] && [ \$(grep -c 'TOTAL RUN TIME' \$j) -eq 1 ]; do
  find ./ -name "\${j%.???}.gbw" ! -name "$BASE_OLD.$i.*.gbw" -exec rm -rf {} \;
  EXIT_CODE=0
  sleep 3
  for i in \$(find ./ -name "$PROJECT.dag.dagman.out" -exec grep "ULOG_SUBMIT for HTCondor Node freq_old-\$(echo \$j | cut -d. -f3)" {} \; | sed 's/^.*(\(.*\)\.0).*/\1/'); do
    condor_rm \$i
  done
  cp \$j \$(echo \$j | sed 's/\(.*\)\..*\.log/\1/').2.log
  break
done
done

echo \$EXIT_CODE
exit \$EXIT_CODE

EOF

chmod +x "$OLD_CHECK"

#########################################################
#CREATING POSTSCRIPT TO CHECK DISTURB FREQUENCY HESSFILES
#########################################################

DISTURB_CHECK=$PROJECT.disturb-check.sh

cat >> $DISTURB_CHECK << EOF
#!/bin/sh

rm *tmp*
EXIT_CODE=1
for j in \$(find ./ -name "\$1.*.log" -exec basename {} \;); do
while [ \$(echo "\$(grep -c 'The current gradient' \${j%.???}.engrad) + \$(grep -c 'SCF CONVERGED' \$j)" | bc) -gt 1 ] && [ \$(grep -c 'TOTAL RUN TIME' \$j) -eq 1 ]; do
  find ./ -name "\${j%.???}.gbw" ! -name "$BASE_DISTURB.$i.*.gbw" -exec rm -rf {} \;
  EXIT_CODE=0
  sleep 3
  for i in \$(find ./ -name "$PROJECT.dag.dagman.out" -exec grep "ULOG_SUBMIT for HTCondor Node freq_disturb-\$(echo \$j | cut -d. -f3)-" {} \; | sed 's/^.*(\(.*\)\.0).*/\1/'); do
    condor_rm \$i
  done
  cp \${j%.???}.engrad \$(echo \$j | sed 's/\(.*\)\..*\.log/\1/').2.engrad
  cp \$j \$(echo \$j | sed 's/\(.*\)\..*\.log/\1/').2.log

NUM_ATOMS=\$(sed '1q;d' \$1.*.xyz)
NUM_ATOMS_N_TWO=\$(echo "\$NUM_ATOMS + 2" | bc)
DOF=\$(echo "\$NUM_ATOMS * 3" | bc)
DOF_N_ONE=\$(echo "\$NUM_ATOMS * 3 + 1" | bc)
DOF_MINUS_ONE=\$(echo "\$NUM_ATOMS * 3 - 1" | bc)
TWICE_DOF=\$(echo "\$NUM_ATOMS * 6" | bc)
TWICE_DOF_MINUS_ONE=\$(echo "\$NUM_ATOMS * 6 - 1" | bc)

while [ \$(find ./ -name "*grad.*.2.log" | wc -l) -eq \${TWICE_DOF} ]; do

rm hessian.txt

FINAL_HESS=\$(echo \$1 | cut -d. -f1).final.hess
FINAL_THERMO=\$(echo \$1 | cut -d. -f1).thermo.txt

TOTAL_COLUMNS=\$(expr \$NUM_ATOMS / 2)

while [ \$(expr \$NUM_ATOMS % 2) -eq 1 ]; do
  TOTAL_COLUMNS=\$(echo "\$(expr \$NUM_ATOMS / 2) + 1" | bc)
  break
done

COLUMNS=6

for column in \$(seq 1 1 \$TOTAL_COLUMNS); do
  rm row_label.txt
while [ \$(echo "\$column * 12" | bc) -gt \$TWICE_DOF ]; do
  COLUMNS=3
  break
done
END_i=\$(echo "\$column * 12 - 1" | bc)
while [ \$END_i -gt \$TWICE_DOF ]; do
  END_i=\$(echo "\$TWICE_DOF - 1" | bc)
  break
done
for i in \$(seq \$(echo "\$column * 12 - 12" | bc) 2 \$END_i); do
  expr \$i / 2 >> row_label.txt
  grep -A\${DOF_N_ONE} 'The current gradient in Eh' \$(echo \$1 | sed 's/\(.*\)\..*/\1/').\$i.2.engrad | tail -n \${DOF} | awk '{print \$1*100}' > grad+.txt
  grep -A\${DOF_N_ONE} 'The current gradient in Eh' \$(echo \$1 | sed 's/\(.*\)\..*/\1/').\$(echo "\$i + 1" | bc).2.engrad | tail -n \${DOF} | awk '{print \$1*100}' > grad-.txt
  awk '{getline t<"grad+.txt"; print t-\$0}' grad-.txt > hessian.\$(echo "\$i+1" | bc).txt
  seq 0 1 \$(wc -l grad+.txt | awk '{print \$1-1}') > column_label.txt
done
cat row_label.txt | pr -ts" " --column \${COLUMNS} >> hessian.txt
while [ \$(echo "\$column * 12" | bc) -gt \$TWICE_DOF ]; do
  paste column_label.txt hessian.\$(echo "\$column * 12 - 11" | bc).txt hessian.\$(echo "\$column * 12 - 9" | bc).txt hessian.\$(echo "\$column * 12 - 7" | bc).txt >> hessian.txt
  break
done
while [ \$(echo "\$column * 12" | bc) -le \$TWICE_DOF ]; do
paste column_label.txt hessian.\$(echo "\$column * 12 - 11" | bc).txt hessian.\$(echo "\$column * 12 - 9" | bc).txt hessian.\$(echo "\$column * 12 - 7" | bc).txt hessian.\$(echo "\$column * 12 - 5" | bc).txt hessian.\$(echo "\$column * 12 - 3" | bc).txt hessian.\$(echo "\$column * 12 - 1" | bc).txt >> hessian.txt
  break
done
done

if [ -f dipoledr.txt ]; then
  rm dipoledr.txt
fi

for i in \$(seq 0 2 \${TWICE_DOF_MINUS_ONE}); do
  grep -m1 'Total Dipole' \$(echo \$1 | sed 's/\(.*\)\..*/\1/').\$i.2.log | awk '{print \$5*100,\$6*100,\$7*100}' | sed 's/ /\n/g' > dipole+.txt
  grep -m1 'Total Dipole' \$(echo \$1 | sed 's/\(.*\)\..*/\1/').\$(echo "\$i + 1" | bc).2.log | awk '{print \$5*100,\$6*100,\$7*100}' | sed 's/ /\n/g' > dipole-.txt
  awk '{getline t<"dipole+.txt"; print t-\$0}' dipole-.txt | xargs echo >> dipoledr.txt
done

rm dipole+.txt
rm dipole-.txt

rm row_label.txt
rm column_label.txt
rm grad+.txt
rm grad-.txt
rm grad_total.txt
rm hess.txt

sed '1,2d' \$(echo \$1 | cut -d. -f1).sp.1.xyz | awk '{print \$1,\$2/0.529177,\$3/0.529177,\$4/0.529177}' > \$(echo \$1 | cut -d. -f1).sp.mass.xyz

for i in \$(seq 1 1 \${NUM_ATOMS}); do
  element=\$(awk "NR==\${i}{print \\\$1}" \$(echo \$1 | cut -d. -f1).sp.mass.xyz)
  mass=\$(grep "\${element} " atomic_mass.txt | awk '{print \$2}')
  sed -i "\${i}s/\${element} /\${element} \${mass} /" \$(echo \$1 | cut -d. -f1).sp.mass.xyz
done

if [ -f job_list.txt ]; then
  rm job_list.txt
fi

for i in \$(seq 0 1 \${DOF_MINUS_ONE}); do
  echo "   \${i}  1  1  1" >> job_list.txt
done

echo -e "\n\\\$orca_hessian_file\n\n\\\$act_atom\n  \$((\${NUM_ATOMS} - 1))\n\n\\\$act_coord\n  2\n\n\\\$act_energy\n     \$(tac \$(echo \$1 | cut -d. -f1).sp.2.log | grep -m1 'FINAL SINGLE POINT ENERGY' | awk '{print \$5}')\n\n\\\$hessian\n\${DOF}\n\n#\n# The atoms: label  mass x y z\n#\n\\\$atoms\n\${NUM_ATOMS}\n\$(cat \$(echo \$1 | cut -d. -f1).sp.mass.xyz)\n\n\\\$actual_temperature\n  0.000000\n\n\\\$dipole_derivatives\n\${DOF}\n\$(cat dipoledr.txt)\n\n\\\$job_list\n\${DOF}\n\$(cat job_list.txt)\n\n\\\$end" > \$(echo \$1 | cut -d. -f1).final.hess

sed -ie "14r hessian.txt" \$FINAL_HESS

wget http://proxy.chtc.wisc.edu/SQUID/tyang29/orca3_vib
chmod +x orca3_vib
./orca3_vib \$FINAL_HESS
sleep 3

rm ./orca3_vib
rm hessian*

#Thermo

#Define fundamental constants

exp=2.71828182846 # Natural exponential
ratio=349.755 # 350 cm-1 = 1 kcal*mol-1
TEMP=298.15 # in Kelvin
P=101325 # in pascal
k_B=0.00198720628 # in kcal*mol-1*K-1
N_A=\$(echo '6.022149*10^23' | bc -l) # Avogadro number
h=\$(echo '3.9903175*10^-10' | bc -l) # in J*s*mol-1
pi=3.1415926

#Calculate the rotational partition function

NUM_ATOMS=\$(grep -A 1 '\\\$atoms' \${FINAL_HESS} | tail -n 1) # number of atoms in the molecule
grep -A \$((\${NUM_ATOMS}+2)) '\\\$atoms' \${FINAL_HESS} | sed '1,2d' > mass_coord.xyz
M=\$(awk '{total = total + \$2}END{print total}' mass_coord.xyz) # mass of the molecule
x_center=\$(awk "{total = total + \\\$2*\\\$3}END{print total/\${M}}" mass_coord.xyz)
y_center=\$(awk "{total = total + \\\$2*\\\$4}END{print total/\${M}}" mass_coord.xyz)
z_center=\$(awk "{total = total + \\\$2*\\\$5}END{print total/\${M}}" mass_coord.xyz)

#echo "center of mass at \${x_center} \${y_center} \${z_center}"

for i in \$(seq 1 1 \${NUM_ATOMS}); do
  m=\$(awk "NR==\${i}{print \\\$2}" mass_coord.xyz)
  x=\$(awk "NR==\${i}{print \\\$3 - \${x_center}}" mass_coord.xyz)
  y=\$(awk "NR==\${i}{print \\\$4 - \${y_center}}" mass_coord.xyz)
  z=\$(awk "NR==\${i}{print \\\$5 - \${z_center}}" mass_coord.xyz)
  echo "\${m}*(\${y}^2+\${z}^2)" | bc -l >> I_xx.txt
  echo "\${m}*\${x}*\${y}" | bc -l >> I_xy.txt
  echo "\${m}*\${x}*\${z}" | bc -l >> I_xz.txt
  echo "\${m}*(\${x}^2+\${z}^2)" | bc -l >> I_yy.txt
  echo "\${m}*\${y}*\${z}" | bc -l >> I_yz.txt
  echo "\${m}*(\${x}^2+\${y}^2)" | bc -l >> I_zz.txt
done

I_xx=\$(awk '{total = total + \$1}END{print total}' I_xx.txt)
I_xy=\$(awk '{total = total + \$1}END{print total}' I_xy.txt)
I_xz=\$(awk '{total = total + \$1}END{print total}' I_xz.txt)
I_yx=\${I_xy}
I_yy=\$(awk '{total = total + \$1}END{print total}' I_yy.txt)
I_yz=\$(awk '{total = total + \$1}END{print total}' I_yz.txt)
I_zx=\${I_xz}
I_zy=\${I_yz}
I_zz=\$(awk '{total = total + \$1}END{print total}' I_zz.txt)

echo -e "\${I_xx},\${I_xy},\${I_xz}\n\${I_yx},\${I_yy},\${I_yz}\n\${I_zx},\${I_zy},\${I_zz}" > I_matrix.txt
I_tot=\$(echo "" | ./eig.py) # in (g*Ang*mol-1)^3

qrot=\$(echo "\${I_tot}" | awk "{print sqrt(\\\$1*\${pi}/(1000)^3/(10^10)^6)*((8*(\${pi}^2)*\${k_B}*\${TEMP}*4184/(\${h}^2))^(3/2))}")

rm I_*.txt

#Calculate the translational partition function

qtran=\$(echo \${M} | awk "{print (2*\${pi}*\\\$1/1000*\${k_B}*\${TEMP}*4184/(\${h}^2))^(3/2)*\${k_B}*\${TEMP}*4184/\${N_A}/\${P}}")

echo \${qtran}

#Calculate vibrational terms

DOF=\$(grep -A 1 'vib' \${FINAL_HESS} | tail -n 1)
grep -A \$((\${DOF} + 1)) 'vib' \${FINAL_HESS} | tail -n \$((\${DOF} - 6)) | awk '\$2 > 0{print \$2}' > vib.txt

ZPE=\$(awk "{total = total + \\\$1}END{print total/\${ratio}/2}" vib.txt)
E_vib=\$(awk "{print \\\$1/\${ratio}*\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B})/(1-\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B}))}" vib.txt | awk '{total = total + \$1}END{print total}')
E_rot=\$(echo "3/2*\${TEMP}*\${k_B}" | bc -l)
E_tran=\$(echo "3/2*\${TEMP}*\${k_B}" | bc -l)
E_tot=\$(echo "\${ZPE}+\${E_vib}+\${E_rot}+\${E_tran}" | bc -l)
TS_ele=\$(echo \${FINAL_HESS} | cut -d_ -f3 | awk "{print \${k_B}*\${TEMP}*log(\\\$1)/log(\${exp})}")
#TS_vib=\$(awk "{print -1*\${k_B}*\${TEMP}*log(1-\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B}))/log(\${exp})+\\\$1/\${ratio}*\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B}/(1-\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B})))}" vib.txt | awk '{total = total + \$1}END{print total}')
TS_vib=\$(awk "{print -1*\${k_B}*\${TEMP}*log(1-\${exp}^(-1*\\\$1/\${ratio}/\${TEMP}/\${k_B}))/log(\${exp})+\\\$1/\${ratio}/(\${exp}^(\\\$1/\${ratio}/\${TEMP}/\${k_B})-1)}" vib.txt | awk '{total = total + \$1}END{print total}')
#if [ \${qrot} -eq 0 ]; then
#  TS_rot=\$(echo "3/2" | awk "{print \${k_B}*\${TEMP}*\\\$1}")
#else
TS_rot=\$(echo "\${qrot}" | awk "{print \${k_B}*\${TEMP}*(log(\\\$1)/log(\${exp})+3/2)}")
#  fi
TS_tran=\$(echo "\${qtran}" | awk "{print \${k_B}*\${TEMP}*((log(\\\$1)/log(\${exp}))+3/2)}")
TS_tot=\$(echo "\${TS_ele}+\${TS_vib}+\${TS_rot}+\${TS_tran}" | bc -l)

if [ -f \${FINAL_THERMO} ]; then
  rm \${FINAL_THERMO}
fi

echo "ZPE                                        \$(printf "%.6g" \${ZPE})" >> \${FINAL_THERMO}
echo "E_vib                                      \$(printf "%.6g" \${E_vib})" >> \${FINAL_THERMO}
echo "E_rot                                      \$(printf "%.6g" \${E_rot})" >> \${FINAL_THERMO}
echo "E_tran                                     \$(printf "%.6g" \${E_tran})" >> \${FINAL_THERMO}
echo "Total correction                           \$(printf "%.6g" \${E_tot})" >> \${FINAL_THERMO}
echo "TS_ele                                     \$(printf "%.6g" \${TS_ele})" >> \${FINAL_THERMO}
echo "TS_vib                                     \$(printf "%.6g" \${TS_vib})" >> \${FINAL_THERMO}
echo "TS_rot                                     \$(printf "%.6g" \${TS_rot})" >> \${FINAL_THERMO}
echo "TS_tran                                    \$(printf "%.6g" \${TS_tran})" >> \${FINAL_THERMO}
echo "Final entropy term                ...      \$(printf "%.6g" \${TS_tot})" >> \${FINAL_THERMO}

#rm vib.txt

#Notify user

echo "$PROJECT is complete" | mail -s "$PROJECT is complete" tyang29@wisc.edu

break

done

break

done
done

echo \$EXIT_CODE
exit \$EXIT_CODE

EOF

chmod +x "$DISTURB_CHECK"

condor_submit_dag $PROJECT.dag

cd $FileDir
