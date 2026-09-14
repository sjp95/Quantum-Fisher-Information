#!/bin/bash
sub=job.sge

#module load compilers/gcc/9.2.0 #||=========== Needed to run c++17 or higher version
#g++ -std=c++17 U_O.cpp -o h -llapack -I /usr/include/lapacke \-L /usr/lib64 #||=========== compilation with c++17

#g++ -std=c++11 U_O.cpp -o h -llapack -I /usr/include/lapacke \-L /usr/lib64

# for U in 4000
#  do
#   for A in 500 2200 2400 3000 3500
#    do
g++ -std=c++11 mag.cpp -o h 
    for i in {0..99}
     do
          #===============================================================================#
          #===============================================================================#
          echo "#!/bin/bash"                                                    >$sub
          echo "#SBATCH  --job-name=SR_"$A"_"$i"    # Job name"                >> $sub
          echo "#SBATCH --partition=serial          # Partition Name (queue)"  >> $sub
          echo "#SBATCH --qos=serial"                                          >> $sub
          echo "#SBATCH --ntasks-per-node=1       # No. of CPU in one node"    >> $sub
          echo "#SBATCH --nodes=1                  # No. of Node for jobs"     >> $sub
          echo "#SBATCH --output=out.log   # Standard output and error log"    >> $sub
          cd ./
          echo "echo "job started at date""                                    >> $sub
          echo "date"                                                          >> $sub
          echo "./h $i "                                                       >> $sub
          echo "echo ----"                                                     >> $sub
          echo "echo "job ended at date""                                      >> $sub
          echo "date"                                                          >> $sub
          #===============================================================================#
          #===============================================================================#
          chmod 777 job.sge
          sbatch  job.sge
  done
#  done
# done


