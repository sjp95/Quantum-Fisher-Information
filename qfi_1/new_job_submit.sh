sub=job.sge
module load compilers/parallel_studio_xe_2018_update1/mkl/default 
module load compilers/parallel_studio_xe_2018_update1/parallel_studio_xe_2018/default
module load openmpi-3.0.0

export MKL_NUM_THREADS=32
export OMP_NUM_THREADS=32
g++ -std=c++11 -O3 st1.cpp -o main \
  -I${MKLROOT}/include \
  -L${MKLROOT}/lib/intel64 \
  -Wl,--start-group \
    -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread \
  -Wl,--end-group \
  -liomp5 -lpthread -lm -ldl

for U in 10 
do

     #********************************************************************************
     echo "#PBS  -N   Slave_rotor_"$U                          > $sub
     echo "#PBS -l nodes=n4:ppn=32"                            >> $sub
     echo "#PBS -j oe "                                        >> $sub
     echo "#PBS -o out.log"                                    >> $sub
     echo "#PBS -e err.log"                                    >> $sub
     echo "module load compilers/parallel_studio_xe_2018_update1/compiler/default"     >> $sub
     echo "module load compilers/parallel_studio_xe_2018_update1/mkl/default"     >> $sub
     echo "cd"    \$PBS_O_WORKDIR                               >> $sub
     echo "date"                                               >> $sub
     echo "./main $U"                                          >> $sub
     echo "date"                                               >> $sub
     #******************************************************************************
     chmod 777 job.sge
     qsub  job.sge

done
