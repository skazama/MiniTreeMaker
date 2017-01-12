dir=`pwd`
run_file="run.list"
pax_ver="pax_v6.0.2"
ver="v0.8"

while read line
do
        input=`echo $line | awk '{print $1}'`

	echo "cd ${dir}"                                                                       > batch/exe_${input}.sh 
	echo "cd ../"                                                                         >> batch/exe_${input}.sh 

        echo "mkdir -p /scratch/$USER/\${PBS_JOBID}"                                          >> batch/exe_${input}.sh
        echo "mkdir -p /scratch/$USER/\${PBS_JOBID}/rootfile"                                 >> batch/exe_${input}.sh

        ls /home/atp/kazama/xenon1t/pax/${pax_ver}/${input}*.root                              > batch/input_${input}.txt

        echo "cp -r bin/ana data jobs/batch/input_${input}.txt /scratch/$USER/\${PBS_JOBID}"  >> batch/exe_${input}.sh
        echo "cd /scratch/$USER/\${PBS_JOBID}"                                                >> batch/exe_${input}.sh 

        echo "./ana input_${input}.txt rootfile/${input}.root"                                >> batch/exe_${input}.sh
        echo "cp rootfile/${input}.root /disk/data1/atp/xenon1t/processed/${pax_ver}/${ver}/" >> batch/exe_${input}.sh

	echo "rm -rf /scratch/$USER/\${PBS_JOBID}"                                            >> batch/exe_${input}.sh

done  < ${run_file}

for i in `ls batch/exe*sh`
do
	echo "qsub -l cput=24:00:00  $i"
done > job_submit.sh
