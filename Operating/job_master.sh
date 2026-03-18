



#SBATCH --ntasks-per-node=1















#SBATCH --cpus-per-task=1















#SBATCH --time=1:59:59





for ((i=101;i<=200;i++));do


        mkdir "td$i"

        cd "td$i"

        cp ../control.inp .

        cp ../job_script.sh .

        sed -i "s/ion_velocity_init_seed=3300/ion_velocity_init_seed=$((i+3300))/" control.inp

	    sed -i "s|#SBATCH --job-name=propane_td1|#SBATCH --job-name=propane_td$i|" job_script.sh

        sbatch job_script.sh	

	cd ..


done
