#PBS -V
#PBS -N hoagie
#PBS -o /dev/null
#PBS -e /dev/null

# going to PBS working directory
echo "going to PBS working directory: ${PBS_O_WORKDIR}"
cd ${PBS_O_WORKDIR}

# running program
echo "running program"
./total_dynamix

# deleting run script
echo "deleting run script"
rm -f /home/andyras/git/dynamix/hoagie.sh
