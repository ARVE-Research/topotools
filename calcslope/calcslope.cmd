### Job name
#PBS -N calcslope_world

### Declare job non-rerunable
#PBS -r n
#PBS -k oe

###  Send email when job aborts or job end
#PBS -m ae

### queue name
#PBS -q debug

### wall time
#PBS -l walltime=00:30:00

### nodes and cores
#PBS -l nodes=1:ppn=20

###############################################################################
#The following stuff will be executed in the first allocated node.            #
#Please don't modify it                                                       #
#                                                                             #
echo $PBS_JOBID : `wc -l < $PBS_NODEFILE` CPUs allocated: `cat $PBS_NODEFILE`
PATH=$PBS_O_PATH
JID=`echo ${PBS_JOBID}| sed "s/.hpc2015-mgt.hku.hk//"`
###############################################################################

echo ===========================================================
echo "Job Start  Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

module load mpc/1.0.3 mpfr/3.1.3 gcc/8.2.0 impi/2019.4 hdf5/impi/1.10.5 netcdf/impi/4.7.1 autoconf/2.69 automake/1.15

export OMP_NUM_THREADS=20

export workdir=/group/esd_kaplan/datasets/topography/topotools/calcslope
export datadir=/group/esd_kaplan/datasets/topography/MERIT-Hydro

$workdir/calcslope 110/115/22/27 $datadir/MERIT-Hydro.nc $workdir/schinav2.nc

echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

exit 0
