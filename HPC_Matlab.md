# Running Matlab code at High Performance Computing (HPC) cluster at NYU.

Step-1: Use Git-Bash in Windows (or terminal in Linux) and [login to hpc](https://wikis.nyu.edu/display/NYUHPC/Logging+in+to+the+NYU+HPC+Clusters) as follow

*    From NYU
```
$ ssh <your-net-ID>@prince.hpc.nyu.edu
<your-net-ID>@prince.hpc.nyu.edu's password:
```
*    Outside NYU
```
$ ssh <your-net-ID>@gw.hpc.nyu.edu
<your-net-ID>@prince.hpc.nyu.edu's password:

$ ssh prince.hpc.nyu.edu
<your-net-ID>@prince.hpc.nyu.edu's password:
```

(Optional): Check your Quota
```
$ myquota

Filesystem   Environment   Backed up?   Allocation       Current Usage
Space        Variable      /Flushed?    Space / Files    Space(%) / Files(%)

/scratch    $SCRATCH       No/Yes        5.0TB/1.0M       0.00GB(0.00%)/4(0.00%)
/beegfs     $BEEGFS        No/Yes        2.0TB/3.0M        0.00GB(0.00%)/0(0.00%)
```

Step-2: Use `$SCRATCH` folder to run your commands from. It offers huge space 4Tb/user but get deleted in 60 days, so move all the content to `$ARCHIVE` after running the code. More details [here](https://wikis.nyu.edu/display/NYUHPC/Clusters+-+Prince).
```
$ cd $SCRATCH
$ pwd
```

Step-3: Clone the code from my repository 
```
$ git clone https://github.com/ishjain/mmWave.git
$ cd mmWave
```

(Optional): Edit the batch file. The content of [mybatch.sbatch](mybatch.sbatch) is pasted here. 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=ishjain@nyu.edu
#SBATCH --output=ish_%j.out
  
module purge
module load matlab/2017b

RUNDIR=$SCRATCH/mmWave/
cd $RUNDIR

matlab -nodisplay -nodesktop -r "run SimulationLOS.m" //Change the filename here
```
(Optional): Check the modules. You will find matlab/2017b
```
$ module avail
```

Step-4: Submit the job
```
$ sbatch mybatch.sbatch
Submitted batch job 4605025
```
(Optional): Check the job queue
```
$ squeue -u <your-net-ID>
JOBID   PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
4605025    c32_38   myTest   ikj211  R       0:14      1 c32-01
```
Step-5: Check the results in local directory
```
$ ls -lastr

```

Step-6: When running the same code many times (atmost 400). In the main file [SimulationLOS.m](SimulationLOS.m), we first get the environment variable `aID = getenv('SLURM_ARRAY_TASK_ID')`. The aID is different for different parallel computations. So, you can use these id to seed the random generator and use them to name the output file. Note: `aID` is an array, so first convert to int for seed.
```
sbatch --array=1-100 mybatch.sbatch
```
This above command runs the main code on 100 parallel machines with `aID` ranging from 1 to 100.
