# Working with job scheduler

<center>![image](./nesi_images/scheduler_image.png){width="500"}</center>

## Introduction to slurm scheduler and directives

An HPC system might have thousands of nodes and thousands of users. How do we decide who gets what and when? How do we ensure that a task is run with the resources it needs? This job is handled by a special piece of software called the scheduler. On an HPC system, the scheduler manages which jobs run where and when. In brief, scheduler is a 

!!! quote ""

    * Mechanism to control access by many users to shared computing resources
    * Queuing / scheduling system for users’ jobs
    * Manages the reservation of resources and job execution on these resources 
    * Allows users to “fire and forget” large, long calculations or many jobs (“production runs”)

!!! circle-info "A bit more on why do we need a scheduler ?"

    * To ensure the machine is utilised as fully as possible
    * To ensure all users get a fair chance to use compute resources (demand usually exceeds supply)
    * To track usage - for accounting and budget control
    * To mediate access to other resources e.g. software licences

    **Commonly used schedulers**
    
    * [x] Slurm
    * PBS , Torque
    * Grid Engine
    
    All NeSI clusters use Slurm (Simple Linux Utility for Resource Management) scheduler (or job submission system) to manage resources and how they are made available to users.

    ![image](./nesi_images/slurm_comms2compute.png)

    <small>Researchers can not communicate directly to  Compute nodes from the login node. Only way to establish a connection OR send scripts to compute nodes is to use scheduler as the carrier/manager</small>


      
## Life cycle of a slurm job

<center>
![image](./nesi_images/batch_system_flow.png){width="1000"}
</center>

- - -
!!! clipboard-list "Commonly used Slurm commands"

    | Command        | Function                                                                                             |
    |:---------------|:------------------------------------------------------------------------------------------------------|
    | `sbatch`       | Submit non-interactive (batch) jobs to the scheduler                                                 |
    | `squeue`       | List jobs in the queue                                                                               |
    | `scancel`      | Cancel a job                                                                                         |
    | `sacct`        | Display accounting data for all jobs and job steps in the Slurm job accounting log or Slurm database|
    | `srun`         | Slurm directive for parallel computing                                                                      |
    | `sinfo`        | Query the current state of nodes                                                                     |
    | `salloc`       | Submit interactive jobs to the scheduler                                                             |

- - - 

??? question "Exercise 5.1 (Optional) - Check the state of the compute cluster " 

    * summary of current states of compute nodes known to the scheduler
    ```bash
    sinfo
    ```

    * similar to above but expanded
    ```bash
    sinfo --format="%16P %.8m %.5a %10T %.5D %80N"
    ```

    * will print a long output as it is one row per compute node in the cluster
    ```bash
    sinfo -N -l
    ```

     * Explore the capacity of a compute node
     ```bash
     sinfo -n wch001 -o "%n %c %m"
     ```


## Anatomy of a slurm script and submitting first slurm job 🧐

As with most other scheduler systems, job submission scripts in Slurm consist of a header section with the shell specification and options to the submission command (`sbatch` in this case) followed by the body of the script that actually runs the commands you want. In the header section, options to `sbatch` should be prepended with `#SBATCH`.


<center>
![image](./nesi_images/anatomy_of_a_slurm_script.png){width="700"}
</center>


!!! square-pen "Commented lines `#`"

    Commented lines are ignored by the bash interpreter, but they are not ignored by slurm. The `#SBATCH` parameters are read by slurm when we submit the job. When the job starts, the bash interpreter will ignore all lines starting with `#`. This is very similar to the shebang mentioned earlier, when you run your script, the system looks at the `#!`, then uses the program at the subsequent path to interpret the script, in our case `/bin/bash` (the program `bash` found in the */bin* directory

-  - - 

??? circle-info "Slurm variables"

    | header          | use                                 | description                                          |
    |:--------------- |:------------------------------------|:-----------------------------------------------------|
    |--job-name 	  | `#SBATCH --job-name=MyJob` 	        |The name that will appear when using squeue or sacct. |
    |--account 	      | `#SBATCH --account=nesi12345` 	    |The account your core hours will be 'charged' to.     |
    |--time 	      | `#SBATCH --time=DD-HH:MM:SS` 	    |Job max walltime.                                     |
    |--mem 	          | `#SBATCH --mem=512MB` 	            |Memory required per node.                             |
    |--cpus-per-task  | `#SBATCH --cpus-per-task=10` 	    |Will request 10 logical CPUs per task.                |
    |--output 	      | `#SBATCH --output=%j_output.out` 	|Path and name of standard output file. `%j` will be replaced by the job ID.         |
    |--mail-user 	  | `#SBATCH --mail-user=me23@gmail.com`|address to send mail notifications.                   |
    |--mail-type 	  | `#SBATCH --mail-type=ALL` 	        |Will send a mail notification at BEGIN END FAIL.      |
    |                 | `#SBATCH --mail-type=TIME_LIMIT_80` |Will send message at 80% walltime.                    |

<br>

??? bell "Assigning values to Slurm variables"
    <center>![image](./nesi_images/sbtach_def_1.png)</center>

    <center>![image](./nesi_images/sbatch_def_2.png)</center>
- - -

??? question "Exercise 5.2"

    Let's put these directives together and compile our first slurm script. Below is a abstract version of the slurm life cycle to assist you with the process

    <center>![image](./nesi_images/slurm_cycle_mini.png)</center>
    
    * First step is to create new working directories inside the existing `~/scripting_workshop/scheduler` 
    ```bash
     cd ~/scripting_workshop/scheduler
    ```

    * confirm the path is correct 
    ```bash
    pwd
    ```

    * create a new directory for this section and change the directory to it - Check for the follow up not `&&`
    ```bash
    mkdir ex_5.2 && cd ex_5.2
    ```
    !!! quote ""
        The meaning of `&&` and `&` are intrinsically different.
      
        * **What is `&&` in Bash?** In Bash—and many other programming languages—`&&` means “AND”. And in command execution context like this, it means items to the left as well as right of && should be run in sequence in this case.
        * **What is & in Bash?** And a single `&`means that the preceding commands—to the immediate left of the &—should simply be run in the background.

    * use a text editor of choice to create a file named firstslurm.sl - we will use nano here
    ```bash
    nano firstslurm.sl
    ```

    * Content of `firstslurm.sl` should be as below. Please discuss as you make progress
    ```bash
    #!/bin/bash 

    #SBATCH --job-name      myfirstslurmjob
    #SBATCH --account       nesi02659
    #SBATCH --time          00:02:00                 #Format is DD-HH:MM:SS
    #SBATCH --cpus-per-task 1
    #SBATCH --mem           512                      #Default unit is Megabytes
    #SBATCH --output        slurmjob.%j.out
    #SBATCH --error         slurmjob.%j.err

    sleep 100

    echo "I am a slurm job and I slept for 100 seconds"
  
    echo "$SLURM_JOB_ID END"
    ```

    * **Save** and **Exit**
    * Submit the script with `sbatch` command
    ```bash
    sbatch firstslurm.sl
    ```
    *  Execute `squeue --me` and `sacct`. Discuss the outputs .i.e.
    ```bash
    squeue --me
    ```
    ```bash
    sacct
    ```
    ??? hand-holding-dollar "`$SLURM_JOB_ID`"
        `$SLURM_JOB_ID` is a Slurm environment variable. 
        - A full list of environment variables for SLURM can be found by visiting the [SLURM page on environment variables](https://slurm.schedmd.com/sbatch.html#SECTION_INPUT-ENVIRONMENT-VARIABLES)
        - These variables are great for recursive operations. 
- - - 
 
!!! magnifying-glass "STDOUT/STDERR from jobs"

    * STDOUT - your process writes conventional output to this file handle
    * STDERR - your process writes diagnostic output to this file handle.
    
    **STDOUT** and **STDERR** from jobs are, by default, written to a file called `slurm-JOBID.out` and `slurm-JOBID.err` in the working directory for the job (unless the job script changes this, this will be the directory where you submitted the job). So for a job with ID 12345 STDOUT and STDERR will be `slurm-12345.out` and `slurm-12345.err`

    - When things go wrong, first step of **debugging** (STORY TIME !) starts with a referral to these files. 



## Assessing resource utilisation (cpu, memory, time)

Understanding the resources you have available and how to use them most efficiently is a vital skill in high performance computing. The three resources that every single job submitted on the platform needs to request are:

* CPUs (i.e. logical CPU cores), and
* Memory (RAM), and
* Time.

***What happens if I ask for the wrong resources?***

---

| Resource         | Asking for too much                                   | Not asking for enough                                                               |
|:---------------  |:------------------------------------------------------|:------------------------------------------------------------------------------------|
| Number of CPUs   | Job may wait in the queue for longer   	           | Job will run more slowly than expected, and so may run out time    |                    
|                  | Drop in fairshare score which determines job priority |                                                                                     |
| Memory           | (above)                                               | Job will fail, probably with `OUT OF MEMORY` error, segmentation fault or bus error|
| Wall time        | (above)                                               | Job will run out of time and get killed                                             |

---

??? question "Exercise 5.3" 


    Let's submit another slurm job and review its resource utilisation

    * Change the working directory to Exercise_5.3
    ```bash
    cd ~/scripting_workshop/scheduler/ex_5.3
    ```

    * Run `ls` command and you should see two files (one .R and one sl) and one directory named slurmout
    ```bash
    ls -F
    ```
      ```bash
      bowtie-test.sl*  input_data/  slurmout/
      ```
    * Review the slurm script bowtie-test.sl with nano and edit the corresponding sections (hint :email)
    ```bash
    sbatch bowtie-test.sl 
    ```

    * use `squeue --me` and `sacct` again to evaluate the job status

    * Once the job ran into completion, use `nn_seff JOBID` command to print the resource utilisation statistics (Replace **JOBID** with the corresponding number)
    ```bash
    $ nn_seff 25222190
    Job ID: 25222190
    Cluster: mahuika
    User/Group: me1234/me1234
    State: COMPLETED (exit code 0)
    Cores: 1
    Tasks: 1
    Nodes: 1
    Job Wall-time:  18.33%  00:00:33 of 00:03:00 time limit
    CPU Efficiency: 93.94%  00:00:31 of 00:00:33 core-walltime
    Mem Efficiency: 1.33%  13.62 MB of 1.00 GB
    ```
    * Now review the content of `.err` and `.out` files in */slurmout* directory

 

!!! quote ""
    **Feeling adventurous 🤠 ?** - Refer to [Supplementary material on slurm profiling](https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/7_supplementary_2/)


- - -

## Compiling slurm scripts for Variant Calling and RNA-seq episodes 

??? question "Exercise 5.4 😬"	
    
    Purpose of this exercise is to compile a slurm submission script based on the script we wrote in [episode 2 - Automating variant calling workflow](https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/2_AutomaticVariantC.html)

    * recommend creating a new directory for the exercise .i.e `ex_5.4`
    * Name of the file is `variant_calling.sl` (note that we have change the extension from `.sh` to `.sl`)
    * In terms of slurm variables

        * name of the job is `variant_calling_workflow`
        * number of CPUS is `2`
        * timelimit `15 minutes`
        * amount of memory in GB `4G`
        * generate  *.err* files and *.out* where both should be re-directed to the directory ***slurmout***
        * an email notification at the end of the job 

    * We don't want to replicate ***input data***  in multiple places .i.e. be conservative in-terms how you use research storage
    * Therefore, use the same reference genome file (assign the filename to variable `genome` and the trimmed read files (assign the path of these files to variable `trimmed`) used in the first episode
    ```bash
    genome=~/scripting_workshop/variant_calling/ref_genome/ecoli_rel606.fasta
    trimmed=~/scripting_workshop/variant_calling/trimmed_reads
    ```
    
- - - 

??? question  "Exercise 5.5 😬"	

    *  Now it's your turn to compile a slurm submission script for the RNA-seq workflow. 😊

---


<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>
