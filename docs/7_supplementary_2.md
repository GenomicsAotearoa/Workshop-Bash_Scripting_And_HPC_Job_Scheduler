# S2 : slurm profiling


Although `nn_seff` command is a quick and easy way to determine the resource utilisation, it relies on **peak** values (data gets recorded every 30 seconds) which doesn't allows us to examine resource usage over the run-time of the job. There are number of in-built/external tools to achieve the latter which will require some effort to understand its deployment, tracing and interpretation. Therefore, we will use **slurm native profiling** to evaluate resource usage over run-time. This is a simple and elegant solution.

???+ question "Exercise S.2.1"

    * Download and decompress the content
    ```bash

    wget -c Exercise_S21.tar.gz https://github.com/DininduSenanayake/nesi_introductory_custom/releases/download/v1.0/Exercise_S21.tar.gz -O - | tar -xz
    ```
    ```bash
    cd Exercise_S21
    ```

    * Run `ls` command and you should see three files (one .R,sl and one .py - We will discuss the purpose of this .py file after submitting the job) and one directory named slurmout
    ```bash
    ls -F
    ```
    ```bash
    example1_arraysum.R  example1_arraysum.sl  profile_plot_Jul2020.py  slurmout/
    ```

    * Review the slurm script with cat Or another text editor and submit with sbatch
    ```bash
    sbatch example1_arraysum.sl 
    ```
    
    * Do take a note of the **JOBID** as we are going to need it for next step. Otherwise, we use `squeue --me` OR `sacct` command as before to monitor the status
    * Also, you can `watch` the status of this job via `$ watch -n 1 -d "squeue -j JOBID"`. 
    *  `watch` command execute a program periodically, showing output fullscreen. Exiting the `watch` screen by done by pressing `Ctrl+x` 

    * Let's create slurm profile graphs

        - collate the data into an HDF5 file using the command. Replace **JOBID** with the corresponding number 

        ```bash
        sh5util -j JOBID
        ```
        ```bash
        sh5util: Merging node-step files into ./job_JOBID.h5
        ```

        - execute the script on .h5 file. We will need one of the Python 3 modules to do this. Ignore the deprecating warning. 
        ```bash
        module purge 
        ```
        ```bash
        module load Python/3.8.2-gimkl-2020a
        ```
        - Replace **JOBID** with the corresponding number
        ```bash
        python profile_plot_Jul2020.py job_JOBID.h5
        ```

       - This should generate a .png file where the filename is in the format of job_23258404_profile.png


    <br>
    ![image](./nesi_images/slurm_profile.png){width="1000"}
    </br>




<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>