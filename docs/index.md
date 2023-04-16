<center>![image](./nesi_images/nesiga_LOGO.png){width="300"}</center>

<center>
# **Introduction to Bash Scripting and HPC Scheduler**
 
</center>


<center>
![image](./nesi_images/slurm_llinux_penguin_dna.png){width="450"}
</center>



!!! check-to-slot "Prerequisites"

    - [x] Familiarity with terminal and basic linux commands
    - [x] Some knowledge on shell environment variables and `for` loops
    - [x] Ability to use a terminal based text editor such as `nano` 
        * [ ] This is not much of an issue as we are using JupyterHub which has a more friendlier text editor.  
    - [x] Intermediate level knowledge on Molecular Biology and Genetics

    **Recommended but not required**

    - [ ] Attend [Genomics Data Carpentry](https://datacarpentry.org/genomics-workshop/) and [RNA-Seq Data Analysis](https://genomicsaotearoa.github.io/RNA-seq-workshop/) workshops

<br>

!!! square-xmark "Some of the things we won't cover in this workshop"

     - Domain specific concepts in
         - Genetics, Genomics and DNA Sequencing 
         - Variant Calling (covers in [Genomics Data Carpentry Workshop](https://datacarpentry.org/genomics-workshop/))
         - RNA sequencing and data analysis ( covers in [RNA-Seq Data Analysis Workshop](https://genomicsaotearoa.github.io/RNA-seq-workshop/))
     

<br>

!!! screwdriver-wrench "Setup"

    Workshop material is designed to run on NeSI Mahuika cluster via Jupyter. Instructions on how to Set/Reset Authentication factors to access NeSI Services and Jupyter Login instructions [can be found here](https://dinindusenanayake.github.io/ganesi_authesetup-login/)

<br>

!!! rectangle-list "Content"
    | **Lesson**                                        | **Overview** | 
    |:---------------------------------------------------|:-------------|
    |[Background :fontawesome-solid-eye:](./0_Background.md){ .md-button .md-button--primary }|Objective of this workshop|
    |[1. Designing a Variant Calling Workflow  :fontawesome-solid-dna:](./1_DesigningVariantC.md){ .md-button .md-button--primary }|Develop and test the steps involved in calling variants|
    |[2. Automating a Variant Calling Workflow  :fontawesome-solid-wand-magic-sparkles:](./2_AutomaticVariantC.md){ .md-button .md-button--primary }|Compile a script based on the above steps|
    |[3. RNA-Seq Mapping And Count Data Workflow  :fontawesome-solid-map-location:](./3_RNAseq.md){ .md-button .md-button--primary }|Develop, test & compile a script for mapping reads and counting transcripts in a RNA-Seq data analysis pipeline|
    |[4. Introduction to HPC  :fontawesome-solid-server:](./4_IntroductiontoHPC.md){ .md-button .md-button--primary }|Introduction to High Performance Computing|
    |[5. Working with Job Scheduler  :fontawesome-solid-calendar-days:](./5_working_with_job_scheduler.md){ .md-button .md-button--primary }|Introduction to HPC Job Schedulers, Slurm Scheduler & life cycle of a Slurm job, Assessing resource utilisation and profiling|
    |[6. Supplementary #1](./6_supplementary_1.md){ .md-button .md-button--primary }||
    |[7. Supplementary #2](./7_supplementary_2.md){ .md-button .md-button--primary }||
    |[8. Supplementary #3](./8_supplementary_3.md){ .md-button .md-button--primary }||


:fontawesome-brands-twitter:{ .twitter }

:fontawesome-solid-lightbulb:{ .lightbulb }