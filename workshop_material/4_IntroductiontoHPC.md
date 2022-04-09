# Introduction to HPC

<p style="text-align:left;">
    <b><a class="btn" href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/2_AutomaticVariantC.html" style="background: var(--bs-green);font-weight:bold">&laquo; 2. Automating a Variant Calling Workflow</a></b>
    <span style="float:right;">
     <b><a class="btn" href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/5_working_with_job_scheduler.html" style="background: var(--bs-green);font-weight:bold">4. Working with job Scheduler &raquo;</a></b>
    </span>
</p>

## Outline
* Do not remove this line (it will not be displayed)
{:toc}

## Introduction to HPC

### Defining high-performance computing

The simplest way of defining high-performance computing is by saying that it is the using of high-performance computers (HPC). However, this leads to our next question what is a HPC .

>>A high-performance computer is a network of computers in a cluster that typically share a common purpose and are used to accomplish tasks that might otherwise be too big for any one computer.

<br>
<p>While modern computers can do a lot (and a lot more than their equivalents 10-20 years ago), there are limits to what they can do and the speed at which they are able to do this. One way to overcome these limits is to pool computers together to create a cluster of computers. These pooled resources can then be used to run software that requires more total memory, or need more processors to complete in a reasonable time.</p>

<p>One way to do this is to take a group of computers and link them together via a network switch. Consider a case where you have five 4-core computers. By connecting them together, you could run jobs on 20 cores, which could result in your software running faster.</p>

### HPC architectures

<p>Most HPC systems follow the ideas described above of taking many computers and linking them via network switches. What distinguishes a high-performance computer from the computer clusters described above is:</p>
<br>

* The number of computers/nodes 
* The strength of each individual computer/node 
* The network interconnect – this dictates the communication speed between nodes. The faster this speed is, the more a group of individual nodes will act like a unit.


### NeSI Mahuika Cluster architecture

NeSI Mahuika cluster (CRAY HPE CS400) system consists of a number of different node types. The ones visible to researchers are:

* Login nodes
* Compute nodes
<br>
<p align="center"><img src="nesi_images/hpc_arch_new_fixalignment.png" alt="drawing" width="700"/></p> 
<br>

<br>
<p align="center"><img src="nesi_images/node_overview.png" alt="drawing" width="500"/></p> 
<br>
In reality

<p align="center"><img src="nesi_images/mahuika_maui_real.png" alt="drawing" width="700"/></p>



---

<p style="text-align:left;">
    <b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/3_RNAseq.html">&laquo; 3. RNA-seq Mapping And Count Data Workflow</a></b>
    <span style="float:right;">
     <b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/5_working_with_job_scheduler.html">5. Working with job Scheduler &raquo;</a></b>
    </span>
</p>

<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>
