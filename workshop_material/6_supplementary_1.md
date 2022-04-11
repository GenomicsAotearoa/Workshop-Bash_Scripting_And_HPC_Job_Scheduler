### Accessing software via modules

On a high-performance computing system, it is quite rare that the software we want to use is available when we log in. It is installed, but we will need to “load” it before it can run.

Before we start using individual software packages, however, we should understand the reasoning behind this approach. The three biggest factors are:

* software incompatibilities
* versioning
* dependencies

One of the workarounds for this issue is Environment modules. A module is a self-contained description of a software package — it contains the settings required to run a software package and, usually, encodes required dependencies on other software packages.

There are a number of different environment module implementations commonly used on HPC systems and the one used in NeSI Mahuika cluster is `Lmod` where the `module` command is used to interact with environment modules.

**Common commands - module**

* View available modules

```bash
#View all modules
$ module avail

# View all modules which match the keyword in their name
$ module avail KEYWORD

# View all modules which match the keyword in their name or description
$ module spider KEYWORD
```
* Load a specific program

    >Note: All modules on NeSI have version and toolchain/environment suffixes. If none is specified, the default version for the tool is loaded. The default version can be seen with the module avail command.

```bash
$ module load MY_APPLICATION
```


* Unload all current modules

```bash
$ module purge
```
>Please **do not** use `$module --force purge`

* Swap a currently loaded module for a different one

```bash
$ module switch CURRENT_MODULE DESIRED_MODULE
```

<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>