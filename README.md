# Long-Read Genomics on the OSPool

This is a start-to-finish tutorial for how to deploy and run R calculations on CHTC's High Throughput Computing (HTC) system.
The goal of this tutorial is to provide a step-by-step example of how to go from running R calculations on your computer using RStudio,
to running *many* such calculations on a High Throughput Computing system.

Using data from the [NOAA Global Historical Climatology Network](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ncdc:C00861/html),
this tutorial generates histograms showing the distribution of daily high and low temperatures across the four meteorological seasons.

Jump to...

- [Tutorial Setup](#tutorial-setup)
- [Example Calculation Using RStudio](#example-calculation-using-rstudio)
- [Transitioning from RStudio](#transitioning-from-rstudio)
- [Logging into CHTC](#logging-into-chtc)
- [Run Example Calculation as a Single Job](#run-example-calculation-as-a-single-job)
- [Run Example Calculation as Multiple Jobs](#run-example-calculation-as-multiple-jobs)
- [Next Steps](#next-steps)
- [Getting Help](#getting-help)
- [Appendix](#appendix-preparing-to-transition-an-r-project-from-your-computer-to-chtc)

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Long-Read Genomics on the OSPool](#long-read-genomics-on-the-ospool)
   * [Tutorial Setup](#tutorial-setup)
      + [Assumptions](#assumptions)
      + [Materials](#materials)
   * [Basecalling Oxford Nanopore long reads using Dorado](#basecalling-oxford-nanopore-long-reads-using-dorado)
      + [Setting up our software environment](#setting-up-our-software-environment)
      + [Data Wrangling and Splitting Reads](#data-wrangling-and-splitting-reads)
         - [For _Simplex_ basecalling](#for-simplex-basecalling)
         - [For _Duplex_ basecalling](#for-duplex-basecalling)
      + [Submitting your basecalling jobs](#submitting-your-basecalling-jobs)
   * [Mapping Sequencing Reads to Genome](#mapping-sequencing-reads-to-genome)
      + [Data Wrangling and Splitting Reads](#data-wrangling-and-splitting-reads-1)
      + [Running Minimap to Map Reads to the Reference Genome](#running-minimap-to-map-reads-to-the-reference-genome)
   * [Next Steps](#next-steps)
      + [Software](#software)
      + [Data](#data)
      + [GPUs](#gpus)
   * [Getting Help](#getting-help)
   * [Appendix: Preparing to transition an R project from your computer to CHTC](#appendix-preparing-to-transition-an-r-project-from-your-computer-to-chtc)
      + [Locate your R scripts](#locate-your-r-scripts)
      + [Locate your other input files](#locate-your-other-input-files)
      + [Check your R scripts for "absolute paths"](#check-your-r-scripts-for-absolute-paths)
      + [Find your version of R](#find-your-version-of-r)
         - [Console](#console)
         - [Packages pane](#packages-pane)
         - [Command](#command)
      + [Identify your R packages](#identify-your-r-packages)
      + [About paths](#about-paths)
      + [About versions](#about-versions)

<!-- TOC end -->

## Tutorial Setup

### Assumptions

This tutorial assumes that you have been using R via the RStudio program installed on your computer.
To participate in the hands-on portion, you'll need an active CHTC account.
You can request such an account here: [go.wisc.edu/chtc-account](https://go.wisc.edu/chtc-account).
To make it easier to copy and paste commands, we recommend that you have this tutorial's GitHub page open in your browser.

> [!TIP]
> It is recommended, though not required, that you complete the "Hello World" guide [Practice: Submit HTC Jobs using HTCondor](https://chtc.cs.wisc.edu/uw-research-computing/htcondor-job-submission) before starting this tutorial.

### Materials

To obtain a copy of the files used in this tutorial, you can

* Clone the repository, with 
  
  ```
  git clone https://github.com/osg-htc/tutorial-long-read-genomics
  ```

  or the equivalent for your device

* Download the zip file of the materials: 
  [download here](https://github.com/CHTC/tutorial-rstudio-to-chtc/archive/refs/heads/main.zip)

We recommend that you create a new R project named "rstudio-to-chtc" and download/copy the files into that directory.

In the RStudio toolbar,

## Basecalling Oxford Nanopore long reads using Dorado

### Setting up our software environment
Before we can begin basecalling our reads, we need to setup our software environment to run Dorado. We are going to setup
our environment using an Apptainer container. 

1. First, let's login to our OSPool Account

    ```
    ssh user.name@ap40.uw.osg-htc.org
    ```

2. Run the following commands on the Access Point (AP) to set our temporary apptainer build directory to `/home/tmp/`.
Once you've built your container, you can delete the contents of this directory to reduce quota usage on `/home`. 

    ```
    mkdir -p $HOME/tmp
    export TMPDIR=$HOME/tmp
    export APPTAINER_TMPDIR=$HOME/tmp
    export APPTAINER_CACHEDIR=$HOME/tmp
    ```

3. We now need to write up a definition file for singularity to build our Dorado container. Copy and paste this block of
text to a new file titled `dorado.def`. You can open up a text editor, such as `vim` or `nano` using a command like: `vim dorado.def`.

    ```
    Bootstrap: docker
    From: nvidia/cuda:12.8.1-cudnn-devel-ubuntu24.04
    
    %post
        DEBIAN_FRONTEND=noninteractive
    
        # system packages
        apt-get update -y
        apt-get install -y \
                build-essential \
                wget
    
        apt install -y python3-pip
       
        # install Dorado and POD5
        cd /opt/
        wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.5-linux-x64.tar.gz
        tar -zxvf dorado-0.9.5-linux-x64.tar.gz
        rm dorado-0.9.5-linux-x64.tar.gz
        
        # install POD5 using pip
        pip install pod5 --break-system-packages
    
    %environment
    
        # set up environment for when using the container
        # add Dorado to $PATH variable for ease of use
        export PATH="/opt/dorado-0.9.5-linux-x64/bin/:$PATH"
    ```

    This definition file uses the Nvidia CUDA Ubuntu 24.04 base image and installs necessary packages to run Dorado and POD5
    on our data.


4. Build your apptainer container on the Access Point (AP) by running the following command:
    ```
    apptainer build dorado.sif dorado.def
   ```
   
5. Move your finalized container image `dorado.sif` to your `OSDF` directory
    
    ```
   mv dorado.sif /ospool/ap40/data/<user.name>/
   ```
   
6. We need to repeat steps 3-5 for samtools. Since samtools has an existing container on Dockerhub, we can use it as the base image to build our container. 
   
    ```
   apptainer build /ospool/ap40/data/<user.name>/samtools.sif docker://biocontainers/samtools:v1.9-4-deb_cv1
   ```

### Data Wrangling and Splitting Reads

Oxford Nanopore sequencing runs general yield POD5 files. Each POD5 file is generated about once an hour throughout the
duration of the sequencing run. This output format does not scale very well, as data output usually plateus after 24-48hrs.
This would mean that POD5 files that are generated from earlier in the sequencing run, will be larger in size compared to
files later in the run. Additionally, this division of data does not allow for _Duplex_ read basecalling. As a result prior
to running Dorado, we must first reorganize the data contained within all the POD5 files. 

---

#### For _Simplex_ basecalling
When basecalling our sequencing data using simplex basecalling mode on Dorado we can subdivide our POD5 files into smaller individual subsets. This subdivision of our files enables us to take advantage of the OSPool's High Throughput Computing (HTC) principles, significantly decreasing the time-to-results for our basecalling. We will use the `POD5` package installed in our `dorado.sif` container—if you need to generate the `dorado.sif` apptainer image, refer to [Setting up our software environment](#Setting-up-our-software-environment). 

1. Create a csv that maps reads in your `pod5_dir` to subset files.
    ```
    pod5 view <path_to_pod5_dir> --include "read_id" | awk 'NR==1 {print "read_id subset_id"; next} {print $0, int((NR-1)/1000)}' > read_subsets.csv
   ```
   _This will generate a CSV table mapping each read_id to a subset_file for basecalling._


2. Using the `read_subsets.csv` mapping file, subset your POD5 reads to the output directory `split_pod5_subsets`
    ```
   pod5 subset <path_to_pod5_dir> --summary read_subsets.csv --columns subset_id --output split_pod5_subsets
   ```
   
3. Create a list of POD5 files to iterate through while basecalling

    ```
   ls split_pod5_subsets > /home/<user.name>/genomics_tutorial/pod5_files
   ```
   
    If you `head` this new file you should see an output similar to this:

    ```
    [user.name@ap40 user.name]$ head /home/<user.name>/genomics_tutorial/pod5_files
    subset_id-0.pod5
    subset_id-100.pod5
    subset_id-101.pod5
    subset_id-102.pod5
    subset_id-103.pod5
    subset_id-104.pod5
    subset_id-105.pod5
    subset_id-106.pod5
    subset_id-107.pod5
    subset_id-108.pod5
    [user.name@ap40 user.name]$ 
   ```

----

#### For _Duplex_ basecalling
When basecalling our sequencing data using simplex basecalling mode on Dorado we can subdivide our POD5 files into smaller individual subsets. This subdivision of our files enables us to take advantage of the OSPool's High Throughput Computing (HTC) principles, significantly decreasing the time-to-results for our basecalling. We will use the `POD5` package installed in our `dorado.sif` container—if you need to generate the `dorado.sif` apptainer image, refer to [Setting up our software environment](#Setting-up-our-software-environment). 

1. Create a csv that maps reads in your `pod5_dir` to subset files.
    ```
    pod5 view <path_to_pod5_dir> --include "read_id, channel" --output summary.tsv
   ```
   _This will generate a TSV table mapping each read_id to each channel for basecalling._

2. Using the `read_subsets.csv` mapping file, subset your POD5 reads to the output directory `split_pod5_subsets`
    ```
   pod5 subset <path_to_pod5_dir> --summary summary.tsv --columns channel --output split_pod5_subsets
   ```
   
3. Create a list of POD5 files to iterate through while basecalling

    ```
   ls split_pod5_subsets > /home/<user.name>/genomics_tutorial/pod5_files
   ```
   
    If you `head` this new file you should see an output similar to this:

    ```
    [user.name@ap40 user.name]$ head /home/<user.name>/genomics_tutorial/pod5_files
    channel-100.pod5
    channel-101.pod5
    channel-102.pod5
    channel-103.pod5
    channel-104.pod5
    channel-105.pod5
    channel-106.pod5
    channel-107.pod5
    channel-108.pod5
    channel-109.pod5
    [user.name@ap40 user.name]$ 
   ```

### Submitting your basecalling jobs


1. Create your Dorado simplex basecalling executable - `/home/<user.name>/genomics_tutorial/executables/basecalling_step2_simplex_reads.sh`

    ```
    #!/bin/bash 
    # Run Dorado on the EP for each POD5 file (non-resumeable)
    #dorado ${dorado_arg_string} ${input_pod5_file} > ${output_bam_file}
    echo "Running dorado with: $1"
    
    # untar your Dorado basecalling models
    tar -xvzf models.tar.gz
    rm models.tar.gz
    echo "completed tar unzip"
    
    args="$1"
    eval "dorado $args"
   ```

2. Create your submit file for Dorado simplex basecalling - `/home/<user.name>/genomics_tutorial/submit_files/basecalling_step2_simplex_reads.sub`

    ```
    +SingularityImage      = "osdf:///ospool/ap40/data/<user.name>/dorado.sif"

    executable		       = ../executables/basecalling_step2_simplex_reads.sh
    arguments		       = "'basecaller --batchsize 16 hac@v5.0.0 --models-directory ./models/ $(POD5_input_file) > $(POD5_input_file).bam'"
    
    transfer_input_files   = osdf:///ospool/ap40/data/<user.name>/split_by_channels/$(POD5_input_file), osdf:///ospool/ap40/data/<user.name>/models.tar.gz

    transfer_output_files  = ./$(POD5_input_file).bam
    output_destination	   = osdf:///ospool/ap40/data/<user.name>/basecalledBAMs/
    
    output                 = ./basecalling_step2/logs/$(POD5_input_file)_$(Process)_basecalling_step2.out
    error                  = ./basecalling_step2/logs/$(POD5_input_file)_$(Process)_basecalling_step2.err
    log                    = ./basecalling_step2/logs/$(POD5_input_file)_$(Process)_basecalling_step2.log
    
    request_cpus           = 1
    request_disk           = 8 GB
    request_memory         = 24 GB 
    request_gpus		   = 1
    
    queue POD5_input_file from /home/<user.name>/genomics_tutorial/pod5_files
   ```

This submit file will read the contents of `/home/<user.name>/genomics_tutorial/pod5_files`, iterate through each line, and assign the value of each line to the variable `$POD5_input_file`. This allows to us programmatically submit _N_ jobs, where _N_ is equal to the number of POD5 subset file we created previously. Each job will have its corresponding POD5 input subset (`subset_id-105.pod5`) and `models.tar.gz` files transferred to the Execution Point (EP). Additionally, we will transfer and start our `dorado.sif` apptainer container image using the `+SingularityImage` attribute on our submit file. 

The submit file will instruct the EP to run our executable `basecalling_step2_simplex_reads.sh` and pass the arguments found in the `arguments` attribute. The `arguments` attribute allows us to customize the parameters passed to _Dorado_ directly on our submit file, without having to edit our executable. 

> [!NOTE]  
> The example submit script above is running the hac@v5.0.0 model for simplex basecalling. You can change this to `duplex sup --models-directory ./models/ $(POD5_input_file) > $(POD5_input_file).bam'`. For additional usage information, refer to the [Dorado User Documentation](https://github.com/nanoporetech/dorado).

3. Submit your set of basecalling jobs

    ```
   condor_submit basecalling_step2_simplex_reads.sub
   ```
   
    You can track the progress of your jobs with the `condor_q` command
    
> [!TIP] 
> You may experience some `held` jobs due to a variety of resource allocation overruns, including using more memory or CPUs than request. We recommend you use the following commands to edit those held jobs and resubmit them. 
>
> ```
> [user.name@ap40 user.name]$ condor_q <cluster.id> -hold
> 12345678.123      user.name       5/6  19:47 Excessive CPU usage. Job used 3 CPUs, while request_cpus=1. Please verify that the code is configured to use a limited number of cpus/threads, and matches request_cpus.
> [user.name@ap40 user.name]$ condor_qedit 12345678.123 requestCpus=4
> [user.name@ap40 user.name]$ condor_release 12345678.123
> Job 12345678.123 released
> [user.name@ap40 user.name]$ 
> ```



## Mapping Sequencing Reads to Genome

### Data Wrangling and Splitting Reads

To get ready for our mapping step, we need to prepare our freshly basecalled reads. You should have a directory with several BAM files, these BAM files need to be sorted and 

1. Generate a list of the BAM files Dorado created during basecalling. Save it as `listOfBasecallingBAMs` in your OSDF directory. 

    ```
   ls /ospool/ap40/data/<user.name>/basecalledBAMs/ > /ospool/ap40/data/<user.name>/listOfBasecallingBAMs
   ```

### Running Minimap to Map Reads to the Reference Genome

1. Indexing our reference genome - Generating `Celegans_ref.mmi`

   1.  Create `minimap2_index.sh` using either `vim` or `nano`
        ```
       #!/bin/bash
       minimap2 -x map-ont -d Celegans_ref.mmi Celegans_ref.fa
       ```
   2. Create `minimap2_index.sub` using either `vim` or `nano`
        ```
        +SingularityImage      = "osdf:///ospool/ap40/data/<user.name>/genomics_tutorial/minimap2.sif"
    
        executable		       = ../executables/minimap2_index.sh
        
        transfer_input_files   = osdf:///ospool/ap40/data/<user.name>/genomics_tutorial/Celegans_ref.fa
    
        transfer_output_files  = ./Celegans_ref.mmi
        output_destination	   = osdf:///ospool/ap40/data/<user.name>/genomics_tutorial/
        
        output                 = ./minimap2/logs/$(Cluster)_$(Process)_indexing_step1.out
        error                  = ./minimap2/logs/$(Cluster)_$(Process)_indexing_step1.err
        log                    = ./minimap2/logs/$(Cluster)_$(Process)_indexing_step1.log
        
        request_cpus           = 4
        request_disk           = 10 GB
        request_memory         = 24 GB 
        
        queue 1
       ```
   3. Submit your `minimap2_index.sub` job to the OSPool
       ```
      condor_submit minimap2_index.sub
      ```
> [!WARNING]  
> Index will take a few minutes to complete, **do not proceed until your indexing job is completed**

2. Map our basecalled reads to the reference _C. elegans_ indexed genome - `Celegans_ref.mmi`
    
   1.  Create `minimap2_mapping.sh` using either `vim` or `nano`
       ```
       #!/bin/bash
       # Use minimap2 to map the basecalled reads to the reference genome
        ./minimap2 -ax map-ont Celegans_ref.mmi "$1" > "mapped_${1}_reads_to_genome.sam"
       
       # Use samtools to sort our mapped reads BAM, required for downstream analysis
       samtools sort "mapped_${1}_reads_to_genome.sam" -o "mapped_${1}_reads_to_genome_sam_sorted.bam"
       ```
   2. Create `minimap2_mapping.sub` using either `vim` or `nano`
       ```
        +SingularityImage      = "osdf:///ospool/ap40/data/<user.name>/minimap2.sif"
    
        executable		       = ../executables/minimap2_index.sh
        arguments              = $(BAM_File)
        transfer_input_files   = osdf:///ospool/ap40/data/<user.name>/Celegans_ref.mmi, osdf:///ospool/ap40/data/<user.name>/basecalledBAMs/$(BAM_FILE)
    
        transfer_output_files  = ./mapped_$(BAM_FILE)_reads_to_genome_sam_sorted.bam
        transfer_output_remaps = "mapped_$(BAM_FILE)_reads_to_genome_sam_sorted.bam=/home/<user.name>/genomics_tutorial/MappedBAMs/mapped_$(BAM_FILE)_reads_to_genome_sam_sorted.bam
        
        output                 = ./minimap2/logs/$(Cluster)_$(Process)_mapping_$(BAM_FILE)_step2.out
        error                  = ./minimap2/logs/$(Cluster)_$(Process)_mapping_$(BAM_FILE)_step2.err
        log                    = ./minimap2/logs/$(Cluster)_$(Process)_mapping_$(BAM_FILE)_step2.log
        
        request_cpus           = 2
        request_disk           = 5 GB
        request_memory         = 10 GB 
        
        queue BAM_File from /home/<user.name>/genomics_tutorial/listOfBasecallingBAMs
       ```
    
        In this step, we **are not** transferring our outputs using the OSDF. The mapped/sorted BAM files are intermediate temporary files in our analysis and do not benefit from the aggressive caching of the OSDF. By default, HTCondor will transfer outputs to the directory where we submitted our job from. Since we want to transfer the sorted mapped BAMs to a specific directory, we can use the `transfer_output_remaps` attribute on our submission script. The syntax of this attribute is:
   
        ```transfer_output_remaps = "<file_on_execution_point>=<desired_path_to_file_on_access_point>``` 
    
   3. Submit your cluster of minimap2 jobs to the OSPool
   
      ```
      condor_submit minimap2_mapping.sub
      ```

## Calling Structural Variants using Sniffles2

### Data Wrangling and Splitting Reads

To get ready for our variant calling step, we need to prepare our freshly mapped BAM files. You should have a directory with several BAM files, these BAM files need to be merged before we can begin. We're going to use our `samtools` container to accomplish this. 

1. Generate a list of the BAM files minimap2/samtools-sort created during the read mapping step. Save it as `listOfBasecallingBAMs` in your OSDF directory. 

    ```
   ls /home/<user.name>/genomics_tutorial/MappedBAMs/ > /home/<user.name>/genomics_tutorial/listOfMappedBAMs
   ```
   
2. Create a tarball of the files in MappedBAMs and save it to `/home/<user.name>/genom/MappedBAMs.tar.gz`

    ```
   cd ~/genomics_tutorial/MappedBAMs/
   tar -czf ../MappedBAMs.tar.gz .
   ```
   
2. Create a `bamMerge.sh` Shell script using `vim` or `nano`
    
    ```
   #!/bin/bash
   
   # Use the tar command to uncompress your MappedBAMs tarball
   tar -xvzf /MappedBAMs.tar.gz
   
   # Use samtools merge method to merge all the mapped and sorted BAM files from the minimap step
   samtools merge -b listOfMappedBAMs -o merged_sorted_mapped_reads_to_Celegans.bam
   ```

3. Create a `bamMerge.sub` HTCondor submission script using `vim` or `nano`

    ```
    +SingularityImage      = "osdf:///ospool/ap40/data/<user.name>/minimap2.sif"
        
    executable		       = ../executables/bamMerge.sh
    
    transfer_input_files   = /home/<user.name>/genomics_tutorial/MappedBAMs.tar.gz, /home/<user.name>/genomics_tutorial/listOfMappedBAMs
    
    transfer_output_files  = ./merged_sorted_mapped_reads_to_Celegans.bam
    output_destination	   = osdf:///ospool/ap40/data/<user.name>/genomics_tutorial/
    
    output                 = ./sniffles/logs/$(Cluster)_$(Process)_merging_step1.out
    error                  = ./sniffles/logs/$(Cluster)_$(Process)_merging_step1.err
    log                    = ./sniffles/logs/$(Cluster)_$(Process)_merging_step1.log
    
    request_cpus           = 8
    request_disk           = 20 GB
    request_memory         = 24 GB 
    
    queue 1
   ```
> [!WARNING]  
> Index will take a few minutes to complete, **do not proceed until your indexing job is completed**

### Structural Variant Calling using Sniffles2

To get ready for our variant calling step, we need to prepare our freshly mapped BAM files. You should have a directory with several BAM files, these BAM files need to be merged before we can begin. We're going to use our `samtools` container to accomplish this.
   
1. Create a `sniffles_sv_calling.sh` Shell script using `vim` or `nano`
    
    ```
   #!/bin/bash
   # Use samtools merge method to merge all the mapped and sorted BAM files from the minimap step
   sniffles -i merged_sorted_mapped_reads_to_Celegans.bam -v merged_sorted_mapped_reads_to_Celegans.vcf
   python3 -m sniffles2_plot -i merged_sorted_mapped_reads_to_Celegans.vcf -o merged_sorted_mapped_reads_to_Celegans_sniffles_plots
   tar -czf ../MappedBAMs.tar.gz .
   ```

2. Create a `bamMerge.sub` HTCondor submission script using `vim` or `nano`

    ```
    +SingularityImage      = "osdf:///ospool/ap40/data/<user.name>/minimap2.sif"
        
    executable		       = ../executables/sniffles_sv_calling.sh
    
    transfer_input_files   = /home/<user.name>/genomics_tutorial/MappedBAMs.tar.gz, /home/<user.name>/genomics_tutorial/listOfMappedBAMs
    
    transfer_output_files  = ./merged_sorted_mapped_reads_to_Celegans.bam, ./merged_sorted_mapped_reads_to_Celegans_sniffles_plots
    output_destination	   = osdf:///ospool/ap40/data/<user.name>/genomics_tutorial/
    
    output                 = ./sniffles/logs/$(Cluster)_$(Process)_svCalling_step2.out
    error                  = ./sniffles/logs/$(Cluster)_$(Process)_svCalling_step2.err
    log                    = ./sniffles/logs/$(Cluster)_$(Process)_svCalling_step2.log
    
    request_cpus           = 12
    request_disk           = 30 GB
    request_memory         = 30 GB 
    
    queue 1

## Next Steps

Now that you've finished this tutorial, you are ready to start transitioning your own R project to be run on the HTC system.
But unless your R project is fairly simple, there are a few more things you'll need to work on to get up and running.

For a full walk-through of how to get started on the HTC system, see our guide [Roadmap to getting started](https://chtc.cs.wisc.edu/uw-research-computing/htc-roadmap).

### Software

This tutorial used a pre-existing container that came with R 4.4.2 and `tidyverse` packages already installed.
If that is all you need, then you're in luck!
Just use the same `container_image` line in your submit file.

If you're like most users, however, then you have additional R packages that you want to use in your scripts. 
To make those packages available for use in your HTC job, we recommend that you build your own container.
While that may sound like a daunting task, we have a lot of documentation and examples to help you get started, 
and the faciliation team is happy to help with any questions or issues.

Our recommendation for most users is to use "Apptainer" containers for deploying their software.
For instructions on how to build an Apptainer container, see our guide [Use Apptainer Containers](https://chtc.cs.wisc.edu/uw-research-computing/apptainer-htc).
If you are familiar with Docker, or want to learn how to use Docker, see our guide [Running HTC Jobs Using Docker Containers](https://chtc.cs.wisc.edu/uw-research-computing/docker-jobs.html).

For examples of containers that you can use or modify, see the [R section of our Recipes GitHub repository](https://github.com/CHTC/recipes/tree/main/software/R/).

This information can also be found in our guide [Overview: How to Use Software](https://chtc.cs.wisc.edu/uw-research-computing/software-overview-htc).

### Data

The ecosystem for moving data to, from, and within the HTC system can be complex, especially if trying to work with large data (> gigabytes).
For guides on how data movement works on the HTC system, see the ["Manage data" section of our HTC guides page](https://chtc.cs.wisc.edu/uw-research-computing/htc/guides.html#manage-data).

### GPUs

If your R project is capable of using GPUs, and you would like to use the GPUs available on the HTC system, see our guide [Use GPUs](https://chtc.cs.wisc.edu/uw-research-computing/gpu-jobs).

## Getting Help

CHTC employs a team of Research Computing Facilitators to help researchers use CHTC computing for their research. 

* **Web guides**: [HTC Computing Guides](https://chtc.cs.wisc.edu/uw-research-computing/htc/guides) - instructions and how-tos for using the HTC system.
* **Email support**: get help within 1-2 business days by emailing [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu).
* **Virtual office hours**: live discussions with facilitators - see the "Get Help" page for current schedule.
* **One-on-one meetings**: dedicated meetings to help new users, groups get started on the system; email [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu) to request a meeting.

This information, and more, is provided in our [Get Help](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) page.

## Appendix: Preparing to transition an R project from your computer to CHTC

### Locate your R scripts

Identify the R scripts that you use to run your calculation. 
Typically you'll have one main R script that is the entry point to your program, and for simple programs this will be the only script.
You can use the "Files" pane to navigate the files in your R project. 

In this tutorial, main script was `example.R`. 
But we also need the script `my_functions.R`, since it is loaded by `example.R`.

For your project, you may have other scripts. 
If you are not sure which or if any of the scripts are needed, take a look at your main R script and see if it references any of the other scripts.
Ideally all of the scripts you use in your calculation are in the same folder (or a subfolder thereof).
If not, you should consider reorganizing your scripts into the project directory. 

### Locate your other input files

Identify the input files besides your R scripts that your calculation needs to function.
To start with, consider what is needed to run a single example calculation. 

In this tutorial, we needed the dataset `.csv` files as input for calculations,
and the files were located in the same directory as the R scripts.

If your input files are not in the same directory as your R scripts, 
you may want to consider consolidating them into the project directory, 
at least for one example calculation.

### Check your R scripts for "absolute paths"

If your R script(s) references or loads other files, or writes outputs to file, you should check if they are using "absolute paths". 
If so, you'll want to rewrite your program to use a "relative" path.
(This is another reason you'll want to consolidate your files into the project directory.)

You will likely need to test that your program still functions as expected.

> [!TIP]
>For more information about "absolute" and "relative" paths, see the note below ([About paths](#about-paths)).

### Find your version of R

There are several ways of finding the version of R that you are using in your project.
Use one or more of them to identify the version, which will be in the pattern `X.Y.Z`.

In the examples below, the version number is `4.4.2`. 
Make sure you note your specific version number.

To minimize the chance of discrepancies, you'll want to use the same version to run your calculations on the HTC system.

#### Console

When you open the console, the very first line contains the version of R, which looks like this:

```
R version 4.4.2 (2024-10-31 ucrt) -- "Pile of Leaves"
```

#### Packages pane

In a box on the right side should be a "Packages" tab that you can click on to open the Packages pane. 
This pane lists packages that are installed (checked box) or that are available to be installed (unchecked box) in your R environment.

Scroll down to the "System Library" section and look for the "base" package, and note the the number in its "Version" column.
This corresponds to the version of R you are using in your environment.

#### Command

You can programmatically identify the version of R that you are using by entering the following command in the R console:

```
R.version.string
```

This will print something like the following:

```
[1] "R version 4.4.2 (2024-10-31 ucrt)"
```

This command can be used wherever you are using R, which makes it useful in scenarios that don't involve RStudio.

### Identify your R packages

Identify the R packages that your project uses, so that later you can reproduce the environment on CHTC.

To start, make a list of the packages that you load in your R scripts, which is generally done using `library('<package_name>')` commands. 
Then, look in the "Packages" pane to identify the corresponding versions of the packages. 
Usually the package names alone is enough, but sometimes the versions of the packages can matter as well ([About versions](#about-versions)). 

> If you'd rather not do this manually, you can install and use a package called `renv` to not only automatically detect the packages you are using,
> but to also create files that can be used to replicate the environment automatically when building a container.
> For more information, see the `renv` recipe in the Recipes repository: [https://github.com/CHTC/recipes/tree/main/software/R/renv](https://github.com/CHTC/recipes/tree/main/software/R/renv).


### About paths

An "absolute path" is used to reference the location of a file in relation to the "root" directory of your computer. 
This is fine when your program is running on your computer, but can break the program if you try to run it on a different computer whose files are organized differently from yours.

A "relative path" is used to reference the location of a file in relation to where the current script is running.
This is useful when you need to run your program on different computers. 

The absolute path to the dataset file on a Windows machine may look like `C:/Users/bbadger/Documents/REPONAME/madison.csv`, 
while on a Mac machine the path may look like `/Users/bbadger/Documents/REPONAME/madison.csv`.

A relative path starts from the current working directory, and defines the location of the file in relation to that.
Such a path may look like `./madison.csv` or `../data/madison.csv`. 
Here, the `.` represents the current directory, while `..` represents the parent directory. 
You can chain together several `..` to go several directories upwards in the file system.

Consider for example the following folder structure:

```
project/
├── data/
│   ├── 2023/
│   │   └── raw_data.csv
│   └── 2024/
│       └── raw_data.csv
└── scripts/
    └── v1/
        └── program.R
```

The script `program.R` can reference the 2024 `raw_data.csv` file using this relative path: `../../data/2024/raw_data.csv`.

### About versions

Most software uses the `Major.Minor.Patch` versioning syntax. 
 
* *Major version number* - A change in this number signals major changes in the software, and commands that worked in the previous version may not work in the new version.
* *Minor version number* - A change in this number signals additional features or enhancements. Commands in previous versions should work fine in later versions, though there may be superficial changes. 
* *Patch version number* - A change in this number signals that bugs have been fixed. There should be no change, superficially or functionally, other than those resulting from correcting the bugs.

**Does the version number matter?**

To a certain extent, yes. 
If your code was written for Major version X, there's no guarantee it will function for a different Major version, so you should continue to use Major version X.
Code written for Minor version Y should function for Minor versions >= Y, but there may be superficial changes you might want to avoid, so it's up to you whether or not to be consistent.
You should always use the latest Patch Version; if there is a discrepancy in your results between two Patch versions, that is (hopefully) because a bug that affected the results has been fixed.
(It is also possible another bug has been introduced - either way, you should investigate the nature of the bug fixes.)
