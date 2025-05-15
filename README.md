# Long-Read Genomics on the OSPool

This tutorial will walk you through a complete long-read sequencing analysis workflow using Oxford Nanopore data on the OSPool high-throughput computing ecosystem. You'll learn how to:

* Basecall raw Nanopore reads using the latest GPU-accelerated Dorado basecaller
* Map your reads to a reference genome using Minimap2
* Call structural variants using Sniffles2
* Breakdown massive bioinformatics workflows into many independent smaller tasks
* Submit hundreds to thousands of jobs with a few simple commands
* Use the Open Science Data Federation (OSDF) to manage file transfer during job submission

All of these steps are distributed across hundreds (or thousands!) of jobs using the HTCondor workload manager and Apptainer containers to run your software reliably and reproducibly at scale. The tutorial is built around realistic genomics use cases and emphasizes performance, reproducibility, and portability. You'll work with real data and see how high-throughput computing (HTC) can accelerate your genomics workflows.

>[!NOTE]
>If you're brand new to running jobs on the OSPool, we recommend completing the HTCondor ["Hello World"](https://portal.osg-htc.org/documentation/htc_workloads/workload_planning/htcondor_job_submission/) exercise before diving into this tutorial.

**Let’s get started!**

Jump to...
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
  * [Structural Variant Calling using Sniffles2](#structural-variant-calling-using-sniffles2)
      + [Submitting our Sniffles2 SV jobs to the OSPool](#submitting-our-sniffles2-sv-jobs-to-the-ospool)
      + [Running Minimap to Map Reads to the Reference Genome](#running-minimap-to-map-reads-to-the-reference-genome)
   * [Next Steps](#next-steps)
      + [Software](#software)
      + [Data](#data)
      + [GPUs](#gpus)
   * [Getting Help](#getting-help)

<!-- TOC end -->

## Tutorial Setup

### Assumptions

This tutorial assumes that you:

* Have basic command-line experience (e.g., navigating directories, using bash, editing text files).
* Have a working OSPool account and can log into an Access Point (e.g., ap40.uw.osg-htc.org).
* Are familiar with HTCondor job submission, including writing simple .sub files and tracking job status with condor_q.
* Understand the general workflow of long-read sequencing analysis: basecalling → mapping → variant calling.
* Have access to a machine with a GPU-enabled execution environment (provided automatically via the OSPool).
* Have sufficient disk quota and file permissions in your OSPool home and OSDF directories.

>[!TIP]
>You do not need to be a genomics expert to follow this tutorial. The commands and scripts are designed to be beginner-friendly and self-contained, while still reflecting real-world research workflows.

### Materials

To obtain a copy of the files used in this tutorial, you can

* Clone the repository, with 
  
  ```
  git clone https://github.com/osg-htc/tutorial-long-read-genomics
  ```

  or the equivalent for your device

* Download the zip file of the materials: 
  [download here](https://github.com/osg-htc/tutorial-long-read-genomics/archive/refs/heads/main.zip)

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
> [!CAUTION]
> These commands must be ran **every time you login and to build a new container**. Building Apptainer containers on the Access Point without first running the commands below places excessive strain on shared storage resources and **violates OSPool usage policies**. Failure to follow these steps may result in restricted access to the system.

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

## Structural Variant Calling using Sniffles2

### Submitting our Sniffles2 SV jobs to the OSPool
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

Now that you've completed the long-read genomics tutorial on the OSPool, you're ready to adapt these workflows for your own data and research questions. Here are some suggestions for what you can do next:

🧬 Apply the Workflow to Your Own Data
* Replace the tutorial datasets with your own POD5 files and reference genome.
* Modify the basecalling, mapping, and variant calling submit files to fit your data size, read type (e.g., simplex vs. duplex), and resource needs.

🧰 Customize or Extend the Workflow
* Incorporate quality control steps (e.g., filtering or read statistics) using FastQC.
* Use other mappers or variant callers, such as ngmlr, pbsv, or cuteSV.
* Add downstream tools for annotation, comparison, or visualization (e.g., IGV, bedtools, SURVIVOR).

📦 Create Your Own Containers
* Extend the Apptainer containers used here with additional tools, reference data, or dependencies.
* For help with this, see our [Containers Guide](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers/).

🚀 Run Larger Analyses
* Submit thousands of basecalling or alignment jobs across the OSPool.
* Explore data staging best practices using the OSDF for large-scale genomics workflows.
* Consider using workflow managers (e.g., [DAGman](https://portal.osg-htc.org/documentation/htc_workloads/automated_workflows/dagman-workflows/) or [Pegasus](https://portal.osg-htc.org/documentation/htc_workloads/automated_workflows/tutorial-pegasus/)) with HTCondor.

🧑‍💻 Get Help or Collaborate
* Reach out to [support@osg-htc.org](mailto:support@osg-htc.org) for one-on-one help with scaling your research.
* Attend office hours or training sessions—see the [OSPool Help Page](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) for details.

### Software

In this tutorial, we created several *starter* apptainer containers, including tools like: Dorado, SAMtools, Minimap, and Sniffles2. These containers can serve as a *jumping-off* for you if you need to install additional software for your workflows. 

Our recommendation for most users is to use "Apptainer" containers for deploying their software.
For instructions on how to build an Apptainer container, see our guide [Using Apptainer/Singularity Containers](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers-singularity/).
If you are familiar with Docker, or want to learn how to use Docker, see our guide [Using Docker Containers](https://portal.osg-htc.org/documentation/htc_workloads/using_software/containers-docker/).

This information can also be found in our guide [Using Software on the Open Science Pool](https://portal.osg-htc.org/documentation/htc_workloads/using_software/software-overview/).

### Data

The ecosystem for moving data to, from, and within the HTC system can be complex, especially if trying to work with large data (> gigabytes).
For guides on how data movement works on the HTC system, see our [Data Staging and Transfer to Jobs](https://portal.osg-htc.org/documentation/htc_workloads/managing_data/overview/) guides.

### GPUs

The OSPool has GPU nodes available for common use, like the ones used in this tutorial. If you would like to learn more about our GPU capacity, please visit our [GPU Guide on the OSPool Documentation Portal](https://portal.osg-htc.org/documentation/htc_workloads/specific_resource/gpu-jobs/).

## Getting Help

The OSPool Research Computing Facilitators are here to help researchers using the OSPool for their research. We provide a broad swath of research facilitation services, including:

* **Web guides**: [OSPool Guides](https://portal.osg-htc.org/documentation/) - instructions and how-tos for using the OSPool and OSDF.
* **Email support**: get help within 1-2 business days by emailing [support@osg-htc.org](mailto:support@osg-htc.org).
* **Virtual office hours**: live discussions with facilitators - see the [Email, Office Hours, and 1-1 Meetings](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) page for current schedule.
* **One-on-one meetings**: dedicated meetings to help new users, groups get started on the system; email [support@osg-htc.org](mailto:support@osg-htc.org) to request a meeting.

This information, and more, is provided in our [Get Help](https://portal.osg-htc.org/documentation/support_and_training/support/getting-help-from-RCFs/) page.