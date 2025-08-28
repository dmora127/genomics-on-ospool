
#### For _Simplex_ basecalling - TO BE REMOVED (ALWAYS USE DUPLEX SPLITS)


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
