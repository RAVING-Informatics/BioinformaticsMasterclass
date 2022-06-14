# Working with files and environments

When starting a new project, it's a good idea to make a new environemnt. This is because the tools we use on the command line have a lot of dependencies. When you install two tools, they might have conflicting dependencies. For example, STRipy needs python3.9.12, but SpliceAI needs python3.8.

To solve this issue, we install them in environments which can't see each other. This also means it's harder to accidently break something

In this workshop, we will create an environment using anaconda, a package and environment manager. Conda also makes it very simple to install bioinformatics tools.

## Creating an environment

We have already installed conda.

Let's make an environment with some tools for viewing and manipulating VCF files.

1. See usage instructions

   ```bash
   conda --help
   ```

   How do we create a new environment?

   Try

   ```bash
   conda create --help
   ```
   
2. See a list of all of the environments

   ```bash
   conda env list
   ```

2. Create new environemnt

   ```bash
   conda create -n bcftools --channel bioconda bcftools_<your_name>
   ```
   This installs bcftools, a tool for manipulating VCF files.
   `-n` is the flag for the environment name, and `--channel` (or `-c`) tells conda to install bcftools from the bioconda channel.
   When prompted, press `y` and enter.

3. Activate your environment

   ```bash
   conda activate bcftools_<your_name>
   ```
   You should see `(bcftools_<your_name>)` at the start of your bash prompt line. This is how you know which environment you are in.

## Working with VCF files

We're now going to use bcftools to work with a joint-called exome VCF

4. View a compressed VCF

   Most VCFs are compressed as they take up a lot of space and are too large to view. They are compressed using gzip ot bgzip, similar to compressed zip folders on a  computer. They generally have the suffix `.gz` and are accompanied by an index file with the suffix `.tbi` or `.csi`. \
   We can view a file on bash simply using `less <my_file>` but this can't view a compressed file. (Press `q` to quit).
   We could unzip the file using `gunzip <my_file>.vcf.gz` to view it but they are generally too big.
   
   ```bash
   cd <your_working_directory>         #Remember you can use tab to autocomplete file, folder, or script names
   cp ../D22_1234.hg38.vcf.gz ../D22_1234.hg38.vcf.gz.tbi .
   ```
   
   Print the header
   ```bash
   bcftools view -h D22_1234.hg38.vcf.gz
   ```
   
   Why don't we just view everything? 
   ```bash
   bcftools view D22_1234.hg38.vcf.gz | wc -l
   ```
   We 'piped' `|` the full output of bcftools view (the whole VCF) to the bash word count `wc` function in line count `-l` mode.
     
5. Splitting a VCF into individual samples

   Often some tools require only one individual per VCF file (eg. linkage analysis), while others can use a cohort VCF.
   We can split a VCF using bcftools, and we can specify if we want the output VCFs to be compressed or not.
   
   There are multiple steps to do this, but they can all be combined into one.
   1. `bcftools query -l <my_file>` is used to list the individual samples within the VCF
   2. `bcftools view -s <my_individual>` is used to extract the specified individual.
   We could run it individually for each, but joint VCFs can have hundreds or thousands of individuals.
   We can combine all of these steps together into a simple bash script.
   
   We can view and edit the bash script with `nano <my_script.sh>`, a simple text editor.
   ```bash
   cp ../split_vcf.sh
   nano split_vcf.sh
   ```
   You will see something like this:

   ```bash
   #! /bin/bash                              #A 'shebang' is at the top of almost all bash scripts and it tells the shell where to find the interpretter for running your script

   for file in old_joint.vcf.gz              #We start a for loop, where we tell it to find just our file of interest (similar to how `ls <my_file>` would return just that file)
   do
           for sample in `bcftools query -l $file`       #We start a loop over the output of ``bcftools query -l $file``, where each loop takes one `sample` (one individual) and passes it to the next line:
           do
                   bcftools view -Oz -s $sample -o $sample.vcf $file       #We use `bcftools view -s $sample` to save our VCF using the inidividual's name defined in the VCF header.
           done                              #Close the first loop
   done                                      #Close the second loop
   ```
   
   
   I last used this on a different file, so we need to replace the filename to that of our joint called VCF.
   Where it says `old_joint.vcf.gz` change it to your VCF's name
   1. Press **Ctrl + X** to exit 
   2. Press **Y** to save 
   3. Press **Enter** to keep the same name 
   
   Now we can run our script:
   ```bash
   bash split_vcf.sh
   ```

6. Challenge 🤓 Use just one line to:
   1. Get the number of output files
   2. Unzip all of the files
   > Hint: everything you will need ot use has already been mentioned, and you don't need to use bcftools
