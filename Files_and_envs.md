# Working with files and environments

When starting a new project, it's a good idea to make a new environemnt. This is because the tools we use on the command line have a lot of dependencies. When you install two tools, they might have conflicting dependencies. For example, STRipy needs python3.9, but SpliceAI needs python3.8.

To solve this issue, we install them in environments which can't see each other. This also means it's harder to accidently break something :)

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

2. Create new environemnt

   ```bash
   conda create -n bcftools --channel bioconda bcftools
   ```
   This installs VCFtools, a tool for manipulating VCF files.
   `-n` is the flag for the environment name, and `--channel` (or `-c`) tells conda to install bcftools from the bioconda channel.
   When prompted, press `y` and enter.

3. Activate your environment

   ```bash
   conda activate bcftools
   ```
   You should see `(bcftools)` at the start of your bash prompt line. This is how you know which environment you are in.
   
4. View a compressed VCF

   Most VCFs are compressed as they take up a lot of space and are too large to view. They are compressed using gzip ot bgzip, similar to compressed zip folders on a  computer. They generally have the suffix `.gz`
   We can view a file on bash simply using `less <my_file>` but this can't view a compressed file. (Press `q` to quit).
   We could unzip the file using `gunzip <my_file>.vcf.gz` but this is slow and unessecary. So we can use bcftools.
   
   ```bash
   cd <your_working_directory>
   cp ../D22_1234.hg38.vcf.gz .
   ```
   
   Print the header
   ```bash
   bcftools view -h D22_1234.hg38.vcf.gz
   ```
   
   Why don't we view everything? 
   ```bash
   bcftools view D22_1234.hg38.vcf.gz | wc -l
   ```
   We 'piped' `|` the full output 
