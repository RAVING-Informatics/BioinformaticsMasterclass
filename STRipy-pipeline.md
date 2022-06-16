# STRipy-pipeline
In this tutorial, we will practice our ability to set-up and run a simple informatics pipeline that runs on CLI.
The pipeline we will use is STRipy, which was created by Andreas Halman with the aim of simplifying the process of genotyping known STR loci. STRipy was initially made as a GUI application, but this fromat has limited throughput capacity, so Andreas developed a pipeline version to query multiple samples and loci at once. 

## Prerequisites
To learn more about STRipy, there is a [website](https://stripy.org/), and an accompanying [paper](https://onlinelibrary.wiley.com/doi/10.1002/humu.24382) for detailed methods. The website contains a curated reference [database](https://stripy.org/) of the known STR expansions, including the gene implicated, repeat region, repeat thresholds, disease association, and the number of repeats in a set of ~2500 unaffected controls.  

## Part 1 - Set up and install STRipy
NOTE: I have already installed stripy on the virtual machine, however I am going to get each of you to do it yourselves in your own environment.

  > **Task 1: Create a subdirectory within the `/data/bioinformatics_masterclass` directory. Keep the name simple e.g. your first name. This directory will act as our own working directory for this tutorial.**
- Ignore this step - it has already been done in the previous tutorial.

```bash
cd /data/bioinformatics_masterclass/
mkdir <enter_name> ; cd <enter_name>
```

Whenever there is a new tool that you want to use, you need to find the repository in which the program code is stored, along with the relevent documentation that describes how and in what environment you can install it. Often, this repository is made using GitHub, and the link to the repository is available in the paper/website associated with the tool/program. Otherwise, you can search for the author's name in GitHub, or run a simple Google search. 

  > **Task 2: Find the repository containing the STRipy-pipeline code and installation instructions**

You will notice that the repository for STRipy-pipeline is actually kept in GitLab, which has functionality and an interface that is brodly the same as GitHub. 

  > **Task 3: Read the README.md file, which gives an overview of the tool and how to set it up**

For now, just read the [Description](https://gitlab.com/andreassh/stripy-pipeline#description) and [Requirements and installation instructions](https://gitlab.com/andreassh/stripy-pipeline#requirements-and-installation-instructions) sections. We can look at pipelne configuration once it is set-up. 

  > **Task 4: Clone the repository**
The first step in installation is to clone the stripy repository to create a local copy of all the files in your own environment. Use the code provided in the `README.md` file for this step. Change into the directory and have a look inside. 

```bash
git clone git@gitlab.com:andreassh/stripy-pipeline.git
cd stripy-pipeline ; ls
```

  > **Task 5: Download dependencies**
Dependencies are all the other programs required for a given program/pipeline to run. You will notice that STRipy has five (5) dependencies listed in the [Requirements and installation instructions](https://gitlab.com/andreassh/stripy-pipeline#requirements-and-installation-instructions) section. The links to each tool's repository are provided. I will go through one by one how to install them. Remember to install each in your own working directory each time. Not in the main `/data` directory. Also don't nest the dependencies within each other. Keep them separate. 

- *Python 3*
There are lots of versions of Python available. Different programs are built on different versions, and for them to run effectively, you need to install the specified verison of Python. STRipy requires Python 3. To install Python, we are going to use Conda. 

1. Create a new conda environment with the required version of python already installed.
```bash
conda create -n stripy_<name> python=3.9.12
```
2. Activate your new conda envionment.
```bash
conda activate stripy_<name>
```
4. Check that the python version is correct and which python is being used (i.e. check that the python executable is in the correct conda envionment). 
```bash
which python ; python --version
```
5. Add the required Python modules for STRipy
```bash
python3 -m pip install -r requirements.txt
```

- *ExpansionHunter*

1. Navigate to the ExpansionHunter releases page (there is a link to this page in the `README.md` file). This page contains precompiled binaries of different ExpansionHunter software versions for different OS environments. The page also has the original source code for each version, if you wish to compile the program yourself. Right click the `ExpansionHunter-v5.0.0-linux_x86_64.tar.gz` file and select `Copy link address`.
2. Use `wget` to copy the file into your workspace. Make sure you are still in your own working directory to avoid rewriting existing directories. 
```bash
wget <insert copied URL>
```
3. Notice the `.tar.gz` extension on the file? This indicates that it has been compressed as a "tar zipped" file. It needs to be decompressed to be used. 
```bash
tar -xvf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
```
4. This will create a directory called `ExpansionHunter-v5.0.0-linux_x86_64/` with a subdirectory `ExpansionHunter/` that contains the program. Have a look inside the directories and subdirectories using `ls`. Also, you no longer need the compressed version, so remove it using `rm`.
```bash
ls ExpansionHunter-v5.0.0-linux_x86_64/ExpansionHunter
rm ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
```
NOTE: The ExpansionHunter executable file is kept in the `bin/` directory. An executable is a file that is used to perform various functions or operations on a computer. They are typically stored in the `bin/` directory by default (but not because they are trash lol).

- REViewer
The installation of REViewer will be identical to that of ExpansionHunter, with one important exception. Can you figure out what it is?

- Samtools
Installation of Samtools is a little different as you would build the package from source. We aren't going to do this step, and instead will use the already installed samtools. But below is the installation code for reference. 
```bash
wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
tar -xvf samtools-1.15.1.tar.bz2 
ls samtools-1.15.1/ ; rm samtools-1.15.1.tar.bz2 
cd samtools-1.15.1/
./configure --prefix=.
make
make install
```

- BWA
Finally, install BWA using instructions on the GitHub page.
```bash
git clone https://github.com/lh3/bwa.git
cd bwa; make
```

  > **Task 6: Add the dependencies to the $PATH**
I am not going to get you to do this in fear of messing up the .bashrc file, and also because there will be too many instances of the same program in the $PATH if we do. But this is a very useful this to do so I will explain what it does but not how to do it. 

Essentially, the $PATH is an environmental variable. Environment variables hold values related to the current environment, like the Operating System or user sessions. $PATH specifies the directories in which executable programs are located on the machine that can be started without knowing and typing the whole path to the file on the command line.

So rather than having to type the full path to the executable file: `/data/samtools-1.13/bin/samtools` you can type `samtools`. This is a great shortcut. 

## Part 2 - Configure STRipy
Navigate to the ["Configure"]([url](https://gitlab.com/andreassh/stripy-pipeline#configuration)) subsection of the STRipy-pipeline page. Have a read of hte diferent pipeline parameters that you can modify and their default values. 
- All of these parameters can be adjusted in using a file called `config.json` which should be in your stripy-pipeline working directory as it was cloned from the git lab repository.

 > **Task 7: Navigate to the config.json file and open it to have a look how it is formatted**
- The first four parameters define where the executable file is stored for each of the dependencies. As mentioned, these would usually be added to the $PATH environmental variable so that you can just type the tool name rather than spelling out the entire path. Because we are not adding the tool to the $PATH, you will need to specify the full path name in the config file. 
- For the remaining parameters, just keep as default. 

## Part 3 - Run STRipy on Test Data
Navigate to the ["Usage"](https://gitlab.com/andreassh/stripy-pipeline#usage) subsection of the STRipy-pipeline page. Have a read.

 > **Task 8: To run STRipy-pipeline we are going to make a run script. Make a new file using `nano` called `run-stripy.sh` within the `stripy-pipeline/` working directory.**
```bash
cd stripy-pipeline #if you're not already there
nano run-stripy.sh
```
This will open a new window in terminal where you can edit the file. Start typing! Use the code in the stripy git lab as a guide (printed below for reference)
```bash
python3 stri.py \
        --genome hg38 \
        --reference reference/hg38.fa \
        --output results/ \
        --locus AFF2,ATXN3,HTT,PHOX2B \
        --input examples/Sample001.bam
```
Some things to keep in mind while you're writing your script:
- Begin the script with `#!/bin/bash` - this is the "she-bang" line or script header; it tells the OS which shell to use. I always use bash by default.  
- Make sure to include the backslash character `\` at the end of each line. This breaks the line of code allowing for the next line continuation. This ensures that all of the parameters specified are linked back to the original command `python3 stri.py`
- Specify the output directory you would like to keep your results in. You will also need to make the directory using `mkdir` command. 
- We use hg38 for our analyses. Hg37 is also an option but the reference .fasta file is in hg38 so stick with that.
- You will notice that you will need a reference genome file in .fasta (.fa). This is already available in the Nimbus instance: `/data/references/Homo_sapiens_assembly38.fasta`
- You will also need some input data. For now, use the example data to check that stripy is working as it should. 
- Specify just a couple of loci as a start. Refer to the [STRipy database]([url](https://stripy.org/database)) for locus options. 
- Use `standard` for the `--analysis` option.

 > **Task 9: Change the permissions on the run script to allow it to run. This is done using the `chmod` command.**
```bash
chmod u+x run-stripy.sh
```
The above command will make the file executable for your user.
You can check the file permissions on a file by running:
```bash
ls -lh run-stripy.sh
```
To understand the output and `chmod u+x`, a good explanation is provided [here]([url](https://kb.iu.edu/d/abdb#:~:text=To%20view%20the%20permissions%20for,in%20a%20directory%20in%20Unix.&text=In%20the%20output%20%20above,a%20file%20or%20a%20directory.)).
 > **Task 10: Run the script!**
```bash
./run-stripy.sh
```
Troubleshoot as need be. 

## Part 3 - Run STRipy on Real Data
All of the data we will be using is in the `/data/bioinformatics_masterclass/bams` folder. These are bam files from WGS of individuals with known, diagnosed STR expansions. See below for a table with the phenotype for each sample. Your task is to use STRipy to determine which locus expansion each individual has and give them a diagnosis!

|     identifier    |     phenotype    |
|-------------------|------------------|
|     S01           |     MND          |
|     S02           |     Ataxia       |
|     S03           |     Ataxia       |
|     S04           |     Ataxia       |
|     S05           |     SBMA         |
|     S06           |     Ataxia       |
|     S07           |     Ataxia       |
|     S08           |     Ataxia       |
|     S09           |     Ataxia       |

Note: Transfer the .HTML files to your local envioronment to view them in a web-browser. 
