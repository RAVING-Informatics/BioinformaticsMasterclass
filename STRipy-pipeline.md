# STRipy-pipeline
In this tutorial, we will practice our ability to set-up and run a simple informatics pipeline that runs on CLI.
The pipeline we will use is STRipy, which was created by Andreas Halman with the aim of simplifying the process of genotyping known STR loci. STRipy was initially made as a GUI application, but this fromat has limited throughput capacity, so Andreas developed a pipeline version to query multiple samples and loci at once. 

To learn more about STRipy, there is a [website]([url](https://stripy.org/)), and an accompanying [paper]([url](https://onlinelibrary.wiley.com/doi/10.1002/humu.24382)) for detailed methods. The website contains a curated reference database of the known STR expansions, including the gene implicated, repeat region, repeat thresholds, disease association, and the number of repeats in a set of ~2500 unaffected controls.  

## Part 1 - Set up and install STRipy
NOTE: I have already installed stripy on the virtual machine, however I am going to get each of you to do it yourselves in your own environment.

  > **Task 1: Create a subdirectory within the `/data/BioinformaticsMasterclass` directory. Keep the name simple e.g. your first name. This directory will act as our own working directory for this tutorial.**

```bash
cd /data/BioinformaticsMasterclass/
mkdir <enter_name> ; cd <enter_name>
```

Whenever there is a new tool that you want to use, you need to find the repository in which the program code is stored, along with the relevent documentation that describes how and in what environment you can install it. Often, this repository is made using GitHub, and the link to the repository is available in the paper/website associated with the tool/program. Otherwise, you can search for the author's name in GitHub, or run a simple Google search. 

  > **Task 2: Find the repository containing the STRipy-pipeline code and installation instructions**

You will notice that the repository for STRipy-pipeline is actually kept in GitLab, which has functionality and an interface that is brodly the same as GitHub. 

  > **Task 3: Read the README.md file, which gives an overview of the tool and how to set it up**

For now, just read the [Description]([url](https://gitlab.com/andreassh/stripy-pipeline#description)) and [Requirements and installation instructions]([url](https://gitlab.com/andreassh/stripy-pipeline#requirements-and-installation-instructions)) sections. We can look at pipelne configuration once it is set-up. 

  > **Task 4: Clone the repository**
The first step in installation is to clone the stripy repository to create a local copy of all the files in your own environment. Use the code provided in the `README.md` file for this step.

```bash
git clone git@gitlab.com:andreassh/stripy-pipeline.git
```

  > **Task 5: Download dependencies**
Dependencies are all the other programs required for a given program/pipeline to run. You will notice that STRipy has five (5) dependencies listed in the [Requirements and installation instructions]([url](https://gitlab.com/andreassh/stripy-pipeline#requirements-and-installation-instructions)) section. The links to each tool's repository are provided. I will go through one by one how to install them. Remember to install each in your own working directory each time. Not in the main `/data` directory. Also don't nest the dependencies within each other. Keep them separate. 

- Python 3
There are lots of versions of Python available. Different programs are built on different versons, and for them to run effectively, you need to install the specified verison of Python. STRipy requires Python 3. To install Python, we are going to use Conda. 

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
- ExpansionHunter
7. Navigate to the ExpansionHunter releases page (there is a link to this page in the `README.md` file). This page contains precompiled binaries of different ExpansionHunter software versions for different OS environments. The page also has the original source code for each version, if you wish to compile the program yourself. Right click the `ExpansionHunter-v5.0.0-linux_x86_64.tar.gz` file and select `Copy link address`.
8. Use `wget` to copy the file into your workspace. Make sure you are still in your own working directory to avoid rewriting existing directories. 
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
The ExpansionHunter executable file is kept in the `bin/` directory. An executable is a file that is used to perform various functions or operations on a computer. They are typically stored in the `bin/` directory by default (but not because they are trash lol).

- REViewer
The installation of REViewer will be identical to that of ExpansionHunter, with one important exception. Can you figure out what it is?

- Samtools
Installation of Samtools is a little different as you will build the package from source. 
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
