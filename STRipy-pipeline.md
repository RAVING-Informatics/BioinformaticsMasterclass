# STRipy-pipeline
In this tutorial, we will practice our ability to set-up and run a simple informatics pipeline that runs on CLI.
The pipeline we will use is STRipy, which was created by Andreas Halman with the aim of simplifying the process of genotyping known STR loci. STRipy was initially made as a GUI application, but this fromat has limited throughput capacity, so Andreas developed a pipeline version to query multiple samples and loci at once. 

To learn more about STRipy, there is a [website]([url](https://stripy.org/)), and an accompanying [paper]([url](https://onlinelibrary.wiley.com/doi/10.1002/humu.24382)) for detailed methods. The website contains a curated reference database of the known STR expansions, including the gene implicated, repeat region, repeat thresholds, disease association, and the number of repeats in a set of ~2500 unaffected controls.  

## Part 1 - Set up and install STRipy
Whenever there is a new tool that you want to use, you need to find the repository in which the program code is stored, along with the relevent documentation that describes how and in what environment you can install it. Often, this repository is made using GitHub, and the link to the repository is available in the paper/website associated with the tool/program. Otherwise, you can search for the author's name in GitHub, or run a simple Google search. 

  > **Task 1: Find the repository containing the STRipy-pipeline code and installation instructions**

You will notice that the repository for STRipy-pipeline is actually kept in GitLab, which has functionality and an interface that is brodly the same as GitHub. 

  > **Task 2: Read the README.md file, which gives an overview of the tool and how to set it up**

For now, just read the [Description](Description) and [Requirements and installation instructions]([url](https://gitlab.com/andreassh/stripy-pipeline#requirements-and-installation-instructions)) sections. We can look at pipelne configuration once it is set-up. 

I have already installed stripy on the virtual machine, however I am going to get each of you to do it yourselves in your own environment.

  > **Task 3: Create a directory within the /data directory that will act as a working directory for this tutorial**

