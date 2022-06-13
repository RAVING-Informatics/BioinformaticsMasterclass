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
conda create -n vcftools --channel bioconda bcftools
```
   This installs VCFtools, a tool for manipulating VCF files.
   `-n` is the flag for the environment name, and `--channel` (or `-c`) tells conda to install bcftools from the bioconda channel.
   When prompted, press `y` and enter.

3. 
