# Bioinformatics Masterclass: Introduction

## **Logging in**

We're logging in to our virtual machine (Nimbus) hosted by Pawsey @ Curtin. 
1. Open MobaXterm
2. Click 'Session'
3. Click 'SSH'
4. Under Remote host, put in the server address `146.118.68.201`
5. Click 'Advanced SSH settings'
6. Tick 'Use private key'
7. Click the blue folder next to it and navigate to the private key file that was sent to you ('mykey2.pem')
8. After getting prompted for 'login as', type 'ubuntu'

You should then see something like
`ubuntu@ravinginformatics2:~$`

After the `$`, we can type our commands.

## **Introducing the shell**

Let's start with some very simple navigation :)

N.B. Unix commands are case sensitive (always lowercase)

1. Find your current working directory

   ```bash
   pwd
   ```

   What do you see? you should get something like `/home/ubuntu`. This is the default landing sirectoy after logging in and is where we can store files and folders. 'Ubuntu' is our username (and the name of the operating system we are using).

   What if I don't know how to use a command?
   `man pwd` and you will see a user manual. Press `q` to quit. If `man` doesn't work, try `pwd --help`.

2. List the files in the current working directory

    ```bash
    ls
    ```

    What do you see?
    Let's move to our data storage directory

3. Change directory

     ```bash
     cd /data
     ```

     Now you can type in `ls` again to see what we have in `/data` and what we've been working on :) 
     Folders are highlighted in blue, and files in white .Let's move to todays work folder. 

     ```bash
     cd bioinformatics_masterclass
     ```

     What would happen if you tried to do `cd /bioinformatics_masterclass`? What's the differrence?

     Now we're in todays working folder. Since we are all working on the same project today, it's a good idea to create our own folder's so we don't ruin each others work.

4. Create a folder

     ```bash
     mkdir <yourname>
     ```

     Replace `<yourname>` with yourname, obviously.
     `cd` into your new directory

     We know how to move forward into a new directory, but how do we move back?
     You could always go back to `/data` with `cd /data` then into your subdirectories, but there's a much easier way

5. Change directory back

     ```bash
     cd ..
     ```

     You can also view (`ls`) the previous directory if you don't want to change.
     To go back more than one directory, simply run `cd ../..` (and so on until you are at the root folder).

     Now we have set up todays workshop folder, let's copy in the files we need.

6. Copy files and folders

     The formula for copying is simple, `cp` followed by file(s), then a destination

     ```bash
     cp <some_file> <some_destination>
     ```

     Let's copy these instructions from `bioinformatics_masterclass/` to your folder (`/data/bioinformatics_masterclass/<yourname>/`)

     ```bash
     cp intro_to_bash.md <yourname>
     ```

     Hint: typing out `into_to_bash.md` takes WAY too long, so we can use a shortcut. Once you've typed the first few letters of the file/folder/function, you can autocomplete it by pressing 'TAB'.

     What if we wanted to copy 100 BAMs, or 500 VCFs? We can use wildcards to save us from using `cp` 600 times.

     `cp *.bam <destination>`

     This will copy anything with the suffix `.bam` to the destination. We could also copy all of the files associated with an individual, eg. `cp D22-1234* <destination>`. `*` is the wildcard for unlimited number of characters. For a single character, we can use `?`

     How do we move a file to a new destination?

7. Moving files and folders

     This works in the same way as `cp`.

     ```bash
     mv <original_name> <destination>
     ```

     Sometimes, when copying or downloading files from the internet, there could be two files with the same name. So how do we rename a file or folder?

8. Renaming files and folders

     ```bash
     mv <original_name> <new_name>
     ```

     But it's not a good idea to use wildcards for this. Overwirtten files are not recoverable.
