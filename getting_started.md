# Getting Started

General instructions on the use of these tutorials

## Summary

This document describes how the directory structure should be set up, how commands can be copied and executed from tutorial instructions, how text files can be used and edited, and how files required for tutorials can be downloaded from the tutorial repository.

## Table of contents

* [Directory structure](#directory_structure)
* [Copying and pasting commands](#copy_pasting)
* [Working with text files](#text_files)
* [Downloading tutorial files](#downloading)

<a name="directory_structure"></a>
## Directory structure

Following the tutorials in this collection may be easiest if you follow the instructions closely. We go through each command and explain what it is doing. All required input data and scripts are available on Saga and you should copy them to your working directory on Saga during the course.

The data and scripts for each exercise can be found in the directories

_/cluster/projects/nn9458k/phylogenomics/week1_

_/cluster/projects/nn9458k/phylogenomics/week2_

Additionally, you will find literature to which we refer in the tutorials in the directory

_/cluster/projects/nn9458k/phylogenomics/literature_

In the first exercise, we will ask you to generate your own working directory within the "_phylogenomics_" directory using your name

_/cluster/projects/nn9458k/phylogenomics/YOURNAME_

We will then ask you to generate subfolder for each exercise to keep things orgainzed, such as

_/cluster/projects/nn9458k/phylogenomics/YOURNAME/bayesian\_phylogenetic\_inference_

A link to each tutorial is provided at the [starting page](https://github.com/ForBioPhylogenomics/tutorials/blob/main/README.md).

To work on Saga you will need a console program like the program called Terminal on Mac OS X and Linux or the Windows Console on Windows. On Windows, Windows Console might be sufficient for the course, but a more versatile option might be an emulator like [Cmder](https://cmder.app) (see [Requirements](requirements.md)). The console program will be important to execute commands on Saga and transfer data to your own computer, which will be required occasionally, for example to look at results with a GUI program. You will recognize console commands in the tutorials by monospace font, gray background, and an outline, like for example this command:

		pwd
		
To execute commands like the above, the "Enter" key always needs to be hit after writing the command. The command `pwd` stands for "_print working directory_", which is exactly what it does.

To access Saga, you have to use the following command:

		ssh YOURUSERNAME@saga.sigma2.no
		
You will be asked for your password, which you will need to type to get access.

To upload data to Saga, you can use the following command:

		scp -R YOURFILE YOURUSERNAME@saga.sigma2.no:/cluster/projects/nn9458k/phylogenomics/YOURNAME/PATH_TO_YOUR_FOLDER_ON_SAGA

Downloading from SAGA works very similar. However, we suggest that you first navigate to the directory on your own computer into which you would like to download the files. Then, you can use the following command for the download:

		scp -R YOURUSERNAME@saga.sigma2.no:/cluster/projects/nn9458k/phylogenomics/YOURNAME/PATH_TO_YOUR_FOLDER_ON_SAGA/YOURFILE .
		
Note that instead of a file name in the above two commands (_YOURFILE_), you could also specify the name of a directory, or you could copy multiple files at once by using the asterisk symbol. For example, this command would download all files with the ending _.txt_:

		scp -R YOURUSERNAME@saga.sigma2.no:/cluster/projects/nn9458k/phylogenomics/YOURNAME/PATH_TO_YOUR_FOLDER_ON_SAGA/*.txt .

<a name="copy_pasting"></a>
## Copying and pasting commands

After you've navigated to the tutorial directory on Saga, it should be possible to execute all commands that are mentioned in the tutorial instructions exactly as they are specified. All one-line commands can be copied from the instructions and pasted into the console window, but care should be taken not to include whitespace symbols before the first letter or after the last letter of the command. Sometimes, it is not evident from the text selection that is used for copying, but whitespace nevertheless appears after the last letter copied to the command line (recognizable as a space between the last letter and the cursor); in that case, the whitespace on the command line should be deleted before executing the command.

For multi-line commands, however, it is safer to copy and paste each line individually, and hit the Enter key each time a line has been copied (again, whitespace before and after the first and last letters of the line should not be copied).
	
* Try to copy and paste this command into the console window, and then execute it:

		for n in {1..10}
		do
			echo ${n}
		done

	This should result in the numbers 1 to 10 being written line by line in the console window; if this is not the resulting output, something must have gone wrong with copy-pasting.

<a name="text_files"></a>
## Working with text files

In some cases, text files will need to be written or edited, and for this, a text editor of some sort will be required. For example, instead of copying many lines of commands into the console window one by one, these could also be copied into a text file, the file could be named, e.g., `script.sh`, and then all commands from that file could be executed jointly the command

		bash script.sh
		
To write and edit text files, either a command-line text editor or one with a graphical user interface could be used. There are many suitable options in the latter category. Programs like TextEdit (on Mac OS X) or Notepad (on Windows) could be used, but only if the default settings are changed so that plain text instead of "rich text" is written. Using Microsoft Word would guarantee trouble sooner rather than later. More convenient than any of these, however, are those that include syntax highlighting, such as [BBEdit](https://www.barebones.com/products/textwrangler/) for Mac OS X, [Notepad++](https://notepad-plus-plus.org) for Windows, or [Sublime Text](https://www.sublimetext.com), [Geany](https://www.geany.org), and [Atom](https://atom.io) that run on all three platforms.

If you want to avoid installing any of these programs, you could also use one of the text editors that are usually available within the console, like Vim or Emacs. No installation should be required for these tools, but their use may not be as intuitive as that of text editors with graphical user interfaces. Short tutorials are available online for [Vim](https://www.howtoforge.com/vim-basics) and [Emacs](https://www.digitalocean.com/community/tutorials/how-to-use-the-emacs-editor-in-linux).

* Use the text editor of your choice to write the following four lines:

		for n in {1..10}
		do
			echo ${n}
		done

* Save the text to a new file named `script.sh`, and then place this script in your working directory on Saga.

* Execute the commands in file `script.sh` with this command on Saga:

		bash script.sh
		
	This should again result in the numbers 1 to 10 being written on separate lines.

<a name="downloading"></a>
## Downloading tutorial files

Links to input files and scripts are included in all tutorials and to run a given tutorial, the linked files should be downloaded and placed in your working directory. Most often, these files are hosted in the same GitHub repository, and will not automatically download when the link is clicked (assuming that the tutorial instructions are opened in a web browser).

* To illustrate the point that linked files from the same GitHub repository are not directly downloaded when the link is clicked, click on this link to an alignment file from the tutorial on substitution model selection:
[`locus_0001.fas`](week1_day1_afternoon/data/locus_0001.fas)

	You should see that clicking this link opens another webpage in which the alignment file is merely embedded. The reason for this is simply the way in which GitHub repositories are organized.
	
* To get to the alignment file itself rather than the webpage in which it is embedded, click the button labeled with "Raw" at the top right of the window, just above the embedded file.

	Depending on the file and the web browser used, this should either open the file in the browser or download it directly to your computer, or a dialog window should appear that allows you to select the download location. If the file opened in the browser, use the menu of the browser (e.g. the File menu and the option "Save As...") to download the file and save it in the tutorial directory. If it downloaded directly, it can most likely be found in the Download directory of your computer, and should then be moved into the tutorial directory.
	
	Perhaps a simpler way to download files is to use right-click when clicking on the "Raw" button; this should open a dialog that allows one to select the download location.
	
	Instead of downloading the file through the browser (which would download it to your own computer), you could also just copy its URL, and then download it through the command line, for example with the `wget` command:
	
		wget https://raw.githubusercontent.com/ForBioPhylogenomics/tutorials/main/week1_day1_afternoon/data/locus_0001.fas
	
* Make sure that the file has been placed in the expected directory, either with the file manager program or the command line. To do so on the command line, you could use this command:

		ls
		
	This command lists all files in the current directory, and if that list includes the name of the alignment file (`locus_0001.fas`), the download worked fine. Also make sure that the web browser or your file system has not automatically added another file extension such as `.txt` to the file name. Depending on the settings of your file manager program, these file extensions may not be shown but nevertheless be there, so unless you are sure that your file manager is not hiding any file extensions, it may be safer to check the file name with on the command line.
	
* As a final test that the file has been downloaded correctly, you could also have a look at the content of the file with this command:

		less locus_0001.fas
		
	This should present the file content in the console window. To return to the command line, type the letter `q`.
	
Note that even though downloading of files may seem trivial, past users of this set of tutorials have often encountered problems that were caused by files that were placed in the wrong locations, with the wrong names, or including the HTML code within which they were embedded on the GitHub webpage.
