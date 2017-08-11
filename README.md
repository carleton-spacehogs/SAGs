# SAGs

## Welcome!

Welcome to the SAGs repo README. This document outlines the contents of folders, files, and documents included in the SAG repo.

As of latest commit, SAGs contains the following folders and files:

1. [bubble_plots](#a): a collection of bubble plots by SAG
2. [dN_dS_scatter_images](#b): scatterplots of dN/dS by ITEP cluster presence/absence
3. [KO_SNV_misc](#c): miscellaneous files for Seq Object data structures: awaiting sorting or deletion.
4. [misc_images](#d): miscellaneous/one-time images
5. [misc_trash](#e): miscellaneous formatting files awaiting deletion.
6. [pa_files](#f): text files related to ITEP/clusterDbanaylsis presence/absence tables.
7. [python_files](#g): scripts, pipelines, custom and classes written in Python 3. Awaiting further organization and documentation.
8. [SAG_data_files](#h): .gff, .ko, .fa assembly files, .tsv contig names (from anvi-script), and .names_map files for each SAG. This is all the data that needs to be integrated for non-ITEP/non-PAML analysis.
9. [variability.zip](#i): zipped anvi'o outputs for SNV and SAAV variability profiles for all 5 SAGs.
10. [papers](#j): the papers I've been reading and pre-existing info on the research

## Folder Overviews
### bubble_plots <a name="a"></a>
### dN_dS_scatter_images <a name="b"></a>
### KO_SNV_misc <a name="c"></a>
### misc_images <a name="d"></a>
### misc_trash <a name="e"></a>
### pa_files <a name="f"></a>
### python_files <a name="g"></a>

####Dependencies
These files are written in Python 3. The following modules are used: \

* pandas
* numpy
* matplotlib (pyplot)
* sklearn (commented out)
* scipy
* pickle

**Note:** The majority of the python modules and functions in this directory are intended for single-purpose use on specific file
types. 
1. SeqObj folder: Contains code for custom Python *SNV*, *SAAV*, *Contig*, *ORF*, objects. Designed to store 
all of the information in the SAG_data_files folder in Python objects to easily sort, store, and parse this information.
The Python object modules are:
* *SNV*, 
* *SAAV* 
* *Contig*
* *ORF*
View the module files for a list of data attributes stored in these objects. \
An array of functions for reading the input files and constructing lists of Sequence objects exist: \
* make_contig_ORF_and_SNV_lists.py
* make_SNVs_and_Contigs.py
* makeSequenceObjects.py
* mergeORFpaData.py
* organizeKO.py
* populateContigLists.py
These functions must be called in order to properly initialize the Sequence objects.


### SAG_data_files <a name="h"></a>
### variability.zip <a name="i"></a>
### papers <a name="j"></a>


### Contributors
Michael Hoffert - Undergrad, hoffertm@carleton.edu \
Rika Anderson - PI, Space Hogs