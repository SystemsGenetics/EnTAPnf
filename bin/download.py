#!/usr/bin/env python3
"""
This is the download script.
"""
import json
import os


#
# dictionary Contains all download tasks where each key is a single download
# task containing a list of shell command to run in order. The key is used as
# the directory name that is created and used for the task.
#
DOWNLOAD_TASKS = {}


#
# io.TextIOWrapper Output log file used to generated log messages whenever a
# shell command of a task completes successfully or fails.
#
logfile = open("log.txt","w")


#
# dictionary Contains the progress of each task where the keys are the same as
# the DOWNLOAD TASKS dictionary but its values is just an integer for the next
# shell command not yet done for that task.
#
sessionState = {}




def defaultState():
    """
    Getter function.

    Returns
    -------
    ret0 : dictionary
           A new default session state where no task is complete.
    """
    ret = {}
    for key in DOWNLOAD_TASKS:
        ret[key] = 0
    return ret




def initializeDownloadTasks():
    """
    Initializes the DOWNLOAD TASKS global dictionary.
    """
    global DOWNLOAD_TASKS
    DOWNLOAD_TASKS = {
        "interproscan": [
            "wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz"
            ,"tar -zxvf interproscan-5.36-75.0-64-bit.tar.gz"
        ]
        ,"interproscan/interproscan-5.36-75.0/data": [
            "wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-14.1.tar.gz"
            ,"tar -zxvf panther-data-14.1.tar.gz"
        ]
        ,"nr": [
            "wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
            ,"gunzip nr"
            ,"diamond makedb --threads 4 --in nr -d nr"
        ]
        ,"refseq": [
            "wget -r -A '*.protein.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/"
            ,"gunzip ./ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.gz"
            ,"cat ./ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.faa > refseq_plant.protein.faa"
            ,"rm -rf ./ftp.ncbi.nlm.nih.gov/"
            ,"diamond makedb --threads 4 --in refseq_plant.protein.faa -d refseq_plant.protein"
        ]
        ,"orthodb": [
            "wget https://v101.orthodb.org/download/odb10v1_level2species.tab.gz"
            ,"gunzip odb10v1_level2species.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz"
            ,"gunzip odb10v1_OG2genes.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_OGs.tab.gz"
            ,"gunzip odb10v1_OGs.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz"
            ,"gunzip odb10v1_all_og_fasta.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_OG_xrefs.tab.gz"
            ,"gunzip odb10v1_OG_xrefs.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_species.tab.gz"
            ,"gunzip odb10v1_species.tab.gz"
            ,"wget https://v101.orthodb.org/download/odb10v1_gene_xrefs.tab.gz"
            ,"gunzip odb10v1_gene_xrefs.tab.gz"
            ,"/Annotater/bin/index_orthodb.py ."
            ,"diamond makedb --threads 4 --in odb10v1_all_og_fasta.tab -d odb10v1_all_og"
        ]
        ,"string": [
            "wget https://stringdb-static.org/download/protein.sequences.v11.0.fa.gz"
            ,"gunzip protein.sequences.v11.0.fa.gz"
            ,"wget https://stringdb-static.org/download/protein.links.full.v11.0.txt.gz"
            ,"gunzip protein.links.full.v11.0.txt.gz"
            ,"wget https://stringdb-static.org/download/protein.info.v11.0.txt.gz"
            ,"gunzip protein.info.v11.0.txt.gz"
            ,"/Annotater/bin/index_string.py --links protein.links.full.v11.0.txt --info protein.info.v11.0.txt --out protein"
            ,"diamond makedb --threads 4 --in protein.sequences.v11.0.fa -d protein.sequences.v11.0"
        ]
        ,"uniprot_sprot": [
            "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
            ,"gunzip uniprot_sprot.fasta.gz"
            ,"wget ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"
            ,"diamond makedb --threads 4 --in uniprot_sprot.fasta -d uniprot_sprot"
        ]
        ,"uniprot_trembl": [
            "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
            ,"gunzip uniprot_trembl.fasta.gz"
            ,"diamond makedb --threads 4 --in uniprot_trembl.fasta -d uniprot_trembl"
        ]
    }




def loadSession():
    """
    Loads the global session state from the session JSON file if it exists, else
    it loads a new default session state.
    """
    global sessionState
    if os.path.isfile("session.json"):
        with open("session.json","r") as ifile:
            sessionState = json.loads(ifile.read())
    else:
        sessionState = defaultState()




def processTask(
    dirName
    ):
    """
    Processes the task with the given key in the global DOWNLOAD TASKS, using
    the global session state to determine where to start in the task. Any task
    that is complete is ignored and this does nothing.

    Parameters
    ----------
    dirName : string
              Key, or directory name, of the task that is processed.
    """
    state = sessionState[dirName]
    tasks = DOWNLOAD_TASKS[dirName]
    if state<len(tasks):
        cwd = os.getcwd()
        if not os.path.isdir(dirName):
            os.makedirs(dirName)
        os.chdir(dirName)
        while state<len(tasks):
            if os.system(tasks[state]):
                logfile.write("FAILURE: task %i of %i for %s: %s\n"%(state+1,len(tasks),dirName,tasks[state]))
                logfile.flush()
                os.chdir(cwd)
                return False
            logfile.write("Successfully completed task %i of %i for %s\n"%(state+1,len(tasks),dirName))
            logfile.flush()
            state += 1
            sessionState[dirName] = state
        os.chdir(cwd)
    return True




def saveSession():
    """
    Saves the global session state to the session JSON file.
    """
    with open("session.json","w") as ofile:
        ofile.write(json.dumps(sessionState,indent=4)+"\n")








### __SCRIPT__ ###
initializeDownloadTasks()
loadSession()
try:
    for dirName in sessionState:
        processTask(dirName)
finally:
    saveSession()
