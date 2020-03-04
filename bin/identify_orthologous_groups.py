#!/usr/bin/env python3

import os
import sys
import functools
import pandas as pd

#note: might want look into collections.defaultDict



def big_scan_gene_data(info, os2_genes):
	"""
	Scans through gene ids in info to match them with their orthologous group ids in the os2_genes 
    file. Reads the whole file into memory before scanning.

    info : A dictionary containing gene ids to be matched with their corresponding orthologous groups.

    os2_genes : A lookup table used to match genes with their corresponding orthologous group ids.

    returns : The updated dictionary, containing the desired gene ids and their corresponding 
              orthologous group ids.
	"""

	df = pd.read_csv(os2_genes, "\t", header=None, names=["orthologous_group_id","gene_id"], index_col=1, dtype={"orthologous_group_id": str, "gene_id": str})
	df.sort_values("gene_id", inplace=True)

	for gene in info:
		try:
			for orthologous_group in df.loc[info[gene]["gene_id"], "orthologous_group_id"]:
				if info[gene].get("orthologous_group_id", None) == None:
					info[gene].update({"orthologous_group_id": [str(orthologous_group)]})
				else:
					info[gene]["orthologous_group_id"].append(str(orthologous_group))
		except:
			pass

	return info



def big_scan_reference_data(info, os2_refs):
	"""
	Scans through reference data for annotations specific to taxonomic ids. Reads the whole file 
	into memory before scanning.

    info : A dictionary containing gene ids, orthologous group ids, orthologous group names, 
    and taxonomic ids.

    os2_refs : A lookup table containing a mapping of taxonomic ids to annotation data.

    returns : The updated dictionary, now containing taxonomic id specific annotations.
	"""

	df = pd.read_csv(os2_refs, "\t", header=None, names=["orthologous_group_id","domain","ipr_num","tax_id"], index_col=3, dtype={"orthologous_group_id": str, "domain": str, "ipr_num": str, "tax_id": str})
	df.sort_values("tax_id", inplace=True)

	for gene in info:
		if info[gene].get("taxonomic_id", None) != None:
			for tax_id in info[gene]["taxonomic_id"]:
				try:
					subset = df.query("tax_id == " + tax_id).values
					for x in range(subset.shape[0]):
						ref_info = subset[x].tolist()
						if info[gene].get("reference_id", None) == None:
							info[gene].update({"reference_id": [str(tax_id)]})
							info[gene].update({"reference_info": [ref_info]})
						else:
							info[gene]["reference_id"].append(str(tax_id))
							info[gene]["reference_info"].append(ref_info)
				except:
					pass

	return info




def big_scan_taxonomic_data(info, os2_tax_id):
	"""
	Scans through orthologous group ids and updates the information dictionary with their 
    orthologous group names and taxonomic ids. Reads the whole file into memory before scanning.

    info : A dictionary containing the genes and their corresponding orthologous group ids.

    os2_tax_id : A lookup table containing a mapping of orthologous group ids and their 
    corresponding taxonomic ids and orthologous group names.

    returns : The updated dictionary, now containing taxonomic ids and orthologous group
    names for the given gene ids and orthologous group ids.
	"""

	df = pd.read_csv(os2_tax_id, "\t", header=None, names=["orthologous_group_id","tax_id", "orthologous_group_name"], index_col=0, dtype={"orthologous_group_id": str, "tax_id": str, "orthologous_group_name": str})
	df.sort_values("orthologous_group_id", inplace=True)

	for gene in info:
		if info[gene].get("orthologous_group_id", None) != None:
			try:
				subset = df.loc[df.index.intersection(info[gene]["orthologous_group_id"])].values
				for x in range(subset.shape[0]):
					tax_id = subset[x][0]
					orthologous_group_name = subset[x][1]
					if info[gene].get("taxonomic_id", None) == None:
						info[gene].update({"taxonomic_id": [str(tax_id)]})
						info[gene].update({"orthologous_group_name": [str(orthologous_group_name)]})
					else:
						info[gene]["taxonomic_id"].append(tax_id)
						info[gene]["orthologous_group_name"].append(str(orthologous_group_name))
			except:
				pass

	return info



	
def mem_size_needed(blast_txt, os2_genes, os2_tax_id, os2_refs):
	"""
	Calculates the estimated maximum amount of memory needed to run this script all the way 
	through.

    blast_txt : A table created through Annotator, contains genes for witch we want to know the 
                groups they are attached to, passed in as a command line argument.

    os2_genes : A lookup table used to match genes with their corresponding orthologous group ids, 
                passed in as a command line argument.

    os2_tax_id : A lookup table containing a mapping of orthologous group ids and their 
                 corresponding taxonomic ids and orthologous group names, passed in as a command 
                 line argument.

    os2_refs : A lookup table containing a mapping of taxonomic ids to annotation data, passed in 
               as a command line argument.

    returns : The estimated max memory needed for this script.
	"""

	max_mem = max([
		os.path.getsize(os2_genes),
		os.path.getsize(os2_tax_id),
		os.path.getsize(os2_refs),
		])

	max_mem += os.path.getsize(blast_txt)

	return max_mem




def output_info(info, output_file):
    """
    Ouputs information gathered about the specified genes inputed from Annotator in a tab dilimited
    format in the order:
    	Gene ID + "\t" + OG ID + "\t" + OG Name + "\t" + Tax ID + "\t" + Annotation Data
    
    info : A dictionary containing gene ids, orthologous group ids, orthologous group names, 
           taxonomic id, and annotation data.

    ouput_file : the desired file to write the information into.

    returns : This method returns nothing.
    """

    output_file.write(
        "Gene_ID"
        + "\t"
        + "OG_ID"
        + "\t"
        + "OG_Name"
        + "\t"
        + "Tax_ID"
        + "\t"
        + "Annotation_Data\n"
        )
	
    for gene in info.keys():
        if info[gene].get("reference_id", None) != None:
            for x in range(len(info[gene]["reference_id"])):

                to_print = [info[gene]["gene_id"] + "\t"]

                for y in range(len(info[gene]["orthologous_group_id"])):
                    if info[gene]["taxonomic_id"][y] == info[gene]["reference_id"][x]:
                        to_print.append(info[gene]["orthologous_group_id"][y] + "\t")
                        to_print.append(info[gene]["orthologous_group_name"][y] + "\t")
                        to_print.append(info[gene]["taxonomic_id"][y] + "\t")
              
                for item in info[gene]["reference_info"][x]:
                    to_print.append( str(item) + "\t")
                
                to_print.append("\n")

                output_file.write(functools.reduce(str.__add__, to_print))




def read_in_gene_ids(info, blast_txt):
	"""
	Scans in the gene ids to look up information about.

    info : An empty dictionary to be populated with gene and orthologous group id information.

    blast_txt : A table created through Annotator, contains genes for witch we want to know the 
                groups they are attached to.

    returns : The updated dictionary, containing the desired gene ids and their corresponding 
              orthologous group ids.
	"""

	for blast_line in blast_txt:
		orthodb_gene_id = blast_line.split("\t")[2]

		if orthodb_gene_id == "Hit_id":
				continue

		info.update({"Gene: " + str(orthodb_gene_id): {"gene_id" : orthodb_gene_id}})
	return info




def small_scan_gene_data(info, os2_genes):
	"""
	Scans through gene ids in info to match them with their orthologous group ids in the os2_genes 
    file.

    info : A dictionary containing gene ids to be matched with their corresponding orthologous groups.

    os2_genes : A lookup table used to match genes with their corresponding orthologous group ids.

    returns : The updated dictionary, containing the desired gene ids and their corresponding 
              orthologous group ids.
	"""

	for gene in info.keys():

		os2_genes.seek(0)
		for gene_line in os2_genes:

			gene_id = gene_line.split("\t")[1].rstrip()
			orthologous_group = gene_line.split("\t")[0]

			if gene_id == info[gene]["gene_id"]:
				if info[gene].get("orthologous_group_id", None) == None:
					info[gene].update({"orthologous_group_id": [orthologous_group]})
				else:
					info[gene]["orthologous_group_id"].append(orthologous_group)
	return info




def small_scan_reference_data(info, os2_refs):
	"""
	Scans through reference data for annotations specific to taxonomic ids.

    info : A dictionary containing gene ids, orthologous group ids, orthologous group names, 
    and taxonomic ids.

    os2_refs : A lookup table containing a mapping of taxonomic ids to annotation data.

    returns : The updated dictionary, now containing taxonomic id specific annotations.
	"""

	for gene in info.keys():
		if info[gene].get("taxonomic_id", None) != None:
			for tax_key in info[gene]["taxonomic_id"]:

				os2_refs.seek(0)
				for ref_line in os2_refs:

					ref_id = ref_line.split("\t")[3].rstrip()
					ref_info = ref_line.split("\t")[0:3]

					if ref_id == tax_key:
						if info[gene].get("reference_id", None) == None:
							info[gene].update({"reference_id": [ref_id]})
							info[gene].update({"reference_info": [ref_info]})
						else:
							info[gene]["reference_id"].append(ref_id)
							info[gene]["reference_info"].append(ref_info)
	return info




def small_scan_taxonomic_data(info, os2_tax_id):
	"""
	Scans through orthologous group ids and updates the information dictionary with their 
    orthologous group names and taxonomic ids.

    info : A dictionary containing the genes and their corresponding orthologous group ids.

    os2_tax_id : A lookup table containing a mapping of orthologous group ids and their 
    corresponding taxonomic ids and orthologous group names.

    returns : The updated dictionary, now containing taxonomic ids and orthologous group
    names for the given gene ids and orthologous group ids.
	"""

	for gene in info.keys():
		if info[gene].get("orthologous_group_id", None) != None:
			for group_id in info[gene]["orthologous_group_id"]:

				os2_tax_id.seek(0)
				for tax_line in os2_tax_id:

					tax_id = tax_line.split("\t")[0]
					tax_key = tax_line.split("\t")[1]
					orthologous_group_name = tax_line.split("\t")[2].rstrip()

					if tax_id == group_id:
						if info[gene].get("taxonomic_id", None) == None:
							info[gene].update({"taxonomic_id": [tax_key]})
							info[gene].update({"orthologous_group_name": [orthologous_group_name]})
						else:
							info[gene]["taxonomic_id"].append(tax_key)
							info[gene]["orthologous_group_name"].append(orthologous_group_name)
	return info

					  


def main():
	"""
	Main takes in arguments from Annotator and assembles annotation information about the desired 
    genes passed in. The function of this program will change based on how much memory the user 
	has available, if there is enough to put the largest table into ram, it will.

    blast_txt : A table created through Annotator, contains genes for witch we want to know the 
                groups they are attached to, passed in as a command line argument.

    os2_genes : A lookup table used to match genes with their corresponding orthologous group ids, 
                passed in as a command line argument.

    os2_tax_id : A lookup table containing a mapping of orthologous group ids and their 
                 corresponding taxonomic ids and orthologous group names, passed in as a command 
                 line argument.

    os2_refs : A lookup table containing a mapping of taxonomic ids to annotation data, passed in 
               as a command line argument.

    returns : This method returns nothing of consequence to the user.
	"""

	info = {}

	mem_needed = mem_size_needed(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

	mem_gigs = (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")) / (1024 ** 3)
	if mem_gigs > mem_needed:

		with open(sys.argv[1], "r") as blast_txt:
			info = read_in_gene_ids(info, blast_txt)

		with open(sys.argv[2], "r") as os2_genes:
			info = big_scan_gene_data(info, os2_genes)

		with open(sys.argv[3], "r") as os2_tax_id:
			info = big_scan_taxonomic_data(info, os2_tax_id)

		with open(sys.argv[4], "r") as os2_refs:
			info = big_scan_reference_data(info, os2_refs)

		with open(sys.argv[5], "w") as output_file:
			output_info(info, output_file)
		
	else:
		with open(sys.argv[1], "r") as blast_txt:
			info = read_in_gene_ids(info, blast_txt)

		with open(sys.argv[2], "r") as os2_genes:
			info = small_scan_gene_data(info, os2_genes)

		with open(sys.argv[3], "r") as os2_tax_id:
			info = small_scan_taxonomic_data(info, os2_tax_id)

		with open(sys.argv[4], "r") as os2_refs:
			info = small_scan_reference_data(info, os2_refs)

		with open(sys.argv[5], "w") as output_file:
			output_info(info, output_file)
	

if __name__ == "__main__":
	main()

