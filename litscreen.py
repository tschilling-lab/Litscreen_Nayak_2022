import requests
import sys
import pandas as pd
import bs4
import time
import numpy as np
from tqdm import tqdm
import os
import argparse

#HOW TO RUN THIS PROGRAM.
#RUN THE PYTHON FILE ON COMMAND LINE WITH THE FIRST ARGUMENT BEING THE FILE NAME OF YOUR EXCEL FILE CONTAINING YOUR INPUT GENE LIST IN ENSEMBL ID FORMAT
#THE SECOND ARGUMENT SHOULD BE THE SEARCH TERM YOU WISH TO USE WITH EACH GENE AGAINST THE PUBMED DATABASE
#THE THIRD ARGUMENT SHOULD BE THE DESIRED NAME OF YOUR OUTPUT EXCEL FILE e.g. "output-excel-file.xlsx" OR WHATEVER YOU WISH TO CALL IT

#After you have run the python file, the program will ask for the path to the ensembl_gene_lists folder and the orthologous_gene_lists folders. These are required to find the orthologous genes and translate ensembl gene ids to conventional gene names
#Lastly, the program will ask for the path to your input gene list excel file. Please enter the path to the FOLDER CONTAINING YOUR INPUT FILE. This will be added to the working path
#After that, the file should run, and your output excel file should be either in the same path as your input file, or whatever your current path is. i think.


gene_list_input = input("Please enter path to the ensembl_gene_lists folder: ")
orth_list_input = input("Please enter path to the orthologous_gene_lists folder: ")
input_list_input = input("Please enter path to your gene list input excel file folder: ")

gene_list_path = gene_list_input
orth_list_path = orth_list_input
input_file_path = input_list_input


human_gl = os.path.join(gene_list_path, "human_genes.txt")
zebrafish_gl = os.path.join(gene_list_path, "zebrafish_genes.txt")
mouse_gl = os.path.join(gene_list_path, "mouse_genes.txt")
rat_gl = os.path.join(gene_list_path, "rat_genes.txt")

zeb2hum = os.path.join(orth_list_path, "zebrafish2human.txt")
zeb2mou = os.path.join(orth_list_path, "zebrafish2mouse.txt")
zeb2rat = os.path.join(orth_list_path, "zebrafish2rat.txt")
hum2zeb = os.path.join(orth_list_path, "human2zebrafish.txt")
hum2mou = os.path.join(orth_list_path, "human2mouse.txt")
hum2rat = os.path.join(orth_list_path, "human2rat.txt")
mou2hum = os.path.join(orth_list_path, "mouse2human.txt")
mou2zeb = os.path.join(orth_list_path, "mouse2zebrafish.txt")
mou2rat = os.path.join(orth_list_path, "mouse2rat.txt")
rat2hum = os.path.join(orth_list_path, "rat2human.txt")
rat2zeb = os.path.join(orth_list_path, "rat2zebrafish.txt")
rat2mou = os.path.join(orth_list_path, "rat2mouse.txt")

human_genelist = pd.read_csv(human_gl)
zebrafish_genelist = pd.read_csv(zebrafish_gl)
mouse_genelist = pd.read_csv(mouse_gl)
rat_genelist = pd.read_csv(rat_gl)

zebrafish2human = pd.read_csv(zeb2hum)
zebrafish2mouse = pd.read_csv(zeb2mou)
zebrafish2rat = pd.read_csv(zeb2rat)
human2zebrafish = pd.read_csv(hum2zeb)
human2mouse = pd.read_csv(hum2mou)
human2rat = pd.read_csv(hum2rat)
mouse2human = pd.read_csv(mou2hum)
mouse2rat = pd.read_csv(mou2rat)
mouse2zebrafish = pd.read_csv(mou2zeb)
rat2human = pd.read_csv(rat2hum)
rat2zebrafish = pd.read_csv(rat2zeb)
rat2mouse = pd.read_csv(rat2mou)

startTime = time.time()
# IMPORTANT change this next line to the name of your ENSEMBL gene ID list. this is your input
orth_df = pd.read_excel(os.path.join(input_file_path, sys.argv[1]))


# these four functions are to search the dataframes for homologous genes for a particular gene input for each different species
# I created a local gene ortho/homology finder for each species for human, mouse, rat, and zebrafish so far
# each function takes as input a string ensembl ID and returns 3 lists, one list for the homologs for each species for your input gene
# These will be called in the local multispecies homolog finder created later to change the overall input dataframe into one containing all the homologs for all genes in the list

def local_homfinder_human(ens_string):
    hl_human2rat = []
    hl_human2mouse = []
    hl_human2zebrafish = []

    hom_indexlist1 = list(np.where(human2rat["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_human2rat.append(human2rat.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(human2mouse["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_human2mouse.append(human2mouse.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(human2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_human2zebrafish.append(human2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_human2rat, hl_human2mouse, hl_human2zebrafish


def local_homfinder_rat(ens_string):
    hl_rat2human = []
    hl_rat2mouse = []
    hl_rat2zebrafish = []

    hom_indexlist1 = list(np.where(rat2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_rat2human.append(rat2human.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(rat2mouse["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_rat2mouse.append(rat2mouse.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(rat2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_rat2zebrafish.append(rat2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_rat2human, hl_rat2mouse, hl_rat2zebrafish


def local_homfinder_mouse(ens_string):
    hl_mouse2human = []
    hl_mouse2rat = []
    hl_mouse2zebrafish = []

    hom_indexlist1 = list(np.where(mouse2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_mouse2human.append(mouse2human.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(mouse2rat["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_mouse2rat.append(mouse2rat.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(mouse2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_mouse2zebrafish.append(mouse2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_mouse2human, hl_mouse2rat, hl_mouse2zebrafish


def local_homfinder_zebrafish(ens_string):
    hl_zebrafish2human = []
    hl_zebrafish2rat = []
    hl_zebrafish2mouse = []

    hom_indexlist1 = list(np.where(zebrafish2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_zebrafish2human.append(zebrafish2human.iloc[hom_indexlist1[a], 1])
        print(hl_zebrafish2human)
    hom_indexlist2 = list(np.where(zebrafish2rat["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_zebrafish2rat.append(zebrafish2rat.iloc[hom_indexlist2[b], 1])
        print(hl_zebrafish2rat)
    hom_indexlist3 = list(np.where(zebrafish2mouse["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_zebrafish2mouse.append(zebrafish2mouse.iloc[hom_indexlist3[c], 1])
        print(hl_zebrafish2mouse)
    return hl_zebrafish2human, hl_zebrafish2rat, hl_zebrafish2mouse

#These four functions call the previous local_homfinder functions and use them to identify all of the homologs for all genes in an
#inputted dataframe.
def local_multi_species_homfinder_zeb(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "rat", "mouse"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_zebrafish(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df

def local_multi_species_homfinder_hum(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["rat", "mouse", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_human(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


def local_multi_species_homfinder_rat(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "mouse", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_rat(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


def local_multi_species_homfinder_mouse(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "rat", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_mouse(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


#function to insert blank rows between all values of dataframe, to pre-allocate rows
#for finding gene names and pubmed counts
def pir(df):
    nans = np.where(np.empty_like(df.values), np.nan, np.nan)
    data = np.hstack([nans, df.values]).reshape(-1, df.shape[1])
    return pd.DataFrame(data, columns=df.columns)


#function for inserting blank rows, restructuring dataframe to fit future data
def restruc_df(test_df):
    test_df = pir(pir(test_df))
    test_df = test_df.append(test_df.iloc[0:3])
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)
    test_df = test_df.drop([0,1,2], axis=0)
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)

    #remove the extra NaNs so that there is only 2 NaN's after every orthologue
    #after this, we can start filling in the NaN's with the actual gene names
    for y in range(3, len(test_df), 3):
        try:
            test_df = test_df.drop(test_df.index[y])
        except IndexError:
            break

    #reset the index
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)

    return test_df

#local ensembl ID finder. After restructuring the dataframe to add spaces between each row item, we call this function
#to do find the actual gene name for each ensembl ID in the dataframe, and add it to the cell directly underneath the respective ensembl_ID
#then output the dataframe.
def local_ens_id_finder(id_table):
    for a in range(0, len(id_table), 3):
        for b in range(0, len(id_table.columns)):
            if type(id_table.iloc[a][b]) == str:

                if id_table.iloc[a][b].startswith('ENSDARG'):

                    genindex = list(np.where(zebrafish_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = zebrafish_genelist.iloc[genindex[0]][1]


                elif id_table.iloc[a][b].startswith('ENSG'):

                    genindex = list(np.where(human_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = human_genelist.iloc[genindex[0]][1]

                elif id_table.iloc[a][b].startswith('ENSRNOG'):

                    genindex = list(np.where(rat_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = rat_genelist.iloc[genindex[0]][1]

                elif id_table.iloc[a][b].startswith('ENSMUSG'):

                    genindex = list(np.where(mouse_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = mouse_genelist.iloc[genindex[0]][1]

    return id_table

#this function finds the number of pubmed articles in for each gene in the dataframe + your search term keyword

def pubcrawl(dataframe, keyword):
    #loop through the dataframe every 3 rows, and for every column
    for x in tqdm(range(1, len(dataframe), 3)):
        for y in range(0 ,len(dataframe.columns)):
            #if the gene name is NaN or a 0, then continue on
            if pd.isnull(dataframe.iloc[x,y]) or dataframe.iloc[x,y] == 0:
                continue
            #search the entrez eutility pubmed API for each gene plus your keyword
            res = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + dataframe.iloc[x, y] + '+'+ keyword)
            ncbi = bs4.BeautifulSoup(res.text)
            count = ncbi.select('Count')
            
            #if nothing is found or there is an error, simply make the number of articles as 0, and continue
            if ncbi.select('errorlist') != []:
                dataframe.iloc[x+1, y] = 0
                continue
            #if nothing is found, make the number of articles as 0 and continue
            if count == []:
                dataframe.iloc[x+1,y] = 0
                continue
            #once you have the number of papers for your search term, conver to string, remove the html tags, and convert to integer
            #then add it to the cell directly underneath the respective gene in the dataframe
            once = str(count[0])
            once = once.replace('<count>', '')
            once = once.replace('</count>', '')
            dataframe.iloc[x+1,y] = int(once)
    # create a last column for totals
    dataframe = dataframe.assign(Total= np.nan)
    # sum up counts for all genes at the end and place in the last column for each row
    for p in range(2, len(dataframe), 3):
        dataframe.iloc[p, len(dataframe.columns) - 1] = dataframe.iloc[p, :].sum(axis=0)

    # add counts to the end of every row to help sort
    dataframe.Total = dataframe.Total.bfill()
    # create a helper column and fill with numbers to help sort but keep dataframe formatting
    dataframe["helper"] = np.arange(len(dataframe))//3
    # sort by helper values in descending order (highest to lowest pubmed article counts)
    dataframe = dataframe.sort_values(['Total', 'helper'], ascending = False)
    # remove the helper column
    dataframe = dataframe.drop(columns = "helper")
    # remove excess counts
    for p in range(2, len(dataframe), 3):
        dataframe.iloc[p - 1, len(dataframe.columns) - 1] = np.nan
        dataframe.iloc[p - 2, len(dataframe.columns) - 1] = np.nan
    dataframe = dataframe.reset_index()
    dataframe = dataframe.drop(dataframe.columns[0], axis=1)
    return dataframe



#HOW TO RUN THIS PROGRAM.
#RUN THE PYTHON FILE ON COMMAND LINE WITH THE FIRST ARGUMENT BEING THE FILE NAME OF YOUR EXCEL FILE CONTAINING YOUR INPUT GENE LIST IN ENSEMBL ID FORMAT
#THE SECOND ARGUMENT SHOULD BE THE SEARCH TERM YOU WISH TO USE WITH EACH GENE AGAINST THE PUBMED DATABASE
#THE THIRD ARGUMENT SHOULD BE THE DESIRED NAME OF YOUR OUTPUT EXCEL FILE e.g. "output-excel-file.xlsx" OR WHATEVER YOU WISH TO CALL IT

def main():
    my_parser = argparse.ArgumentParser(description="This program takes as input a zebrafish, human, rat, or mouse ENSEMBL ID gene list and a search term, the returns all orthologues of those genes,"
                                                    "and searches each gene plus your keyword into pubmed, "
                                                    "returning the number of articles published for that particular search term.")
    my_parser.add_argument('Input Excel File', metavar='Input-File', type=str, help='the path to the gene list excel file')
    my_parser.add_argument('Search Term(s)', metavar='Search-Term', type=str,
                           help='the search term to search pubmed alongside each gene in the list. If your search term is more than one word, please enclose the entire search term in quotes')
    my_parser.add_argument('Output Excel File Name', metavar='Output-File-Name', type=str,
                           help='the name of the excel file you would like to output')
    args = my_parser.parse_args()
    #Check if the file has all zebrafish, human, rat, or mouse Ensemble IDs and is fully homogenous with these IDs.
    #If they do, then run the program with the correct homolog finder function for the species.
    if orth_df['Name'].str.startswith('ENSDARG').all():
        local_multi_species_homfinder_zeb(orth_df)
        test_df = restruc_df(orth_df)
        test_df = local_ens_id_finder(test_df)
        test_df = pubcrawl(test_df, sys.argv[2])
        test_df.to_excel(sys.argv[3])
        print('The program took {0} seconds !'.format(time.time() - startTime))
    elif orth_df['Name'].str.startswith('ENSG').all():
        local_multi_species_homfinder_hum(orth_df)
        test_df = restruc_df(orth_df)
        test_df = local_ens_id_finder(test_df)
        test_df = pubcrawl(test_df, sys.argv[2])
        test_df.to_excel(sys.argv[3])
        print('The program took {0} seconds !'.format(time.time() - startTime))
    elif orth_df['Name'].str.startswith('ENSRNOG').all():
        local_multi_species_homfinder_rat(orth_df)
        test_df = restruc_df(orth_df)
        test_df = local_ens_id_finder(test_df)
        test_df = pubcrawl(test_df, sys.argv[2])
        test_df.to_excel(sys.argv[3])
        print('The program took {0} seconds !'.format(time.time() - startTime))
    elif orth_df['Name'].str.startswith('ENSMUSG').all():
        local_multi_species_homfinder_mouse(orth_df)
        test_df = restruc_df(orth_df)
        test_df = local_ens_id_finder(test_df)
        test_df = pubcrawl(test_df, sys.argv[2])
        test_df.to_excel(sys.argv[3])
        print('The program took {0} seconds !'.format(time.time() - startTime))
    else:
        print("The input excel file does not contain homogenous zebrafish, human, rat, or mouse ENSEMBL IDs or is not supported by this program.")
        print("Please input a valid excel file.")
main()