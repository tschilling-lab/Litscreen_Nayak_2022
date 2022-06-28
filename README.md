# Litscreen v1.0

Litscreen finds all orthologous genes in different species for every gene (currently only ENSEMBL ID inputs are accepted) in a user input gene list, converting the ENSEMBL IDs into conventional gene names, and then performing a search for a user input keyword (for example "ribosomal RNA" or "cartilage development") against the Pubmed database via the [Pubmed API](https://www.ncbi.nlm.nih.gov/books/NBK25501/) for each gene and ortholog with the user inputted search terms.

This is a python command line program, called litscreen.py, and takes 3 arguments as input:   

1) The file name of the excel file containing the input gene list in ENSEMBL ID format. The excel file should be formatted as shown in the following screenshot (Additional examples of input can be seen in the Sample Input Folder on this Github Repository):

![Sample Input](/Screenshots/Cartilage-Input.PNG)

2) The search term(s) you wish to search with each gene against the Pubmed API/Pubmed Research Database

3) The desired name of your output excel file, e.g. "output-file.xlsx"

A sample input would look like this, where the input ensembl gene list excel file to be searched is called test-ens-gene-list.xlsx, the search term for pubmed is "cartilage" and the output excel file name is "cartilage-output.xlsx" :

    python litscreen.py test-ens-gene-list.xlsx "cartilage" "cartilage-output.xlsx"


Additionally, the program currently requires that you have downloaded the "ensembl_gene_lists" folder (and all files within), as well as the "orthologous_gene_lists" folder (and all files within). After entering the arguments for the program, it will ask for the paths to these folders, and the path to the parent folder containing your input file as well.

you should input the paths without quotes, such as shown below:

    Please enter path to the ensembl_gene_lists folder: C:\Users\USER\Desktop\litscreen\ensembl_gene_lists\
    Please enter path to the orthologous_gene_lists folder: C:\Users\USER\Desktop\litscreen\orthologous_gene_lists\
    Please enter path to your gene list input excel file folder: C:\Users\USER\Documents\input_files\


 Once these have all been input, the program should run normally, and output a file similar to the following screenshot (Additional examples of output can be seen in the Sample Output Folder on this Github Repository):

![Sample Output 1](/Screenshots/Cartilage-Output-1.PNG)

![Sample Output 2](/Screenshots/Cartilage-Output-2.PNG)

and we can validate that the number of articles written for this particular search (the search term in this case is "cartilage") by manually searching a gene and the search term on Pubmed, and looking at the number of articles returned on the web page.

Example 1, compared with genes boxed on first output screenshot:
![Validation 1](/Screenshots/Validation-1.PNG)

Example 2:
![Validation 2](/Screenshots/Validation-2.PNG)

Screenshots of Program Working:


After running program in command line, entering all arguments and paths correctly, the program should display a progress bar via the tqdm module, shown below:

![Working-Program](/Screenshots/program_working.PNG)

The Input file is in the correct path, shown below:

![program-input](/Screenshots/program_input.PNG)

After completion, the program should display a completed progress bar and the time taken to complete the program, shown below:

![Completed-Program](/Screenshots/program_completed.PNG)

The output excel file should now have been created in the working directory:

![program-output](/Screenshots/program_output.PNG)

## Known Issues and Potential Future Directions

1) When Genes in the user inputted excel file are not found in the local data storage text file, the output is returned under that gene as GENE NOT FOUND IN DATABASE, see example below:

![Not Found](/Screenshots/Genes_not_found.PNG)

An improvement could be made here where the program would use ENSEMBL REST API requests only for genes not found in the local database, so that potentially useful genes are not missed in the search.

2) The program only takes ENSEMBL IDs as input currently, and only finds orthologous between zebrafish, human, rat and mouse. Interconversions between different gene ID formats can be considered in future, as well as the addition of orthologous genes from other organisms (e.g. Drosophila, C. Elegans, etc.)

3) Sometimes, the number of articles returned does not correspond to the number of articles written if you were to validate against the Pubmed website, and this is due to differences in what the Entrez Eutils API returns versus what the Pubmed website returns. However, I usually see that even if there is a discrepancy in the quantity of articles written, I have not seen a discrepancy where the API returns 0 articles, and the website returns a positive number of articles, so the program should still give some sort of useful information regardless of if the exact number of articles is 100% accurate.

4) Occasionally, if a search that returns 0 results on the Pubmed website will return an absurdly large number of articles (in the thousands to tens of thousands) in the program. This is again due to an issue with the Entrez Eutils API. A filter could be written to remove these abnormal results in future.

5) Despite being useful, the program could be more user friendly if turned into a web application with a User Interface to allow users who are not familiar with command line/programming usage to extract useful genes from their data.










