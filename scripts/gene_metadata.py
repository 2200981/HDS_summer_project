#main function 
def _main():
    #output file 
    gene_data_file = open("gene_metadata.txt", "w")
               
  
    # Getting the name of genes required and saving this into a list
    with open('required_genes.tsv', 'r') as input_genes:
        # to skip first line/header of file
        next(input_genes)
        gene_names = []
        for line in input_genes:
            (gene_name, locus_tag) = line.strip().split('\t')
            #gene_names = gene_name.split('+')
            gene_name = gene_name.strip()
            gene_names.append(gene_name) 
        
                    
    print(gene_names)               
    # reading in the gff file  
    with open('Mycobacterium_tuberculosis_H37Rv_gff_v4.gff', 'r') as gene_infos:
        # to skip first line/header of file
        next(gene_infos)
        official_names = []
        for gene_info in gene_infos:
            (nc, myco, output, start_coordinate, end_coordinate, dot, pos_neg, num, metadata, *_) = gene_info.strip().split('\t')
            # each line looks like: NC_000962.3	Mycobrowser_v4	tRNA	2794176	2794249	.	+	0	Locus=MTB000032;Name=argW;
            #extract gene name from metadata
            official_name = metadata.split(';')[1].replace('Name=', '')
            official_names.append(official_name)   
            if official_name in gene_names:
                print(start_coordinate, end_coordinate, pos_neg, official_name)
                gene_data_file.write(start_coordinate + "\t" + end_coordinate + "\t" + pos_neg + "\t" + official_name + "\n")

    
    
    gene_data_file.close()
