#loading libraries 
import argparse
import logging
import os
import sys
import Bio
from Bio import SeqIO
from Bio.Seq import Seq

#main function 
def _main():
    # output files
    annotated_file = open("parsed.annotated.txt", "w")
    non_annotated_file = open("parsed.not_annotated.txt", "w")

    # Saving the full path of genes
    gene_fasta_files = dict()
    with open('input_genes.tsv', 'r') as input_genes:
        for line in input_genes:
            (gene_name, sequence_file) = line.strip().split('\t')
            gene_fasta_files[gene_name] = sequence_file
            print("gene_fasta_files[" + gene_name + "] --> " + gene_fasta_files[gene_name])

    # Reading mutations
    with open('input_mutations.txt', 'r') as input_mutations:
        for line in input_mutations:
            (mut, drug, reference, resistance) = line.strip().split("\t")
            mutations = mut.split('+')
            checked = False
            annotated = False
            for mutation in mutations:
                # extract mutation's gene name
                (gene_name, change) = mutation.strip().split('@')
                # dealing with genes that are a range as they are being flagged as a problem:
                if '-' in gene_name:
                    print(mutation)
                    # extract mutation's reference allele and position
                    ref = list(change)[0]
                    print("Range of genes given")
                    checked = True 
                    #mutations where genes have been provided as a range (intergenic regions -> cannot annotate)
                    non_annotated_file.write(mutation + "\t" + drug + "\t" + reference +  " genes are in a range" + "\n" )
                    
                
                if gene_name != 'na' and '-' not in gene_name:
                    print(mutation)
                    #removing any whitespace in gene name 
                    gene_name = gene_name.replace(" ", "")
                    # extract mutation's reference allele and position
                    ref = list(change)[0]
                    # determine whether nucleotide or amino acid change
                    type = 'unknown'
                    if ref.islower() and ref != 'x':
                        type = 'nucleotide'
                    if ref.isupper():
                        type = 'amino_acid'
                    if ref == 'x':
                        type = 'indel'
                    if type == 'unknown':
                        print('Unknown reference type: ' + ref + ' in ' + change)

                    # extract positition
                    digits = [int(s) for s in list(change) if s.isdigit()]
                    pos = ''
                    for digit in digits:
                        pos += str(digit)
                    
                    # make sure gene's corresponding FASTA file exists
                    if gene_name not in gene_fasta_files:
                        logging.error(f"{gene_name} (mutation {change}) not found in {input_genes}")
                    


                    # read gene's FASTA file and translate to amino acid sequence
                    first_record = next(SeqIO.parse(gene_fasta_files[gene_name], "fasta"))
                    gene_BioSeq_nt = first_record.seq
                    gene_BioSeq_aa = gene_BioSeq_nt.translate()

                    # compare mutation's reference allele to gene's FASTA
                    ref_fasta = ''
                    if type == 'amino_acid':
                        print(pos)
                        if int(pos) > len(gene_BioSeq_nt):
                            print("Mutation location incorrect length")
                            checked = True
                        else:
                            try:
                                ref_fasta = gene_BioSeq_aa[int(pos)-1]
                                print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                                checked = True
                            except IndexError as e:
                                print("Index Error")

                    if type == 'unknown':
                        print(pos)
                        if int(pos) > len(gene_BioSeq_nt):
                            print("Mutation location incorrect length")
                            checked = True
                        else:
                            try:
                                ref_fasta = gene_BioSeq_aa[int(pos) - 1]
                                print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                                checked = True
                            except IndexError as e:
                                print("Index Error")

                    if type == 'nucleotide':
                        print(pos)
                        if int(pos) > len(gene_BioSeq_nt):
                            print("Mutation location incorrect length")
                            checked = True
                        else:
                            try:
                                ref_fasta = gene_BioSeq_nt[int(pos)-1]
                                ref = ref.upper()
                                print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                                checked = True 
                            except IndexError as e:
                                print("Index Error")
                        
                        
                    if type == 'indel': 
                        print(pos)
                        if int(pos) > len(gene_BioSeq_nt):
                            print("Mutation location incorrect length")
                            checked = True
                        else:
                            try:
                                ref_fasta = gene_BioSeq_nt[int(pos)-1]
                                #want the reference to actually be the final letter 
                                #ref = list(change)[-1]
                                change = change.replace("_ins_", "").replace("_del_", "").replace("x", "").replace(pos, "")
                                #somestimes multiple letters are given, want to use the first one as reference 
                                ref = list(change)[0]
                                ref = ref.upper()
                                print('change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                                checked = True 
                            except IndexError as e:
                                print("Index Error")
                    
                    
                    if ref == ref_fasta or '_ins_' in mutation:   #otherwise print into separate file
                        #insertions = cannot 'annotate' so accepting as truth?
                        annotated_file.write(mutation + "\t" + drug + "\t" + reference + "\t" + resistance + "\n")
                        #print('mutation ' + mutation + ' pos ' + str(pos) + ' change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                
                    #with insertions need to assume they're correct, because ref might equal ref_fasta or not. either way, cannot tell         

                    else:
                        print('mutation ' + mutation + ' pos ' + str(pos) + ' change ' + change + ' ref ' + ref + ' ref_fasta ' + ref_fasta)
                        if not checked: 
                        #keeping track of mutations that cause an index error, used exceptions to bypass them but keeping track of occurances with "not checked"
                            non_annotated_file.write(mutation + "\t" + drug + "\t" + reference + "\t" + resistance + " not checked" + "\n" )
                        else: 
                            non_annotated_file.write(mutation + "\t" + drug + "\t" + reference + "\t" + resistance + "\n")
                            



    non_annotated_file.close()
    annotated_file.close()
