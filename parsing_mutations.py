#############################################################################################################
# functions
#############################################################################################################
import re 

def extract_digit_from_string(string):
    digit_list = [str(s) for s in list(string) if s.isdigit()]
    digit = ''.join(digit_list)
    return digit


def extract_nondigit_from_string(string):
    nondigit_list = [str(s) for s in list(string) if not s.isdigit()]
    nondigit = ''.join(nondigit_list)
    return nondigit


def DnaCheck(sequence):
    return all(base.upper() in ('A', 'C', 'T', 'G') for base in sequence)


def check_aa_mutation(mutation):
    """ Returns true if amino acid mutation meets expected format (e.g. V720G), else false """
    right_format = False
    items = list(mutation)
    ref_aa = items[0]
    pos_aa = ''.join(items[1:len(items)-1])
    if pos_aa.isdigit():
        if ref_aa.isupper():
            right_format = True
    return right_format

def reformat_nt_mutation(mutation):
    return mutation.lower()

def check_nt_mutation(mutation):
    """ Returns true if nucleotide mutation meets expected format (e.g. G2299T), else false """
    right_format = False
    items = list(mutation)
    ref_nt = items[0]
    pos_nt = ''.join(items[1:len(items)-1])
    nucleotides = ["A", "T", "C", "G", "a", "t", "c", "g"]
    if pos_nt.isdigit():
        if ref_nt in nucleotides:
            right_format = True
    return right_format


def reformat_promoter_mutation(mutation):
    """ reformat promoter mutations to meet TBprofiler format (e.g. [-11 C>A] to -11C>A) """
    items = mutation.replace("upstream ", "").split(' ')
    pos_nt = items[0].replace("[", "")
    change_nt = items[1].split("]")[0]
    new_mutation = pos_nt + change_nt
    return new_mutation


def check_promoter_mutation(mutation):
    """ Returns true if promoter mutations meets TBprofiler format (-11C>A) """
    right_format = False
    items = mutation.split(">")
    if len(items) == 2:
        alt_nt = items[1]
        items2 = list(items[0])
        pos_nt = ''.join(items2[0:len(items2)-1])
        pos_nt = pos_nt.replace("-", "")
        ref_nt = items2[-1]
        nucleotides = ["A", "T", "C", "G"]
        right_format = True
                    
    return right_format


def reformat_ins_nt_mutation(mutation):
    """ reformat insertions with format like [Ins A nt38–39] """
    items = mutation.replace("[Ins ", "").replace("[Ins", "").split(']')
    items2 = items[0].split(' ')
    seq_nt = items2[0].lower()  # nucleotide are represented with lower case in this script
    tmp = items2[1].replace('nt', '')
    pos_nt = tmp.split('-')[0].split('–')[0]
    new_mutation = 'x' + str(pos_nt) + '_ins_' + seq_nt.replace('insertion', '')
    return new_mutation


def reformat_del_nt_mutation(mutation):
    """ reformat insertions with format like [Del gg nt18-19] """
    items = mutation.replace("[Del ", "").replace("[Del", "").split(']')
    items2 = items[0].split(' ')
    seq_nt = items2[0].lower()  # nucleotide are represented with lower case in this script
    tmp = items2[1].replace('nt', '')
    pos_nt = tmp.split('-')[0].split('–')[0]
    new_mutation = 'x' + str(pos_nt) + '_del_' + seq_nt.replace('deletion', '')
    return new_mutation


def reformat_ins_mutation(mutation):
    """ reformat insertions with format like [Ins139G] """
    tmp = mutation.replace(" ins", "").replace(" Ins", "")#.split(']')[0]
    if '_' in mutation:
        mutation.split('_')[1]
        pos_nt = extract_digit_from_string(tmp)
        seq_nt = extract_nondigit_from_string(tmp).lower()
        new_mutation = 'x' + str(pos_nt) + '_ins_' + seq_nt
    else:   
        pos_nt = extract_digit_from_string(tmp)
        seq_nt = extract_nondigit_from_string(tmp).lower()
        new_mutation = 'x' + str(pos_nt) + '_ins_' + seq_nt
    return new_mutation


def reformat_del_mutation(mutation):
    """ reformat insertions with format like [Del336C] """
    tmp = mutation.replace("[Del ", "").replace("[Del", "").split(']')[0]
    pos_nt = extract_digit_from_string(tmp)
    seq_nt = extract_nondigit_from_string(tmp).lower()
    new_mutation = 'x' + str(pos_nt) + '_del_' + seq_nt
    return new_mutation


def reformat_insertion_mutation(mutation):
    """ reformat insertions with format like [G193 insertion] or [2734AACTT insertion] """
    tmp = mutation.replace(" insertion", "").replace(" position", "").replace(" ins", "").replace("ins", "").replace(" Ins", "").replace("Ins", "").replace(" -", "").replace(" ", "").replace("of", "").replace("_ins_", "")
    pos_nt = extract_digit_from_string(tmp)
    seq_nt = extract_nondigit_from_string(tmp).lower().replace("_", "")
    new_mutation = 'x' + str(pos_nt) + '_ins_' + seq_nt
    if '_' in mutation:
        tmp2 = tmp.split('_')[0]
        #sometimes, starts with nuclotides then has the numbers then ins eg 'CC816_817 insertion '
        tmp3 = extract_digit_from_string(tmp2)
        new_mutation = 'x' + str(tmp3) + '_ins_' + seq_nt
    return new_mutation


def reformat_deletion_mutation(mutation):
    """ reformat deletions with format like [CGCTGGGC371–378 deletion] or [491AACGGGG deletion] """
    tmp = mutation.replace(" deletion", "").replace("Deletion", "").replace("del", "").replace("Del", "").replace(" position", "").replace(" in", "").replace(" codon", "").replace(" del", "").replace("del", "").replace(" Del", "").replace(" -", "").replace(" ", "").replace("of", "").replace("_del_", "")
    pos_nt = extract_digit_from_string(tmp)
    seq_nt = extract_nondigit_from_string(tmp).lower().replace("_", "")
    new_mutation = 'x' + str(pos_nt) + '_del_' + seq_nt
    if '_' in mutation:
        tmp2 = tmp.split('_')[0]
        #sometimes, starts with nuclotides then has the numbers then del eg 'CC816_817 del'
        tmp3 = extract_digit_from_string(tmp2)
        new_mutation = 'x' + str(tmp3) + '_del_' + seq_nt
    return new_mutation

def reformat_arrow_mutation(mutation):
    tmp = mutation.replace(" at", "")
    tmp2 = tmp.split('→')
    #pos_nt = [int(i) for i in tmp2 if i.isdigit()]
    #pos_nt = [int(s) for s in re.findall(r'\d+', tmp2)]
    #pos_nt2 =''.join(map(str,pos_nt))
    pos_nt = extract_digit_from_string(tmp2)
    seq_nt = extract_nondigit_from_string(tmp2[0]).lower()
    seq_nt2 = extract_nondigit_from_string(tmp2[-1]).lower()
    new_mutation = seq_nt + pos_nt + seq_nt2
    return new_mutation

def check_indel_mutation(mutation):
    """ Returns true if mutation meets expected insertion format (e.g. x2734_ins_aactt), else false """
    right_format = False
    items = mutation.split('_')
    labels = ['ins', 'del']
    if len(items) == 3:
        pos_nt = items[0].replace('x', '')
        label = items[1]
        ins_nt = items[2]
        if pos_nt.isdigit():
            if label in labels:
                if DnaCheck(ins_nt):
                    right_format = True
    return right_format

#############################################################################################################
# Main function
#############################################################################################################

def _main():
    # output files
    non_parse_file = open("non_parsed.txt", "w")
    parsed_file = open("parsed.txt", "w")
    
    with open('raw_mutations_database.tsv', 'r') as mut_lines:
        # to skip first line/header of file
        next(mut_lines)
        for mut_line in mut_lines:
            (gene, mutation, drug, reference, resistance, *_) = mut_line.strip().split('\t')
            
            #want to keep only cases where resistance = R 
            #removing whitespace
            resistance = resistance.strip()
            
            #keeping only cases where resistance is 'R', ie removing any uncertain or 's'
            if resistance == 'R': 
                original_mutation = mutation
                new_mutation = ""
                checked = False     # to check whether mutation has been checked in this scripts
                parsed = False     # whether mutation has been successfully parsed


            # standardising gene names for the annotations -> using gene name instead of locus name, but if novel gene, then official gene name is the locus name 
            ## for all genes, use the name of genes as they are in the annotation file
            
            #removing whitespace at the end of the gene name as sometimes in google sheets an extra whitespace is added at the end 
                gene = gene.rstrip()
    
            
            ## 23SRna == rrl 
                if gene == '23S rRNA (rr/)': 
                    gene = gene.replace("23S rRNA (rr/)", "rrl")
                
                if gene == '23S rRNA':
                    gene = gene.replace("23S rRNA", "rrl")
            
            ## mmpR5 and mmpR == Rv0678, Novel gene so no official name in annotation, use locus name     
                if gene == 'mmpR5':
                    gene = gene.replace("mmpR5", "Rv0678") 
                
                if gene == 'mmpR':
                    gene = gene.replace("mmpR", "Rv0678")
                #spelling of the gene 
                if gene == 'rv0678':
                    gene = gene.replace("rv0678", "Rv0678")
                
           ## Rv1979 is supposed to be Rv1979c but there are some which are already in the correct format
                if gene == 'Rv1979':
                    gene = gene.replace("Rv1979", "Rv1979c") 
                
            ## Rv0676c is supposed to be mmpL5 
                if gene =='Rv0676c':
                    gene = gene.replace("Rv0676c", "mmpL5") 
                 
                if gene =='Rv0676c (mmpL5)':
                    gene = gene.replace("Rv0676c (mmpL5)", "mmpL5")      
                
            ## Rv0677c (mmpS5) should simply be mmpS5 
                if gene == 'Rv0677c (mmpS5)':
                    gene = gene.replace("Rv0677c (mmpS5)", "mmpS5") 
                
                if gene == 'Rv0677c' :
                    gene = gene.replace("Rv0677c", "mmpS5")
                
            ## Rv1304 should be atpB 
                if gene == 'Rv1304 (atpB)':
                    gene = gene.replace("Rv1304 (atpB)", "atpB")
                
            ## Rv2535c should be pepQ 
                if gene == 'Rv2535c' :
                    gene = gene.replace("Rv2535c", "pepQ")
                
                if gene == 'Rv2535c (pepQ)':
                    gene = gene.replace("Rv2535c (pepQ)", "pepQ")
            
            ## fbiD should be Rv2983
                if gene == 'fbiD' :
                    gene = gene.replace("fbiD", "Rv2983")
            
                
                     
            # Pre-editing tasks of mutation variable
            #   Removing leading white space
                mutation = mutation.lstrip()
            #   Replacing multiple whitespaces with only one: will make parsing mutation easier later on
                mutation = ' '.join(mutation.split())
            # some mutations start with 'c.' and others start with 'p.', need to remove this 
                mutation = mutation.replace("c.", "").replace("p.", "")
            
            
            #there seems to be 5 types of mutations:
            # 1) Nucleotide substitutions : they start with and end with a nucloetide with numbers in between them eg C110A
            # 2) promotors -> start with '-', some are already in correct format eg -11 C>A
            # 3) insertions and deleteions (nucleotides) 
            # 4) AA sustitutions in short form -> assumption unlikely to start with and end with the bases eg M146T
            # 5) AA in long format eg Arg126His

            
           
            # 1) dealing with nucleotide substitutions 
            # assumption: if a mutation starts with a base and ends with a base, it is a nucleatide 
                nucleotides = ["A", "C", "T", "G", "a", "c", "t", "g"]
                if (mutation.startswith(tuple(nucleotides))) and (mutation.endswith(tuple(nucleotides))):
                    print(mutation)
                    checked = True
                #if nucleotide mutation, want it to have lower cases to diffrentiate from AA 
                    if check_nt_mutation(mutation):
                        new_mutation = reformat_nt_mutation(mutation)
                        print("\t --> Nucleotide mutation " + new_mutation)
                        parsed = True    
                        
            # 2) nuclotide promoters that are already in correct format eg '-11 C>A'
                if (mutation.startswith('-')) and ('>' in mutation) :
                    print(mutation)
                    checked = True
                    new_mutation = mutation
                    if check_promoter_mutation(new_mutation):
                    #making the promotors in lower case to show they are nucleotide 
                        new_mutation = mutation.lower()
                        parsed = True
                        print("\t --> Promoter mutation " + new_mutation)
                    else:
                        print("\t --> Wrong promoter nucleotide: not parsed")
                    

            # 3a) dealing with insertions
                if 'ins' in mutation.lower():
                    checked = True
                # within insertions starting with '[Ins', a subset contain ' nt'
                    if ' nt' in mutation:
                        new_mutation = reformat_ins_nt_mutation(mutation)
                    else:
                        new_mutation = reformat_insertion_mutation(mutation)
                # finally check correct insertion format before saving
                    print(mutation)
                    print(new_mutation)
                    if check_indel_mutation(new_mutation):
                        print("\t --> Insertion mutation " + new_mutation)
                        parsed = True
                    else:
                        print("\t --> Wrong insertion format: not parsed")

            # 3b) dealing with deletions (same patterns expected as for insertions)
                if 'del' in mutation.lower():
                    checked = True
                # within deletions starting with '[Del', a subset contain ' nt'
                    if ' nt' in mutation:
                        new_mutation = reformat_del_nt_mutation(mutation)
                    else:
                        new_mutation = reformat_deletion_mutation(mutation)
                # finally check correct deletion format before saving
                    print(mutation)
                    print(new_mutation)
                    if check_indel_mutation(new_mutation):
                        print("\t --> Deletion mutation " + new_mutation)
                        parsed = True
                    else:
                        print("\t --> Wrong deletion format: not parsed")
                    
            #dealing with AA mutations in long form format 
                amino_acid = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
                if any(amino in mutation for amino in amino_acid):
                    checked = True 
                    print(mutation)
                #some of these, start with 'p.', need to remove this from the mutation
                    new_mutation = mutation.replace(".p", "") 
                    print("\t --> Full AA mutation " + new_mutation)
                    parsed= True
            
            # dealing with AA in short format 
                short_form = ["R", "N", "D", "E", "Q", "H", "I", "L", "K", "M", "F", "P", "S", "W", "Y", "V"]
                if any(short in mutation for short in short_form):
                    checked = True 
                    print(mutation)
                    if check_aa_mutation(mutation):
                        new_mutation = mutation 
                        print("\t --> short AA mutation " + new_mutation)
                        parsed= True
                
                

            # and finally save mutation or new_mutation in the right output file
                if not checked:
                    non_parse_file.write(mut_line.strip() + '\t' + 'mutation_not_checked' + '\n')
            
                else:
                    if parsed:
                    # saving parsed mutations
                        parsed_mutation = gene + "@" + new_mutation 
                    #parsed_file.write(original_mutation + "\t" + parsed_mutation + "\t" + drug + "\t" + reference + '\n')
                        parsed_file.write(parsed_mutation + "\t" + drug + "\t" + reference + '\t' + resistance + '\n' )
                
                    else:
                        non_parse_file.write(mut_line.strip() + '\t' + 'mutation_not_parsed' + '\n')


    parsed_file.close()
    non_parse_file.close()