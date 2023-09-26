def _main():
    with open("tb_profiler_mutations.txt", "w") as tb:
        with open("parsed.annotated.txt", "r") as annotations:
            for annotation in annotations:
                (Mutation, Drug, Reference, Resistance) = annotation.strip().split('\t')

                list = Mutation.split("@")
                Gene = list[0]
                Mut = list[1]
                profiler = ""
                import re

                # substitute all letters to 3 letter abbreviations for amino acid substitutions
                if "A" in Mut[0]:
                    Mut = Mut.replace("A", "Ala")
                elif "R" in Mut[0]:
                    Mut = Mut.replace("R", "Arg")
                elif "N" in Mut[0]:
                    Mut = Mut.replace("N", "Asn")
                elif "D" in Mut[0]:
                    Mut = Mut.replace("D", "Asp")
                elif "C" in Mut[0]:
                    Mut = Mut.replace("C", "Cys")
                elif "E" in Mut[0]:
                    Mut = Mut.replace("E", "Glu")
                elif "Q" in Mut[0]:
                    Mut = Mut.replace("Q", "Gln")
                elif "G" in Mut[0]:
                    Mut = Mut.replace("G", "Gly")
                elif "H" in Mut[0]:
                    Mut = Mut.replace("H", "His")
                elif "I" in Mut[0]:
                    Mut = Mut.replace("I", "Ile")
                elif "L" in Mut[0]:
                    Mut = Mut.replace("L", "Leu")
                elif "K" in Mut[0]:
                    Mut = Mut.replace("K", "Lys")
                elif "M" in Mut[0]:
                    Mut = Mut.replace("M", "Met")
                elif "F" in Mut[0]:
                    Mut = Mut.replace("F", "Phe")
                elif "P" in Mut[0]:
                    Mut = Mut.replace("P", "Pro")
                elif "S" in Mut[0]:
                    Mut = Mut.replace("S", "Ser")
                elif "T" in Mut[0]:
                    Mut = Mut.replace("T", "Thr")
                elif "W" in Mut[0]:
                    Mut = Mut.replace("W", "Trp")
                elif "Y" in Mut[0]:
                    Mut = Mut.replace("Y", "Tyr")
                elif "V" in Mut[0]:
                    Mut = Mut.replace("V", "Val")
                elif "rrl" in Gene:
                    Mut = Mut
                else:
                    print("Error")

                if "Glylu138" in Mut:
                    Mut = Mut.replace("Glylu", "Glu")
                else:
                    Mut = Mut

                # substituting the end nucleotide in the mutation  p.His722R
                if "A" in Mut[-1]:
                    Mut = Mut.replace("A", "Ala")
                elif "R" in Mut[-1]:
                    Mut = Mut.replace("R", "Arg")
                elif "N" in Mut[-1]:
                    Mut = Mut.replace("N", "Asn")
                elif "D" in Mut[-1]:
                    Mut = Mut.replace("D", "Asp")
                elif "C" in Mut[-1]:
                    Mut = Mut.replace("C", "Cys")
                elif "E" in Mut[-1]:
                    Mut = Mut.replace("E", "Glu")
                elif "Q" in Mut[-1]:
                    Mut = Mut.replace("Q", "Gln")
                elif "G" in Mut[-1]:
                    Mut = Mut.replace("G", "Gly")
                elif "H" in Mut[-1]:
                    Mut = Mut.replace("H", "His")
                elif "I" in Mut[-1]:
                    Mut = Mut.replace("I", "Ile")
                elif "L" in Mut[-1]:
                    Mut = Mut.replace("L", "Leu")
                elif "K" in Mut[-1]:
                    Mut = Mut.replace("K", "Lys")
                elif "M" in Mut[-1]:
                    Mut = Mut.replace("M", "Met")
                elif "F" in Mut[-1]:
                    Mut = Mut.replace("F", "Phe")
                elif "P" in Mut[-1]:
                    Mut = Mut.replace("P", "Pro")
                elif "S" in Mut[-1]:
                    Mut = Mut.replace("S", "Ser")
                elif "T" in Mut[-1]:
                    Mut = Mut.replace("T", "Thr")
                elif "W" in Mut[-1]:
                    Mut = Mut.replace("W", "Trp")
                elif "Y" in Mut[-1]:
                    Mut = Mut.replace("Y", "Tyr")
                elif "V" in Mut[-1]:
                    Mut = Mut.replace("V", "Val")
                elif "*" in Mut[-1]:
                    Mut = Mut
                elif "rrl" in Gene:
                    Mut = Mut
                else:
                    Mut = Mut
                    #print("Error: Mut[-1] contains unrecognised format ->>" + Mut[-1] + "\n" + "Mut is : " + Mut)

                if Gene == 'rrl':
                    loc = re.findall('\d+', Mut)
                    nuc = Mut.split(''.join(loc))

                    profiler = ''.join(Gene) + "\t" + "r." + ''.join(loc) + ''.join(nuc[0]) + ">" + ''.join(nuc[1]) + "\t" \
                               + Drug + "\t" + Reference + "\t" + Resistance
                    print(profiler, file=tb)
                       
                #needing to account for nucleotide mutations 
                elif Mut[0].islower() and Gene != 'rrl':
                    profiler = ''.join(Gene) + "\t" + "c." + Mut + "\t" \
                               + Drug + "\t" + Reference + "\t" + Resistance
                    print(profiler, file=tb)
                    

                else:
                    loc1 = re.findall('\d+', Mut)
                    aa = Mut.split(''.join(loc1))

                    profiler = ''.join(Gene) + "\t" + "p." + ''.join(aa[0]) + ''.join(loc1) + ''.join(aa[1]) + "\t" \
                               + Drug + "\t" + Reference + "\t" + Resistance
                    print(profiler, file=tb)

    tb.close()
