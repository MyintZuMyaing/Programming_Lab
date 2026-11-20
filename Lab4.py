# Input CDS FASTA file
#input_file = "Toxo_Ser_Thr_kinase.fasta"

# Output protein FASTA file
#output_file = "Toxo_Ser_Thr_kinase.txt"

# make a codon table
bases = "tcag".upper()
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
#print (codon_table["ATG"])

# codon_table is a dict containing codons and amino acids
# print it out to check you understand the format

# Translate each CDS sequence and save the protein sequences

Gene_Dic = {}

# Open the Text file and read its contents
with open ("GenesWithSignalPeptide_SignalP.txt","r") as my_file:
    #As the above table has header in line 1, the below command is needed to start from line 2 
    header = my_file.readline()
    for line in my_file:
        line = line.rstrip('\n')
        gene_id, peptide_start, peptide_end, transcipt_id= line.split('\t')
        #print (gene_id)
        
        #break to print only the first gene-id
        
        Gene_Dic[gene_id] = [peptide_start, peptide_end]
print (Gene_Dic["TGME49_200240"])
        

with open ("Toxo_Ser_Thr_kinase.fasta","r") as my_fastafile:
    for line in my_fastafile:
            line = line.rstrip('\n')
            if line.startswith(">"):
                gene_id = line[1:]
                #print (gene_id)
                #break
                
            else:
                DNA_seq = line
                #print (DNA_seq)
                #break
                
                if gene_id in Gene_Dic:
                    #print (gene_id)
                    #print (DNA_seq)
                    
                    #codons =[]
                    protein_seq=""
                    for i in range(0, len(DNA_seq), 3):
                        # the current value of the range is the start
                        # the end is i+3
                        codon = DNA_seq[i:i+3]
                        #print (codon)
                        #codons.append(codon)
                        #print (codons)
                        
                        amino_acid= codon_table[codon]
                        protein_seq= protein_seq+amino_acid
                        #print (gene_id)
                        #print (protein_seq)
                        
                        
                        # Extract the signal peptide
                        signal_peptide = protein_seq[0:29]
#print (signal_peptide)
                    
                        if signal_peptide in protein_seq:
                            print (signal_peptide)