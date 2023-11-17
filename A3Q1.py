import scipy 
import numpy as np 
import pandas as pd
import io 
import os

def read_vcf(population):

    #parse vcf file  

    #create dictionary to store each value in pertaining list 
    vcf_data = {'CHROM': [], 'POS': [], 'ID': [], 'REF': [], 'ALT': [], 'QUAL': []}

    with open(population, 'r') as vcf_file:
    
        for l in vcf_file:

            #extract SNP columns from header
            if l.startswith('#'): 
                header_fields= l.strip().split('\t')

                #take columns after all other fields in dictionary 
                snp_columns= header_fields[6:]

                #initialize to empty list 
                for col in snp_columns: 
                    vcf_data[col]= []


            elif not l.startswith("#"):

                #extract tab separated fields 
                fields= l.strip().split('\t')
                vcf_data['CHROM'].append(fields[0])
                vcf_data['POS'].append(int(fields[1]))
                vcf_data['ID'].append(fields[2])
                vcf_data['REF'].append(fields[3])
                vcf_data['ALT'].append(fields[4])
                vcf_data['QUAL'].append(fields[5])

                for i, col in enumerate(snp_columns):
                    vcf_data[col].append(fields[6+i].split("|"))

    df= pd.DataFrame(vcf_data)
    return df, snp_columns

def genotype_tally(snp_lst, ref_lst):

    #different genotypes 
    #g1 = homozygous for reference allele
    #g2 = heterozygous
    #g3 = homozygous for alternate allele 
    g1, g2, g3= 0, 0, 0
    
    for i in range(len(snp_lst)): 
        if i in ref_lst:

            #0|0
            if (snp_lst[i][0] == '0') and (snp_lst[i][1] == '0'):
                g1 += 1

            # 0|1 or 1|0 
            elif snp_lst[i][0] != snp_lst[i][1]: 
                g2 += 1
            
            #1|1
            elif (snp_lst[i][0] == '1') and (snp_lst[i][1] == '1'):
                g3 += 1
    
    return [g1, g2, g3]


def gwas(df, snp_columns, genotypes):

     #parse txt file 
    txt_data = {'IND': [], 'DISEASE': []}   #disease = 1 if there is a disease 
    
    txt_file= open(genotypes, "r")
    
    for l in txt_file: 
        if l.startswith('IND'): 

            fields= l.strip().split('\t')

            txt_data['IND'].append(fields[0][3:])   #only append IND number 
            txt_data['DISEASE'].append(fields[1])

    txt_file.close()


    #separate individuals (IND) into two lists; healthy and sick 
    healthy= []
    sick= []

    disease_list= txt_data['DISEASE']
    
    for i in range(len(disease_list)): 
        if disease_list[i] == '0':
           healthy.append(i)
        else:
           sick.append(i)

    print("healthy: ", healthy)
    print("sick: ", sick)

    
    for col in range((len(snp_columns))):
      
        #extract each key 
        ind_key= snp_columns[col] #ie. snp_columns[1] = IND1

        #IND num 
       # ind_num= ind_key[3:]
  
        #extract each list of chromosomes 
        snp_lst= df[ind_key]


        
        #this gives the chromosomes in list form under the COLUMN of IND1, for example 
        #snp_list[i] = list corresponding to each snp 
        #so snp_list[1] would give the corresponding ROW of chromosomes in list form (for SNP1, for example)


    #extract genotype numbers for each snp 

    #for col in range(0, 1000):
 
        for i in range(len(snp_lst)):
       

            healthy_row =[] 
            sick_row = [] 

        """
        table: [    [(SNP1)[healthy],[sick]]
                    [(SNP2)[healthy],[sick]]
                                .
                                .
                                .
                    [(SNP1000)[healthy],[sick]]
        ]
        
        """

        if i in healthy: 
            healthy_row= genotype_tally(snp_lst, healthy)
        elif i in sick:
            sick_row= genotype_tally(snp_lst, sick )
        

    #add both to snp table for that snp 
    geno_list = [healthy_row, sick_row]
    
    print(geno_list)

  
   

 

    """#print(vcf_data['IND1'][1]) gives the row, ie. snp_list[1]


        #
        each individual has 1000 SNPS
        for each SNP create a table 
        the "healthy" list is the controls, the "disease" list is the cases 
        so for EACH SNP (iterate thru SNPS) calculate G1/G2/G3 
            do this for healthy and sick lists 
        
        - healthy G1/G2/G3 = first row
        - sick G1/G2/G3 = second row 
        - append to same index in table 


        #

        #initialize snp_table 
        snp_table= np.zeros((1000, 2,3), dtype = int)
       
        #
        [ [[healthy1], [sick1]], [[healthy2], [sick2]], ... , [[healthy1000], [sick1000]] ]
        -> healthy1= [10, 23, 12]

        - outer 1000
        ie. table[0] = [[healthy1], [sick1]]

        - inner1 2
        ie. table[0][0] = [healthy1]

        - inner2 3
        ie. table[0][0][0] = 10 
    

        geno_results= genotype_tally(snp_columns, snp_lst, healthy, sick)

        for i in range(len(snp_table)):
            snp_table[i] = geno_results

    #print(snp_table)
    """



   
   
def main():

    vcf_data, snp_columns = read_vcf("gwas_population.vcf")
    gwas(vcf_data, snp_columns,  "gwas_phenotypes.txt")

if __name__ == "__main__":
    main()