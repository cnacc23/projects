from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np 

#function to find the number of homo and heterozygous genotypes in an snp 
def tally_genotypes(snp_list):
    g0, g1, g2= 0, 0, 0

    for i in range(len(snp_list)):
       
       #homozygous cases
        if snp_list[i][0] == snp_list[i][1]:

            #homozygous for ref allele (0|0)
            if snp_list[i][0] == '0':
                g0 += 1

            #homozygous for alt allele (1|1)
            elif snp_list[i][0] == '1':
                g2 += 1

        #heterozygous case (0|1 or 1|0)
        else:
            g1 += 1

    return[g0, g1, g2]

#function that calculates the conditional probability of two events
##CHECK 
def cond_pbt(a, b, total):
    return ((a+b)/total) / (b/total)


#function that filters a list based on indexes 
def filter_list(condition, raw):

    #this list will contain the filtered contents 
    filtered= []

    for i in range(len(condition)):
        index = condition[i]
        filtered.append(raw[index])

    return filtered



def read_vcf(population, genotypes ):


    #parse txt file 
    txt_data = {'IND': [], 'DISEASE': []}   #disease = 1 if there is a disease 

    txt_file= open(genotypes, "r")

    for l in txt_file: 
        if l.startswith('IND'): 

            fields = l.strip().split('\t')

            txt_data['IND'].append(fields[0][3:])   #only append IND number 
            txt_data['DISEASE'].append(fields[1])

    txt_file.close()


    #separate individuals (IND) into two lists; healthy and sick 
    healthy = []
    sick = []

    disease_list = txt_data['DISEASE']

    for i in range(len(disease_list)): 
        if disease_list[i] == '0':
            healthy.append(i)
        else:
            sick.append(i)


    #parse vcf file  

    #create dictionary to store each value in pertaining list 
    vcf_data = {'CHROM': [], 'POS': [], 'ID': [], 'REF': [], 'ALT': [], 'QUAL': []}

    with open(population, 'r') as vcf_file:
    
        for l in vcf_file:

            #extract SNP columns from header
            if l.startswith('#'): 
                header_fields= l.strip().split('\t')

                #take columns after all other fields in dictionary 
                ind_columns= header_fields[6:]
                

                #initialize to empty list 
                for col in ind_columns: 
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

                for i, col in enumerate(ind_columns):
                    vcf_data[col].append(fields[6+i].split("|"))
            
    

   #extract alleles for each individiual  
    
    gwas_table = []
    ind_alleles = []
    
    for col in range(len(ind_columns)):
        
        ind_key= ind_columns[col] # gives 'IND1' if col == 1

        #append each individiual's alleles to a matri
        ind_alleles.append(vcf_data[ind_key])

        #contains the tally numbers for each genotype 
    
        #this will give alleles for each snp!!
        # i represents snp num 
        # j represents individual
       
    

    #ind_alleles contains 10000 snps, each 1000 nums long'
    for i in range(10000):    #snp number 
        curr_snp = [] 
        curr_snp_healthy = []
        curr_snp_sick = []

        for j in range(len(ind_alleles)):   #individual 
    
            if j <= 1000:
                #sort each snp based on healthy and sick individuals 
                if j in healthy: 
                    curr_snp_healthy.append(ind_alleles[i][j])
                else:
                    curr_snp_sick.append(ind_alleles[i][j])

            #sum different genotypes for each snp based on if healthy/sick
            curr_snp.append(tally_genotypes(curr_snp_healthy))
            curr_snp.append(tally_genotypes(curr_snp_sick))
       
            gwas_table.append(curr_snp)

    print(len(gwas_table))
    return gwas_table
   

#calculate the p-values for each snp and disease
def p_values(table):

    uncorrected, corrected, num_tests = 0, 0, 0
    uncorrected_pvals, corrected_pvals, significant_results, odd_ratio_het_list, odd_ratio_alt_list = [], [], [], [], []

    #bonferroni correction; len table = num_tests
    bonferroni = 0.05/num_tests

   #num_tests * p < 0.05

    num_tests = 0
    for i in range(len(table)):
        for j in range(2):

            #get rid of 0 columns 
            if table[i][j][0] == 0 or table[i][j][1] == 0:
                num_tests += 1
        
        #calculate pvalue for each snp using chi2_contingecy 
        pvals = chi2_contingency(table[i])[1]
        uncorrected_pvals.append(pvals)
        

    #extract number of uncorrected p-values < 0.05
    for i in range(len(uncorrected_pvals)):
        if uncorrected_pvals[i] < 0.05: uncorrected += 1
        
   
    #now correct p-values
    #bonferroni correction; len(table) = num tests 
    corrected_pvals = multipletests(uncorrected_pvals, method= "bonferroni")[1]
    
    #extract significant corrected pvals    
    sig_corrected = []

    for i in range(len(uncorrected_pvals)):
        if corrected_pvals[i] < bonferroni: 
            significant_results.append(i)
            sig_corrected.append(corrected_pvals[i])
            corrected += 1

    # extract significant uncorrected and corrected pvals 
    sig_uncorrected = filter_list(significant_results, uncorrected_pvals)
    
    #odds_ratio 

    for i in range(len(table)):

        #num sick individuals per snp 
        sick = sum(table[i][1])


        #for each snp, go through healthy and sick lists of genotypes
        for j in range(len(table)):

            #homo for reference allele 
            hom_ref = (table[i][0][0] + table[i][1][0]) 

            #hetero
            het = (table[i][0][1] + table[i][1][1])                                                                                         

            #homo for alt allele
            hom_alt = (table[i][0][2] + table[i][1][2]) 

            #calculate both pbts
            odds_ratio_het = cond_pbt(sick, het, 10000) / cond_pbt(sick, het, 10000)
            odds_ratio_alt = cond_pbt(sick, hom_alt, 10000) / cond_pbt(sick, hom_ref, 10000)
    
    
        #append pbts to lists 
        odd_ratio_het_list.append(odds_ratio_het)
        odd_ratio_alt_list.append(odds_ratio_alt)


    #create table of significant snp values 
    snp_data = {"SNP ID": significant_results, 
                "UNCORRECTED P-VALUE": sig_uncorrected , 
                "CORRECTED P-VALUE": sig_corrected, 
                "DISEASE ODDS RATIO FOR HETEROZYGOUS INDIVIDUALS": filter_list(significant_results, odd_ratio_het_list), 
                "DIESEASE ODDS RATIO FOR HOMOZYGOUS INDIVIDUALS": filter_list(significant_results, odd_ratio_alt_list)}
    
    df = pd.DataFrame(snp_data)
    print(df)
    print("uncorrected pvalues: ", uncorrected, " corrected pvalues: ", corrected)

   
def main():

    gwas_table= read_vcf("gwas_population.vcf", "gwas_phenotypes.txt")
    #p_values(gwas_table)
    
if __name__ == "__main__":
    main()
