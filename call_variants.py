from collections import defaultdict

ref_genome =open('TreeToReads/example/mini_ref.fasta', 'r')
lines = ref_genome.readlines()
var_sites_matrix = open('TreeToReads/example_out/var_site_matrix', 'r')

    #expect (number between 1 & |sample| or name ) and (number in mutations locations list)
def get_base(sample, location):
	base = ''
	#if "original", then we take the base from the reference genome
	if sample == "original":
		lineNum = 1+location / 70
		characterNum = location % 70
		line = lines[lineNum]
		base = line[characterNum]
	#otherwise, we take the base from the varsitesmatrix file (output of treetoreads)
	else:
		for line in var_sites_matrix:
			k = ""
			ln = line.split()
			if sample == ln[0] and str(location) == ln[2]:
				base = ln[1]
	return base 

def call_variant(ref_base, sample_doubt, sample_list, mutation_doubt, mutation_list):
    call = ref_base
    non_ref = ""
    snp_matrix = []

    cutoff_dist = 1000 #in kilobases but 1000 bases chosen for simulated data. User defined
    proximal_mutations = []
    other_samples = []
    for sample in sample_list:
        if sample != sample_doubt:
            other_samples.append(sample)
    for i in range (len(mutation_list)):
        mutation = mutation_list[i] # a no. ie mutation location
        dist = mutation - mutation_doubt
        if (dist < 0):
            dist = dist * -1
        
        if (dist <= cut_off and dist != 0 ):
            proximal_mutations.append(mutation)
    for i in range (len(proximal_mutations)):
        mutation = proximal_mutations[i]
        snp_matrix.append([mutation]) #append an list meant to represent SNP at this location (across the samples). First element in the list is the position of the mutation between 0 and length of the genome. The element at position n (n between 1 and no. of samples) from the second onward will be either 0 or 1 depending on whether the sample has a non variant or variant allele at this location. 
        for j in range (len(other_samples)):
            sample = other_samples[j] # a string name of sample
            sample_base = get_base(sample, mutation)#base in this sample at this mutation location
            if sample_base == ref_base:
                snp_matrix[i][j+1] = 0 #if non variant

            else:
                if non_ref = "":
                    non_ref = sample_base #assign non ref char value
                snp_matrix[i][j+1] = 1 #if variant



    basis_snp = 0 #to store the SNP that the in doubt mutation is most closely related to (r-sq linkage disequilibrium). Our variant call should infer base based on this basis_snp

    snp_doubt = [mutation_doubt] #create list representation of in doubt SNP

    for i in range (len(other_samples)):#creating representation of this snp in doubt
        sample = other_samples[i]
        sample_base = get_base(sample, mutation_doubt)
        if sample_base == ref_base:
            snp_doubt[i+1] = 0 #if non variant

        else:
            snp_matrix[i+1] = 1#if variant


    r_sq_cutoff = .96 #user defined and requires statistical deliberation
    found_basis = False
    for snp in snp_matrix:
        r_sq = get_r_square(snp_doubt, snp)
        if r_sq > r_sq_cutoff:
            basis_snp = snp[0]#a snp to base the call on
            found_basis = True
            break

    if (found_basis is not True):#searches through second order SNP combinations
        l = len(snp_matrix)
        curr_snp_pos = l #to store pos at which new snps have to be instroduced to snp_matrix
        for i in range (l-1):
            for j in range (i+1, l): #loop combination searches all combinations of 2 snps to create new constructed snps
                snp_matrix.append([(snp_matrix[i][0], snp_matrix[j][0])])#first of 2 new spns added 
                for k in range (1, len(other_samples)):
                    if  (snp_matrix[i][k] and (not snp_matrix[j][k])):#suggests first snp happened above a fork between the second and mutation_doubt
                        snp_matrix[curr_snp_pos][k] = 1
                    else:
                        snp_matrix[curr_snp_pos][k] = 0

                curr_snp_pos = curr_snp_pos + 1#increment pos 

                snp_matrix.append([(snp_matrix[j][0], snp_matrix[i][0])])#second of 2 new snps
                for k in range (1, len(other_samples)):
                    if  (snp_matrix[j][k] and (not snp_matrix[i][k])):#suggests second snp happened above a fork seperating the first and mutation doubt
                        snp_matrix[curr_snp_pos][k] = 1
                    else:
                        snp_matrix[curr_snp_pos][k] = 0
                    
                curr_snp_pos = curr_snp_pos + 1 #increment pos



        basis_snp_tup = ()#the basis for calling is on a tuple of snps
        for i in range (l, len(snp_matrix)):

            r_sq = get_r_square(snp_doubt, snp_matrix[i])
            if r_sq > r_sq_cutoff:
                basis_snp_tup = snp_matrix[i][0]#second order cosntructed  snp provides strong basis for a call
                found_basis = True
                break

        if found_basis is True:
            guess = get_base(sample_doubt, basis_snp_tup[0]) #use tuple discovered as basis for call

    else:

        guess = get_base(sample_doubt, basis_snp) #use snp discovered as basis for call


    return guess #return call
            

def get_r_square (snp1, snp2):

    #TODO r_sq implementation
    #Functionality yet to be implemented
    #we suggest R's Linkage Disequilibrium package as a way of getting LD values 
    #given snp1 and snp2 are headed by the mutation location they represent followed by  binary strings of length no. of samples. Each position i is 1 or 0 depending on whether the i th sample has a variant or non variant allele at the mutation location. For eg, if samples 1, 3, 5, and 8 have non ref allele and 2,4,6,7 do not at location 789 then the snp will look like a list [789,1,0,1,0,1,0,0,1]
    #returns an LD value between 0 and 1



if __name__ == '__main__':

        
        #list of mutation locations 
        #list of samples 
        #iterations = 100
        #calls_correct value to store number of times the algorithm correctly calls the base for given sample and location
        #for i in range(iterations): loop to arbitrary amount for testing
                #choose arbitrary mutation x location from list
                #pick one sample at random s
                #ground truth = get_base(sample, location)
                #ref_base = get base("original", location)
                #if call_variant(ref_base, s, sample_list, x, mutation_list) == ground_truth:
                        #calls_correct += 1

       #fraction_correct = calls_correct/iterations
       #percentage_correct = fraction_correct * 100
       #print "algorithm estimated call correctly" + str(percentage_correct) + "% of the iterations"

