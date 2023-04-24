# Obtain the DNA sequences that code for the amino acid sequence and do NOT have any stop codons DNA sequences

# define the genetic code as a dictionary
genetic_code = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'],
    'D': ['GAT', 'GAC'],
    'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'],
    'M': ['ATG'],
    'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    '*': ['TAA', 'TAG', 'TGA']
}


dp_array = []

# define the amino acid sequence
aa_sequence = 'MRKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYKRPAANDENYAASV'
# aa_sequence = 'MRKGEELFTG'

if len(aa_sequence) == 1:
    # Check if it is a stop codon
    if aa_sequence[0] == '*':
        dp_array =  []
    else:
        dp_array = genetic_code[aa_sequence[0]]

valid = True
# dp_array only needs to store the previous DNA sequences
for idx in range(len(aa_sequence) - 1):
    if not dp_array:
        # Get the DNA sequences that code for the previous amino acid
        prev_dna = genetic_code[aa_sequence[0]]
        # Get the DNA sequences that code for the current amino acid
        curr_dna = genetic_code[aa_sequence[1]]
    else:
        prev_dna = dp_array[-1]
        curr_dna = genetic_code[aa_sequence[idx]]

    curr_dp = set()
    # Do sliding window on all the combinations of the previous and current DNA sequences to see if they contain the stop codon
    for prev in prev_dna:
        for curr in curr_dna:
            # If they don't overlap, then we cannot use this combination
            if dp_array and prev[-3:] != curr:
                continue
            # Get rid of the overlap if there is any
            total_seq = prev[-3:] + curr
            # If the total sequence contains the stop codon, then we cannot use this combination
            if 'TAA' in total_seq or 'TAG' in total_seq or 'TGA' in total_seq:
                continue
            # If the total sequence does not contain the stop codon, then we can use this combination
            # Add the current DNA sequence to the curr_dp array
            curr_dp.add(total_seq)
    if not curr_dp:
        valid = False
        break
    dp_array.append(curr_dp)

ex_array = dp_array
number_of_valid = 50 # number of valid DNA sequences we want to find
# Combine two strings if their last 3 and first 3 characters match
def combine_strings(s1, s2):
    if s1[-3:] == s2[:3]:
        return s1 + s2[3:]
    else:
        return None
# Repeatedly combine strings until the array is empty
while len(ex_array) > 1:
    combined = set()
    # Combine strings from the first two sets in the array
    for s1 in ex_array[0]:
        for s2 in ex_array[1]:
            new_str = combine_strings(s1, s2)
            if new_str and new_str not in combined:
                combined.add(new_str)
            if len(combined) > number_of_valid:
                break
    # Replace the first two sets with the combined set
    ex_array = [combined] + ex_array[2:]

# Return the final set of combined strings
result = list(ex_array[0])