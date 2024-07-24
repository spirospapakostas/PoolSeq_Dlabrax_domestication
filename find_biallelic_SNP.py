#################################################

# ============ FIND BIALLELIC LOCI ============ #
# ================== STEP1/2 ================== #

#################################################
# set to zero the counts if, per 4 nucleotides (one sample), counts are less or equal to the following thresholds
# thresholds: <1% for total_counts>100 and <1.25% for total_counts>100 (total counts == counts of all 4 nucleotides)
# ============================================= #

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

TOT_POS = range(748392)
A_POS = range(1, 97, 4)
THRESHOLD = 1
THRESHOLD_LOW = 1.25
ZERO_DIV_POSITIONS = []


def set_zeros(SUM, count, nucl_type):  # nucl_type is "a", "c", "g", or "t"
    
        if SUM <= 100 and ((count/SUM * 100) < THRESHOLD_LOW) and nucl_type=="a":
            return 0
        elif SUM > 100 and ((count/SUM * 100) < THRESHOLD) and nucl_type=="a":
            return 0
        elif SUM <= 100 and ((count/SUM * 100) < THRESHOLD_LOW) and nucl_type=="c":
            return 0
        elif SUM > 100 and ((count/SUM * 100) < THRESHOLD) and nucl_type=="c":
            return 0
        elif SUM <= 100 and ((count/SUM * 100) < THRESHOLD_LOW) and nucl_type=="g":
            return 0
        elif SUM > 100 and ((count/SUM * 100) < THRESHOLD) and nucl_type=="g":
            return 0
        elif SUM <= 100 and ((count/SUM * 100) < THRESHOLD_LOW) and nucl_type=="t":
            return 0
        elif SUM > 100 and ((count/SUM * 100) < THRESHOLD) and nucl_type=="t":
            return 0
        else:
            return count

        
def alleles_to_zero(filename):
    for position in TOT_POS:
        print(position, end=" ")
        for A in A_POS:
            try:
                a = filename.iloc[position, A]
                c = filename.iloc[position, A+1]
                g = filename.iloc[position, A+2]
                t = filename.iloc[position, A+3]
                SUM = a+c+g+t
                filename.iloc[position, A] = set_zeros(SUM, a, "a")
                filename.iloc[position, A+1] = set_zeros(SUM, c, "c")
                filename.iloc[position, A+2] = set_zeros(SUM, g, "g")
                filename.iloc[position, A+3] = set_zeros(SUM, t, "t")
                SUM = 0
            except ZeroDivisionError:
                ZERO_DIV_POSITIONS.append(position)
                continue

alleles_to_zero(genotypes)

###################################################################################
#===============================STEP 2/2==========================================#
# find biallelic loci
# distinguish between monomorphic SNPs (everywhere only the same nucleotide has non-zero counts)
# and distinguish between SNPs with >=3 alleles
# if zero_counts > 2 across pops (then the SNP is monomoprhic)
# if zero_counts < 2 across pops (then the SNP is triallelic or multiallelic
#=================================================================================#

TOT_POS = range(748392)
A_POS = range(1, 97, 4)
DROPPED = []


def filt_mono_SNP(filename):
    for position in TOT_POS:
        print(position, end=" ")
        a_list = []
        c_list = []
        g_list = []
        t_list = []
        nucl_sums = []
        for A in A_POS:
            a = filename.iloc[position, A]
            c = filename.iloc[position, A+1]
            g = filename.iloc[position, A+2]
            t = filename.iloc[position, A+3]
            a_list.append(a)
            c_list.append(c)
            g_list.append(g)
            t_list.append(t)
        nucl_sums.append(sum(a_list))
        nucl_sums.append(sum(c_list))
        nucl_sums.append(sum(g_list))
        nucl_sums.append(sum(t_list))

# =========== code section for keeping only biallelic =========== #
         if nucl_sums.count(0) > 2 or nucl_sums.count(0) < 2:
             DROPPED.append(position)
         else:
             continue
     return filename.drop(index=DROPPED)

 result_biallelic = filt_mono_SNP(genotypes) 

# =========== code section for keeping the monomorphic =========== #
         if nucl_sums.count(0) > 2:
             DROPPED.append(position)
         else:
             continue
     return filename.drop(index=filename.index.difference(DROPPED))

 result_monomorphic = filt_mono_SNP(genotypes)

# =========== code section for keeping the triallelic =========== #
        if nucl_sums.count(0) < 2:
            DROPPED.append(position)
        else:
            continue
    return filename.drop(index=filename.index.difference(DROPPED))
            
result_triallelic = filt_mono_SNP(genotypes)
