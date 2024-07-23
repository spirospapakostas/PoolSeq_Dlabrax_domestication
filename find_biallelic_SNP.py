#################################################

# ============ FIND BIALLELIC LOCI ============ #
# ================== STEP1/2 ================== #

#################################################
# set to zero the counts if, per 4, are less or equal to the following thresholds
# set zero counts first to split the process in a simplified way
# <1% for depth>100 and <1.25% for depth>100
# 1.25% is allowing more than 1 count at depth>100
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