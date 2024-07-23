import random
import numpy as np
from collections import Counter

# Global variables
SIMULATIONS = 100000
SAMPLE_POOL = 25
DEPTH = 80

# Initialize a list to store the number of representatives in each simulation
num_representatives = []

# Run the simulations
for i in range(SIMULATIONS):
    sample = [random.choice(range(1, SAMPLE_POOL)) for i in range(DEPTH)]
    count = len(Counter(sample))
    num_representatives.append(count)

# Calculate the average number of representatives
average_representatives = np.mean(num_representatives)

# Calculate the lower 95% confidence interval
confidence_interval = np.percentile(num_representatives, 2.5)

# Print the results
print("Average number of representatives: {:.2f}".format(average_representatives))
print("Lower 95% confidence interval: {:.2f}".format(confidence_interval))