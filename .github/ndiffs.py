import re
import numpy as np

def extract_numbers_from_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", content)
    return list(map(float, numbers))

def compute_least_squares_difference(file1, file2):
    numbers1 = extract_numbers_from_file(file1)
    numbers2 = extract_numbers_from_file(file2)

    if len(numbers1) != len(numbers2):
        print ('num values differ. Using min of ',len(numbers1), len(numbers2) )
        #nn = np.min( [len(numbers1), len(numbers2)] )
        #raise ValueError("Files must contain the same number of numerical values")
    nn = np.min( [len(numbers1), len(numbers2)] )
    # Compute the least squares difference
    differences = np.array(numbers1[:nn]) - np.array(numbers2[:nn])
    least_squares_difference = np.sqrt(np.sum(differences**2))/float(nn)
    return least_squares_difference

file1 = 'out15' #'path/to/your/first_file.txt'
file2 = 'out15.ref' #'path/to/your/second_file.txt'

result = compute_least_squares_difference(file1, file2)
print(f"Least Squares Difference: {result}")
