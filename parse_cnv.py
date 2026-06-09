import numpy as np

def parse_cnv(file_path, variables):
    with open(file_path, 'r', encoding='utf-8', errors="ignore") as f:
        lines = f.readlines()

    data_start_index = 0
    for i, line in enumerate(lines):
        if "*END*" in line:
            data_start_index = i + 1  # Start after the *END* line
            break

    # Load data into NumPy arrays, assuming the data is space/tab-separated
    data = np.loadtxt(lines[data_start_index:], dtype=float)
    
    # Initialize dictionary to store column data
    data_dict = {}

    # Assign each variable name to the corresponding column data
    for i, var in enumerate(variables):
        if i < data.shape[1]:  # Ensure we don't exceed column count
            data_dict[var] = data[:, i]
        else:
            print(f"Warning: No column {i} for variable '{var}'")

    return data_dict, data





