import re

def read_cnv_header(file_path):
    latitude = None
    longitude = None
    cast_num=None
    variables = []

    with open(file_path, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            
            # Extract cast number (marked with '* ** cast number')
            if "* ** cast number" in line:
                match = re.search(r"[-+]?\d+", line)  # Match floats without decimal points
                if match:
                    cast_num = float(match.group())
                    
            # Extract latitude (marked with '* ** lat')
            if "* ** lat" in line and "decimal degrees" in line:
                match = re.search(r"[-+]?\d*\.\d+", line)  # Match floats with decimal points
                if match:
                    latitude = float(match.group())

            # Extract longitude (marked with '* ** lon')
            if "* ** lon" in line and "decimal degrees" in line:
                match = re.search(r"[-+]?\d*\.\d+", line)  # Match floats with decimal points
                if match:
                    longitude = float(match.group())

            # Extract variable names (arked with "# name")
            if "# name " in line:
                parts = line.strip().split("=")
                if len(parts) > 1:
                    variables.append(parts[1].strip())
                    

    return cast_num, latitude, longitude, variables

