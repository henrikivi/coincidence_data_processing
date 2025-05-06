"""Functions to load coincidence data for EIEI
"""
import os
import csv

def load_times(filename, order, headers):
    """
    Function to load times from a pairs file
    file should be two columns of t1 and t2 with white space inbetween
    if header lines present set number using headers variable
    inputs
    filename - path to file
    order - if 0 then t2 > t1 if 1 then t1 > t2
    headers - how many header lines in file

    Outputs
    t1, t2 - lists of the times

    """
    t1 = []
    t2 = []
    with open(filename, "r") as infile:
        for _ in range(headers):
            next(infile)
        for line in infile:
            # read in
            temp = line.split()
            t1.append(int(temp[0]))
            t2.append(int(temp[1]))
    if order == 1:
        temp = t2
        t2 = t1
        t1 = temp

    return t1, t2


def load_eiei_data(filestub, t1long=False):
    """Function load results from a EIEI data run
    reads in the singles, pairs and triples.
    Is written for EIEI but would be easy to generalise to other coincidence measurements

    Parameters
    ----------
    filestub : str
        the file name with out the file type at end

    t1long : bool
        is data sorted as t1 longest time or t2 longest. Default is False - t1 shortest
    Returns
    -------
    dict
        A dictionary containing the results, separated as singles, pairs, triples and a dictionary with experimental information
    """
    singles = []
    pairs = []
    triples = []
    dat_info = {}
    with open(filestub + ".mas", "r") as infile:
        # read first line
        txt = infile.readline()
        if txt == "Version 2.0\n":
            # This is electron impact file type
            dat_info["mactype"] = 1
            # next line is file name
            dat_info["file"] = infile.readline()
            # Read next few lines as headers
            for __ in range(27):
                temp = infile.readline().split(",")
                dat_info[temp[0]] = temp[1]
            # set up singles spectrum
            for line in infile:
                singles.append(int(line))
    # Now get pairs
    with open(filestub + ".prs", "r") as infile:
        # First line is number of pairs - this is knonwn
        infile.readline()
        for line in infile:
            temp = line.split()
            prs = []
            for ch in temp:
                prs.append(int(ch))
            pairs.append(prs)
    # Now get triples
    with open(filestub + ".trp", "r") as infile:
        # First line is number of triples - this is knonwn
        infile.readline()
        for line in infile:
            temp = line.split()
            trp = []
            for ch in temp:
                trp.append(int(ch))
            triples.append(trp)
    # Generate the time axis for future use
    time_axis = range(int(dat_info["First time:"]), int(dat_info["Last time:"]))
    results = {
        "singles": singles,
        "pairs": pairs,
        "triples": triples,
        "info": dat_info,
        "time_axis": time_axis,
    }
    return results

def convert_value(value):
    """
    Convert an input string to int, float, or bool if possible. Used in the parameter loading function.
    """
    value = value.strip()  # Remove extra spaces
    
    #if value == "":
    #    return value
    if value.lower() == "true":  
        return True
    elif value.lower() == "false":  
        return False
    try:
        if "." in value:  
            return float(value)  # Convert to float if it has a decimal point
        return int(value)  # Convert to int otherwise
    except ValueError:
        return value  # Keep as string if conversion fails

def load_params(file_name):
    """
    Loads the peak parameters as lists for automated analysis from a csv file.

    Takes argument file_name (str), which is the name of the csv file (without extension)
    in the same directory as the script to be executed

    Csv file needs to be set up as three separate csv tables: Singles, pairs and triples.
    Below is a schematic of how the csv file should look:

    #Singles
    #Ion,m/z,LHS,RHS,BLHS,BRHS
    [SINGLES_PARAMS_HERE]
    #Pairs
    #Frag. 1,Frag. 2,m/z(1),m/z(2),x1,y1,x2,y2,Box Gradient,Box shape
    [PAIRS_PARAMS_HERE]
    #Triples
    #Frag. 1,Frag. 2,Frag. 3,m/z(1),m/z(2),m/z(3),t1,t2,x1,y1,x2,y2,Box Gradient,Box shape
    [TRIPLES_PARAMS_HERE]

    The two hashed lines before each parameter set are important, 
    the code will not work if they are removed!
    
    Input parameters in csv file:
    Ion & Frag. X: str
        Identity of the fragment ion (e.g. H+)
    m/z(X): int or float
        mass-to-charge ratio of the fragment (X)
    LHS & RHS: int
        Left and right bounds of the peak
    BLHS & BRHS: int
        Left and right bounds for background calculation
    x1 & y1 & x2 & y2: int
        coordinates of the box for pairs to be counted
    t1 & t2: int
        time range for triples selection
    Box gradient: int or float
        gradient of the aforementioned box
    Box shape: boolean
        if True, flat-topped box. If False, regular box
    DTC: boolean
        Dead time correction. If True, correction is applied
        and a csv file is created with a time of flight difference plot with 
        the fragment m/z:s as file names

    Returns
    ------------
    sections[0]: list
        List of singles parameters as a 2d list
    sections[1]: list
        List of pairs parameters as a 2d list
    sections[2]: list
        List of triples parameters as a 2d list
    """
    sections = [[], [], [], []]  # Four lists for Singles, Pairs, Triples, shift calibration
    section_index = -1  # Start before the first section
    skip_next = False  # Flag to skip the next row (headers)

    with open(f"{file_name}.csv", mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)

        for row in reader:
            # Remove leading/trailing spaces, do not ignore empty columns
            row = [r.strip() if r.strip() else "" for r in row]
            
            # Remove trailing empty columns
            while row and row[-1] == "":
                row.pop()
            
            # Skip empty rows (entirely blank)
            if not row:
                continue

            # Skips the first row after a section header (the actual column headers)
            if skip_next:
                skip_next = False  # Reset flag
                continue  # Skip the header row

            # If the row starts with "#", it's a section header
            if row[0].startswith("#"):
                section_index += 1  # Move to the next section
                skip_next = True  # Mark that the next row is a header
                continue  # Skip section header row

            # Convert values before storing
            converted_row = [convert_value(value) for value in row]

            # Append data rows to the correct section (if valid)
            if 0 <= section_index < len(sections):
                sections[section_index].append(converted_row)

    return sections[0], sections[1], sections[2], sections[3]  # Singles, Pairs, Triples, Calibration

def save_results(single_results, pairs_results, triples_results, filename="results.csv", folder="."):
    """Write the results to a CSV file in a specified folder.

    Args:
        single_results: Data for single results. (generated by analyse_files.py-script)
        pairs_results: Data for pairs results. (generated by analyse_files.py-script)
        triples_results: Data for triples results. (generated by analyse_files.py-script)
        filename (str): Name of the output CSV file.
        folder (str): Directory where the output CSV file will be saved.
    """
    # Ensure the folder exists
    os.makedirs(folder, exist_ok=True)

    # Create the full path for the CSV file
    file_path = os.path.join(folder, filename)

    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file)

        # Write single results
        writer.writerow(["Singles:"])  # Header for single results
        writer.writerow(["m/z","Counts","Background"])
        for result in single_results:
            writer.writerow(result)
        writer.writerow([])  # Blank line for separation

        # Write pairs results
        writer.writerow(["Pairs:"])  # Header for pairs results
        writer.writerow(["m/z(1)","m/z(2)","Counts","Gradient","\u0394Gradient"])
        for result in pairs_results:
            writer.writerow(result)
        writer.writerow([])  # Blank line for separation

        # Write triples results
        writer.writerow(["Triples:"])  # Header for triples results
        writer.writerow(["m/z(1)","m/z(2)","m/z(3)","Counts","Gradient","\u0394Gradient"])
        for result in triples_results:
            writer.writerow(result)

    print(f"Results saved to {file_path}")