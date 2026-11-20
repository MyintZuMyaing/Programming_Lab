# CW1 - Assessment
# Extract split reads from a SAM alignment file to identify intron junctions in genes
# Count supporting reads and assign junctions to genes
 
# Import standard Python libraries
import re
import sys
from collections import defaultdict
from io import StringIO

OUTPUT_FILE = "3278653.txt"

# Example CIGAR: 63M146N37M
# Regular expression for extracting CIGAR operations
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNS])")

# Example location: TGME49_chrVIII:6,793,066..6,795,596(-)
# Regular expression for extracting gene locations
GENE_PATTERN = re.compile(r"([^:]+):([\d,]+)\.\.([\d,]+)\([+-]\)")

# Extract all junctions from a CIGAR string
# Separate the CIGAR string into numbers and operations
# Calculate junction positions when it encounters N
# M and D consume reference positions, while I and S do not
def find_read_junctions(start_position, cigar):

    reference_position = start_position
    found_junctions = []

    # A read may contain more than one N
    for length_text, operation in CIGAR_PATTERN.findall(cigar):

        length = int(length_text)

        if operation == "N":

            # N indicates a skipped reference region
            junction_start = reference_position
            junction_end = junction_start + length

            found_junctions.append((junction_start, junction_end))
            reference_position = junction_end

        elif operation == "M" or operation == "D":

            # M and D move forward along the reference genome
            reference_position += length

        else:

            # Ignore I and S when calculating reference coordinates
            continue

    return found_junctions


# Process one SAM record and identify uniquely mapped split reads
# Ignores header rows, unmapped reads, and reads without an N operation
# Header lines do not describe mapped reads
# Assessment description states that NH is always the last column
def read_sam_record(line):

    # Skip SAM headers and empty lines
    if line.startswith("@") or not line.strip():
        return []

    # Split the SAM record into tab-separated columns
    columns = line.rstrip("\r\n").split("\t")

    # Check that the record has the required columns
    if len(columns) < 12:
        return []

    # Only include reads with NH:i:1 that aligned to one loaction are considered
    if columns[-1] != "NH:i:1":
        return []

    # Extract the chromosome name from the 3rd SAM column (zero-based index 2)
    chromosome = columns[2]
    
    # Extract the CIGAR string from the 6th SAM column (zero-based index 5)
    cigar = columns[5]

    # Skip unmapped reads or reads without introns.
    if chromosome == "*" or "N" not in cigar:
        return []

    # Starting alignment position from the 4th SAM column (zero-based index 3)
    start_position = int(columns[3])

    junctions = find_read_junctions(start_position, cigar)

    # Include chromosome in the junction identifier
    record_junctions = []

    for junction_start, junction_end in junctions:

        record_junctions.append((chromosome, junction_start, junction_end))

    return record_junctions


# Count the number of reads supporting each unique splice junction
# Demonstrate nested loops
def count_junctions(sam_file):
    
    # Create a dictionary with a default count of zero
    # Store read counts by splice junction, starting each new junction at zero
    supporting_reads = defaultdict(int)

    # Outer loop processes each line in the SAM file
    for line in sam_file:

        # Inner loop processes each splice junction within the alignment and adds one to its read count
        for junction in read_sam_record(line):

            supporting_reads[junction] += 1

    return supporting_reads


# Extract chromosome names and gene coordinates using a regular expression
# Remove commas from gene positions
def read_gene_location(location_text):

    # Trim surrounding whitespace and check that the entire location matches the pattern
    match = GENE_PATTERN.fullmatch(location_text.strip())

    # Return None if the location is in an unexpected format
    if match is None:
        return None

    chromosome = match.group(1)
    
    # Remove commas from gene start and end positions before converting to integers
    gene_start = int(match.group(2).replace(",", ""))
    gene_end = int(match.group(3).replace(",", ""))

    return chromosome, gene_start, gene_end


# Group splice junctions by chromosome
# Sort each group by start position and read count
def organise_junctions(junction_counts):

    # Create a dictionary containing a list for each chromosome
    # Store junctions by chromosome, creating an empty list for each new chromosome
    junctions_by_chromosome = defaultdict(list)

    # Retrieve each junction and its read count
    for junction, count in junction_counts.items():

        chromosome, junction_start, junction_end = junction

        junctions_by_chromosome[chromosome].append((junction_start, junction_end, count))

    # Sort each chromosome's junctions by start position, smallest first
    # Higher read counts first when two junctions have the same start
    # The lambda is a small function that tells .sort() how to order each junction tuple (start, end, count)
    for chromosome in junctions_by_chromosome:

        junctions_by_chromosome[chromosome].sort(key=lambda item: (item[0], -item[2]))

    return junctions_by_chromosome


# Find junctions completely within a gene on the same chromosome
def junctions_for_gene(chromosome, gene_start, gene_end, grouped):

    matching_junctions = []

    # Get junctions on the gene's chromosome, or use an empty list if it is absent
    possible_junctions = grouped.get(chromosome, [])

    for junction_start, junction_end, count in possible_junctions:

        # Include junctions fully within the gene
        # Start exactly at gene_start or end exactly at gene_end
        if gene_start <= junction_start and junction_end <= gene_end:

            matching_junctions.append((junction_start, junction_end, count))

    return matching_junctions


# Write four tab-separated columns with blank lines between genes
# Gene ID, junction start, junction end, and read count
def write_gene_results(genes_file, junction_counts, output_file):

    grouped = organise_junctions(junction_counts)

    # Skip gene-table header
    next(genes_file, None)

    wrote_a_gene = False

    for gene_line in genes_file:

        # Skip blank lines and lines containing only whitespace
        if not gene_line.strip():
            continue

        # Remove the line-ending characters from the right end of the row
        # Split the row into tab-separated fields
        columns = gene_line.rstrip("\r\n").split("\t")

        # Skip rows with fewer than three columns
        if len(columns) < 3:
            continue

        # Read the gene ID from column 1 (zero-based index 0)
        gene_id = columns[0]

        # Parse the gene location from column 3 (zero-based index 2)
        gene_location = read_gene_location(columns[2])

        # Skip rows with incorrectly formatted gene locations
        if gene_location is None:
            continue

        chromosome, gene_start, gene_end = gene_location

        matching_junctions = junctions_for_gene(chromosome, gene_start, gene_end, grouped)

        # Skip genes with no matching junctions
        if not matching_junctions:
            continue

        # Add a blank line before the next gene
        if wrote_a_gene:
            output_file.write("\n")

        wrote_a_gene = True

        # Write one output row for each junction that matches the current gene
        for junction_start, junction_end, count in matching_junctions:

            output_file.write(f"{gene_id}\t{junction_start}\t{junction_end}\t{count}\n")


# Assertion tests

# Create sample SAM records to build the test input for assertion tests
def sample_sam_record(name, chromosome, position, cigar, nh=1):

    # Create a list containing the first six SAM fields
    fields = [name, "0", chromosome, str(position), "255", cigar]

    # Add the remaining required fields and the NH tag for the alignment count
    fields += ["*", "0", "0", "*", "*", f"NH:i:{nh}"]

    # Join the fields with tabs and add a newline
    return "\t".join(fields) + "\n"


# Test CIGAR processing and intron calculations
# Check if junction coordinates are calculated correctly for S, D, I, and Ns 
def test_cigar_extraction():

    # Soft clipping should not alter junction coordinates
    assert find_read_junctions(100, "5S10M5N10M") == [(110, 115)]

    # Deletions should move the reference position
    assert find_read_junctions(100, "8M2D5N10M") == [(110, 115)]

    # Insertions should not move the reference position
    assert find_read_junctions(100, "10M2I5N10M") == [(110, 115)]

    # A read can contain multiple junctions
    assert find_read_junctions(100, "10M5N10M6N10M") == [(110, 115), (125, 131)]

    # Check the coordinate example from the assignment
    assert find_read_junctions(3000000, "63M146N37M") == [(3000063, 3000209)]


# Test read counts and chromosome-specific junctions
# Check if reads supporting the same junction are counted together
# Check if the SAM header and multimapping read are skipped
# Check if identical coordinates on different chromosomes are counted separately
def test_sam_counting():

    # StringIO creates a file-like object in memory
    sam_text = "@HD\tVN:1.6\n"

    sam_text += sample_sam_record("a", "chr8", 100, "10M5N10M")
    sam_text += sample_sam_record("b", "chr8", 100, "10M2I5N10M")
    sam_text += sample_sam_record("c", "chr8", 100, "10M5N10M6N10M")

    # Multimapping read should be ignored
    sam_text += sample_sam_record("d", "chr8", 100, "10M5N10M", 2)

    # A junction on another chromosome is counted separately
    sam_text += sample_sam_record("e", "chr9", 100, "10M5N10M")

    counts = count_junctions(StringIO(sam_text))

    expected_counts = {("chr8", 110, 115): 3, ("chr8", 125, 131): 1, ("chr9", 110, 115): 1}

    # Check actual results against expected results
    assert dict(counts) == expected_counts, "Incorrect junction counts"


# Test gene boundaries, chromosome matching and output format
# Check if gene locations with commas are parsed
# Check if incorrectly formatted locations return None
# Check if junctions are matched using chromosome and gene boundaries
# Check the exact tabs and blank lines in the output
def test_gene_output():

    assert read_gene_location("chr8:1,000..2,000(-)") == ("chr8", 1000, 2000)

    assert read_gene_location("incorrect location") is None

    genes_text = (
        "Gene\tTranscript\tLocation\n"
        "geneA\tt1\tchr8:100..300(+)\n"
        "geneB\tt2\tchr8:400..500(-)\n"
        "geneC\tt3\tchr9:100..300(+)\n")

    junction_counts = {
        ("chr8", 150, 200): 3,
        ("chr8", 150, 210): 1,
        ("chr8", 450, 500): 2,
        ("chr9", 150, 200): 7}

    # Create an empty in-memory text file to capture the results
    output = StringIO()

    # Read the sample gene table and write the results to the in-memory file
    write_gene_results(StringIO(genes_text), junction_counts, output)

    # Include all expected junctions and blank lines between genes
    expected_output = (
        "geneA\t150\t200\t3\n"
        "geneA\t150\t210\t1\n"
        "\n"
        "geneB\t450\t500\t2\n"
        "\n"
        "geneC\t150\t200\t7\n")

    # Verify that the actual output exactly matches the expected output
    assert output.getvalue() == expected_output, "Gene matching or output formatting is incorrect"


# Run assertion tests before processing the real files
def run_tests():

    test_cigar_extraction()
    test_sam_counting()
    test_gene_output()

    print("All assertion tests passed.")


# Main function controls input, testing, processing and error handling
# It coordinates the workflow from checking input arguments to writing the results
def main():

    # Require exactly two input filenames on the command line
    if len(sys.argv) != 3:

        sys.exit("Usage: python myScript.py mySamFile.sam myInputTable.txt")

    # Get input filenames from command line
    sam_filename = sys.argv[1]
    genes_filename = sys.argv[2]
        
    # Run assertion tests
    run_tests()

    # Try-except block
    
    # Process the input files and write the results
    try:

        # Open the SAM file for reading alignments and counting junctions
        with open(sam_filename, "r") as sam_file:
            junction_counts = count_junctions(sam_file)

        # Open the gene table for reading and output file for writing
        with open(genes_filename, "r") as genes_file, \
            open(OUTPUT_FILE, "w") as output_file:

            # Match junctions to genes and write the results
            write_gene_results(genes_file, junction_counts, output_file)

    # Catch either of these exceptions from the try block
    # A missing file or a permission problem can raise an OSError
    # A failed integer conversion can raise a ValueError
    except (OSError, ValueError) as error:

        # Stop the program with an error exit status
        sys.exit(f"Error while processing files: {error}")

    print(f"Results written to {OUTPUT_FILE}")


# Start the program
# Run main() when the script is executed directly
if __name__ == "__main__":
    main()