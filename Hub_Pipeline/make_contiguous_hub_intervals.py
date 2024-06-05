#### --------------------------------------------------------------------------------------------------------------------------//
#### Execution: python make_contiguous_hub_intervals.py Control_reconnected_EDGE.tsv Control_hub_by_interval.txt
#### Overview: This is a python helper script to the "hub_pipeline.sh" script that finds the largest genomic interval of a given hub
# within a list based on its enh/pro coordinates to eventually create a new list of hubs with contiguous genomic intervals
#
#### Args: .tsv file of hubs in a given condition annotated with interaction count, enh/pro count, and list of enh/pro intervals (chrom_start_end); output filename
#
#### Final Output: .tsv file that contains a datatable of all hubs defined by their largest contiguous genomic interval
# with interaction and enh/pro counts
#### -------------------------------------------------------------------------------------------------------------------------//

#### Import necessary python libraries
import sys
import linecache

#### Open input file and count number of lines/cliques
inputFile = open(sys.argv[1], "r") 
content = inputFile.read()
count_one_list = content.split("\n")
count_one = 0; count_two = 0

for line in count_one_list:
    count_one+=1

#### Move to beginning of input file and create file to hold output with header
inputFile.seek(0)
output_name = sys.argv[2]
outfile = open(output_name, 'w')
outfile.write("chr" + '\t' + "start" + '\t' + "end" + '\t' + "hub_id" + '\t' + "interaction_count" + '\t' + "enh/pro_count" + "\n")

#### Generate contiguous genomic interval for hubs in given condition by reading each enh/pro coordinate 
# from the list of all enh/pro for a given hub and searching for the most upstream and downstream enh/pros within the list
for i in range(1, count_one):
    # Read specific line, split it into fields, and store relevant data
    line = linecache.getline(sys.argv[1], i)
    fields = line.split('\t')
    name = fields[0]; edge = fields[1]; node = fields[2]; items = fields[3]
    intervals = items.strip("\n").split(' ')
    
    # Initialize loop variables with dummy values
    start = 10000000000000; end = 0; j = 1

    # Search for largest interval of hub by iterating through enh/pro coordinates 
    for coordinate in intervals:
        region = coordinate.split('_'); chrom = region[0]; starti = int(region[1]); endi = int(region[2])

        if starti < start:
            start = starti
        
        if endi > end:
            end = endi
        
        j+= 1
    
    # Write the contiguous genomic interval to the output file
    outfile.write(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + name + '\t' + str(edge) + '\t' + str(node) + "\n")

# Close remaining files
inputFile.close()
outfile.close()

