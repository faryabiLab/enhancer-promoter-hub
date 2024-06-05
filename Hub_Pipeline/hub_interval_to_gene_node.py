#### --------------------------------------------------------------------------------------------------------------------------//
#### Execution: python hub_interval_to_gene_node.py hub_list.txt gene_list.bed output_name.txt
#### Overview: This is a python helper script to the "hub_pipeline.sh" script that takes in a hub (union)list and appends 
# two new columns: 1) annotations for hub status (i.e. lost/gain) and 2) a list of genes (from the provided .bed file) 
# that overlap with the hub genomic interval.
# 
#### Args: .txt list of (union) hubs, .bed list of desired (expressed) gene intervals, and output filename
#
#### Final Output: .txt list that is identical as the input list but has an additional column with gene annotations for each hub
#### -------------------------------------------------------------------------------------------------------------------------//

#### Import necessary python libraries
import string
import sys
import linecache

#### Open input files and count number of lines in both
union_file = open(sys.argv[1], "r") 
gene_file = open(sys.argv[2], "r") 
content = union_file.read()
contentb = gene_file.read()
count_one_list = content.split("\n")
count_two_list = contentb.split("\n")
count_one = 0; count_two = 0

for line in count_one_list:
    count_one+=1

for line in count_two_list:
    count_two+=1

#### Move to beginning of each file and create new file to hold output with header
union_file.seek(0)
gene_file.seek(0)

output_name = sys.argv[3]
output = open(output_name, 'w')

if sys.argv[4] == "single":
  print("Annotating individual condition hub list(s)...")
  output.write("chr" + '\t' + "start" + '\t' + "end" + '\t' + "Hub_id" + '\t' + "Edge_Count" + '\t' + "Node_Count" + '\t' + "Genes" + "\n")
  single_run = True
else:
  print("Annotating union hub list...")
  output.write("chr" + '\t' + "start" + '\t' + "end" + '\t' + "Control_Hub_id" + '\t' + "Control_Edge" + '\t' + "Control_Node" + '\t' + "Case_Hub_id" + '\t' + "Case_Edge" + '\t' + "Case_Node" + '\t' + "log2_FC" + '\t' + "Genes" + '\t' + "Functional_Status" + '\t' + "Union_Status" + "\n")
  single_run = False

#### Define list var to hold intersecting genes for each hub
genes = []

#### Iterate through each hub in merged list, and categorize it as normal/lost/gain/denovo lost/denovo gain based on interaction number delta
for i in range(1, count_one):
    
    # Read specific line, split it into fields, and store relevant data
    aline = linecache.getline(sys.argv[1], i)
    afields = aline.split('\t')
    
    if single_run:
      achrom = afields[0]; astart = int(afields[1]); aend = int(afields[2]); DMSO_Clique = afields[3]; DMSO_Edge = int(afields[4]); DMSO_Node = int(afields[5])
    else:
      achrom = afields[0]; astart = int(afields[1]); aend = int(afields[2]); DMSO_Clique = afields[3]; DMSO_Edge = int(afields[4]); DMSO_Node = int(afields[5]); SHOCK_Clique = afields[6]; SHOCK_Edge = int(afields[7]); SHOCK_Node = int(afields[8]); log2_FC = afields[9]; status = afields[10]
      # Annotate hub as de novo lost or gained (or neither) based on interaction count
      if SHOCK_Edge < 6:
          func_status = 'lost'
          if SHOCK_Edge == 0:
              func_status = 'dnlost'
      elif DMSO_Edge < 6:
          func_status = 'gain'
          if DMSO_Edge == 0:
              func_status = 'dngain' 
      else:
          func_status = ''

    # Compare this hub to gene intervals in file2 to determine if intersection
    for k in range (1, count_two):

        # Read intersecting gene in file2, split it into fields, and store relevant data
        bline = linecache.getline(sys.argv[2], k)
        bfields = bline.split('\t')
        bchrom = bfields[0]; bstart = int(bfields[1]); bend = int(bfields[2]); gene1 = str(bfields[3]); gene = gene1.strip("\n")
        
        # Define booleans for hub and gene interval overlap
        intersect = (achrom == bchrom) and ((bstart >= astart and bstart <= aend) or (bend >= astart and bend <= aend) or ((bstart <= astart) and (bend >= aend)))

        # If intersection, add gene to list of genes that are contained within this hub
        if (intersect):
            genes.append(gene)
                
    # When done iterating through gene file, add all intersecting gene names to outputfile and reset genes list iterator
    gene_list_string = ','.join(genes)
    end = 0
    if single_run:
      output.write(achrom + '\t' + str(astart) + '\t' + str(aend) + '\t' + DMSO_Clique + '\t' + str(DMSO_Edge) + '\t' + str(DMSO_Node) + '\t' + ','.join(genes) + '\n')
      genes = []
    else: 
      slog2_FC = str(log2_FC).strip()
      output.write(achrom + '\t' + str(astart) + '\t' + str(aend) + '\t' + DMSO_Clique + '\t' + str(DMSO_Edge) + '\t' + str(DMSO_Node) + '\t' + SHOCK_Clique + '\t' + str(SHOCK_Edge) + '\t' + str(SHOCK_Node) + '\t' + str(slog2_FC) + '\t' + ','.join(genes) + '\t' + str(func_status) + '\t' + str(status))
      genes = []

#### Close remaining files
union_file.close()
gene_file.close()
output.close()
