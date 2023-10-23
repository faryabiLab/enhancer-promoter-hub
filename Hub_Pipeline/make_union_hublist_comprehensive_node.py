#### --------------------------------------------------------------------------------------------------------------------------//
#### Execution: python make_union_hublist_comprehensive_node.py Control_hub_by_interval_sort.txt Case_hub_by_interval_sort.txt output.txt
#### Overview: This is a python helper script to the "hub_pipeline.sh" script that generates the unionlist of hubs 
# between the control and case conditions by mapping hubs between the conditions onto one another on the basis of overlapping
# genomic intervals. The interaction count and node count of the union hubs in each condition is then recalculated as necessary. 
# 
#### Dependency: requires bedtools 
# 
#### Args: 2 sorted .txt datatables of hubs in control and case conditions (respectively) annotated with interaction and node counts, 
# and sorted first by ascending chromosome and then by ascending start coordinate with no header, and output filename
#
#### Final Output: .txt datatable containing all union hubs as well as their attributes in each condition (i.e. interaction count, node count, 
# merge status, etc.)
#### -------------------------------------------------------------------------------------------------------------------------//

#### Define clique (i.e. hub) objects
class clique:
	def __init__(self, id, chr, start, end, edge, node, condition):
		self.id = id
		self.chr = chr
		self.start = start
		self.end = end
		self.edge = edge
		self.node = node
		self.condition = condition

	def __str__(self):
		return f"{self.id}\t{self.chr}\t{self.start}\t{self.end}\t{self.edge}\t{self.node}"

class union_clique:
	def __init__(self, id_1, id_2, chr, start, end, edge_1, edge_2, node_1, node_2, FC, status):
		self.id_1 = id_1
		self.id_2 = id_2
		self.chr = chr
		self.start = start
		self.end = end
		self.edge_1 = edge_1
		self.edge_2 = edge_2
		self.node_1 = node_1
		self.node_2 = node_2
		self.FC = FC	
		self.status = status
	def __str__(self):
		return f"{self.chr}\t{self.start}\t{self.end}\t{self.id_1}\t{self.edge_1}\t{self.node_1}\t{self.id_2}\t{self.edge_2}\t{self.node_2}\t{self.FC}\t{self.status}\n"

import math
import sys
import linecache
import os; import subprocess

#### Open input files and count numbers of lines/cliques
file_one = open(sys.argv[1], "r") 
file_two = open(sys.argv[2], "r") 
bedfile = open("int_file.bed", "w")
unionfile = open(sys.argv[3], 'w')
content = file_one.read()
contentb = file_two.read()
count_one_list = content.split("\n")
count_two_list = contentb.split("\n")
count_one = 0; count_two = 0

for line in count_one_list:
	count_one+=1

for line in count_two_list:
	count_two+=1
print("Line number is: ", count_one, count_two)

#### Move to beginning of input files
file_one.seek(0)
file_two.seek(0)

#### Define lists that will hold clique objects from each condition
cliqlist = []
unioncliqlist = []

#### Iterate through each clique in list, create clique objects for condition 1, & write to .bed file
for i in range(1, count_one):
	
	# Read specific line, split it into fields, and store relevant data
	aline = linecache.getline(sys.argv[1], i)
	afields = aline.split('\t')
	achrom = afields[0]; astart = int(afields[1]); aend = int(afields[2]); aname = afields[3]; aedge = int(afields[4]); anode = int(afields[5])
	
	if int(aedge) == 0:
		continue

	bedfile.write(achrom + '\t' + str(astart) + '\t' + str(aend) + "\n")
	cliqlist.append(clique(aname, achrom, astart, aend, aedge, anode, "control"))

####  Iterate through each clique in list, create clique objects for condition 2, & write to .bed file
for i in range(1, count_two):
	
	# Read specific line, split it into fields, and store relevant data
	aline = linecache.getline(sys.argv[2], i)
	afields = aline.split('\t')
	achrom = afields[0]; astart = int(afields[1]); aend = int(afields[2]); aname = afields[3]; aedge = int(afields[4]); anode = int(afields[5])

	if int(aedge) == 0:
		continue

	bedfile.write(achrom + '\t' + str(astart) + '\t' + str(aend) + "\n")
	cliqlist.append(clique(aname, achrom, astart, aend, aedge, anode, "case"))

#### Run bedtools merge on sorted .bed file containing intervals of all cliques (including redundancies)
# ASSUMES BEDTOOLS IS ALREADY LOADED INTO ENVIRONMENT
bedfile.close()
cmd_1 = 'sort -k1,1 -k2,2n int_file.bed > int_file_sorted.bed'
cmd_2 = 'bedtools merge -i int_file_sorted.bed > bedtools_output.bed'
os.system(cmd_1); os.system(cmd_2)

#### Open merged interval .bed file
merge_file = open("bedtools_output.bed", "r")

#### Define variables to count number of 'special case' cliques between conditions
edge_num1 = 0; edge_num2 = 0
lost_num = 0; gain_num = 0
status = ""

#### Iterate through merged file and create union clique objects/list
for line in merge_file:
	# Divide interval into chr, start, and end for comparison
	interval = line.split("\t"); mchr = interval[0]; mstart = int(interval[1]); mend = int(interval[2])

	# Define vars to take on final values and iterators
	intersect_id1 = ""; intersect_id2 = ""
	intersect_edge1 = 0; intersect_edge2 = 0
	intersect_node1 = 0; intersect_node2 = 0  
	a = 0; b = 0 # iterators to inform clique type (i.e. gain // lost)

	# Iterate through cliquelist and create unionclique objects based on overlap
	for c1 in cliqlist:
		intersect = (c1.chr == mchr) and ((mstart >= c1.start and mstart <= c1.end) or (mend >= c1.start and mend <= c1.end) or ((mstart < c1.start) and (mend > c1.end)))	
		if (intersect):
			if (c1.condition == "control"):
				# Create single string to hold 2+ clique ids
				intersect_id1 = f"{c1.id}_{intersect_id1}"

				# Find new intersecting edge number and log2FC
				intersect_edge1 += int(c1.edge)
				intersect_node1 += int(c1.node)
				a += 1

			elif (c1.condition == "case"):
				# Create single string to hold 2+ clique ids
				intersect_id2 = f"{c1.id}_{intersect_id2}"

				# Find new intersecting edge number and log2FC
				intersect_edge2 += int(c1.edge)
				intersect_node2 += int(c1.node)
				b += 1
	
	# Count/increment numbers of each special cliques
	if (a == 1) and (b == 1):
		status = "normal"

	else: 
		if a == 0:
			gain_num += 1
			status = "gained"
			log2_FC = " "
			intersect_id1 = "N/A"
		elif (a > 1) and (b == 1):
			edge_num1 += 1
			status = "merge_in_case"
		elif b == 0:
			lost_num += 1
			status = "lost"
			log2_FC = " "
			intersect_id2 = "N/A"
		elif (b > 1) and (a == 1):
			edge_num2 += 1
			status = "fragment_in_case"
		elif (b > 1) and (a > 1):
			status = "complex"
		elif (a == 0) and (b == 0):
			print("ERROR: INTERVALS WITH NO CLIQUE OVERLAP REPORTED")
			break
	
	# Calculate log2FC in edge number and create union_clique object
	if not (status == "lost" or status == "gained"):
		log2_FC = math.log2(int(intersect_edge2)) - math.log2(int(intersect_edge1)) 
	
	intersect_id1r = intersect_id1.strip('_'); intersect_id2r = intersect_id2.strip('_')
	
	u1 = union_clique(intersect_id1r,intersect_id2r,mchr,mstart,mend,intersect_edge1,intersect_edge2,intersect_node1,intersect_node2,log2_FC,status)
	unioncliqlist.append(u1)


####  Write unioncliqlist to file 
for c2 in unioncliqlist:

	# Eliminate cliques from list that have less than 6 edges in both condition 1 AND condition 2
	if (c2.edge_1 < 6) and (c2.edge_2 < 6):
		continue
	else: 
		unionfile.write(str(c2))

# Close remaining files
file_one.close()
file_two.close()
merge_file.close()
unionfile.close()

