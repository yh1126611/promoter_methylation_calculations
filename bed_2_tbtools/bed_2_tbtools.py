# Converts .bed file to BINStat.tab.xls format (input for TBtools circos track)
# Written by Young Ho Lee, yh1126@snu.ac.kr

# Custom parameters
bed_filename = input("Enter name of input bed file (Extension: .bed):")
tbtools_filename = input("Enter name of output file (Recommended extension: .BINStat.tab.xls):")

bed = open(bed_filename, "r") # If you are using the script on local, enter name of input .bed file here directly.
tbtools = open(tbtools_filename, "w") # If you are using the script on local, enter name of output .BINStat.tab.xls file here directly.
window = 10000 # Enter window size (normally set to 10,000 as default by TBtools)

current_chromosome = ""
start = 0
end = 0
mover = 1
counter = 0

for line in bed:
    start = int(line.split("\t")[1])
    end = int(line.split("\t")[2])

    if current_chromosome != line.split("\t")[0]: # ★ Case I: mover is in another chromosome
        if current_chromosome == "":
            mover = 1
            counter = 0
            current_chromosome = line.split("\t")[0]
        else:
            tbtools.write(str("\t".join([str(current_chromosome), str(mover - (mover % window)), str(mover), str(counter), "\n"])))
            mover = 1
            counter = 0
            current_chromosome = line.split("\t")[0]
    if current_chromosome == line.split("\t")[0]:
        while mover < start: # ★ Case II: mover is before start
            if mover%window == 0:
                tbtools.write(str("\t".join([str(current_chromosome), str(mover - window), str(mover), str(counter), "\n"])))
                counter = 0
            mover = mover + 1
        while mover >= start and mover < end: # ★ Case III: mover is in b/w start & end
            counter = counter + 1
            if mover%window == 0:
                tbtools.write(str("\t".join([str(current_chromosome), str(mover - window), str(mover), str(counter), "\n"])))
                counter = 0
            mover = mover + 1
tbtools.write(str("\t".join([str(current_chromosome), str(mover - (mover % window)), str(mover), str(counter), "\n"])))

bed.close()
tbtools.close()
