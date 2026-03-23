
# extract ROI from TetR/TetA plasmids matches

# read in sequence

with open('SPARK_1464_C1_2.fa', 'r') as file:
    # save the first line as 'id'
    id = next(file)
    # read the rest of the file in as the sequence
    seq = file.read()
    
# extract ROI - round the values to outer 100 bp
lower_bound = 122100

upper_bound = 128200

roi = seq[lower_bound:upper_bound]

# need to remove newlines already in it 

# remove newline characters from roi for better formatting 
roi = roi.replace('\n', '')

# save roi of each plasmid for use in clinker alignment and aliview tree generation

with open('roi_test/SPARK_1464_C1_2_roi.fa', 'w') as file:
    # header
    file.write(id)
    # sequence
    for i in range(0, len(roi), 70):
        chunk = roi[i:i + 70]
        line = chunk + '\n'
        file.write(line)
    file.close()

