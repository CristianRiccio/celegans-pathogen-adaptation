# Export all the 'gene' lines in the gtf file to another file
# gtf = open('output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head50.gtf')
gtf = open(snakemake.input.gtf)
# gtfGenes = open('output/genome/celegans/ref/ensembl/release95/annotation/gtf/Caenorhabditis_elegans.WBcel235Head50gene.gtf', 'w')
gtfGenes = open(snakemake.output.gtfFeature, 'w')

# Create a variable that stores strings of the header in a list
header = []

for line in gtf:
    lineStrip = line.strip()
    # Store the first character of the line in a variable
    firstChar = lineStrip[0]
    # If the first character of the line is a '#', then the line is part of the header and is added 
    # to the header variable
    if firstChar == '#':
        header.append(line)
    if firstChar != '#':
        lineSplit = line.split('\t')
        feature = lineSplit[2]
        print('The feature on this line is: ', feature)
        if feature == 'gene':
            gtfGenes.write(line)
print('The header is ', header)
