import cairo as c
import re
import argparse
#JAMES - FIX ARGPARSE BEFORE TURNING IT IN!!!!! 
def get_arguments():
	ap = argparse.ArgumentParser(description="Motif Marker is used to find motif sequences within a given DNA sequence. DNA bases in lower-case are assumed to be introns and capitalized bases are exons.")
	ap.add_argument("-f", "--file", type=str, help="FASTA input file")
	ap.add_argument("-m", "--mfile", type=str, help="Motif input file - 1 motif per line")
	return ap.parse_args()

args = get_arguments()
ifile = args.file
mfile = args.mfile

def generate_map(gene_objects, motif_objects, motif_count):  #draw data to image
    mapcount = len(gene_objects) + 1
    w = 2200
    h = 250 + mapcount * 200
    #Create canvas
    canvas = c.ImageSurface(c.FORMAT_RGB24, w, h)
    ctx = c.Context(canvas)
    #Set background
    ctx.rectangle(0, 0, w, h)
    ctx.set_source_rgb(0,0,0)
    ctx.fill()
    #Title 
    titlesrc = c.LinearGradient(0, 10, 0, 40)
    titlesrc.set_extend(c.EXTEND_REPEAT)
    titlesrc.add_color_stop_rgb(0.0, 0.8, 0.8, 1)
    titlesrc.add_color_stop_rgb(0.5, 1, 0.8, 1) 
    ctx.select_font_face("Arial", c.FONT_WEIGHT_BOLD)
    ctx.set_font_size(60)
    ctx.move_to((w/2)-150, 70)
    ctx.text_path("MOTIF MARKER")
    ctx.set_source(titlesrc)
    ctx.fill()

    #Motif color map - if more than 6 motifs, add colors here!
    colormap = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]]
                 #Red    Green   Blue    Yellow  Pink   Turquoise 
    #Legend (Supports max of 6 motifs!)
    ctx.rectangle(675, 100, 950, 200)
    ctx.set_source_rgb(1,1,1)
    ctx.set_line_width(2)
    ctx.stroke()
    ctx.set_line_width(10)
    ctx.set_font_size(28)
    for i in range(motif_count):
        ctx.set_source_rgb(colormap[i][0],colormap[i][1],colormap[i][2])
        if i < 3: #column 1 
            legend_Y = 100+(50*(i+1))
            x = 575 #used to center the legend
            ctx.move_to(150+x, legend_Y+10)
            ctx.show_text(motif_objects[i].ID)
            ctx.move_to(300+x, legend_Y)
            ctx.line_to(500+x, legend_Y)
            ctx.stroke()
        else:  #column 2
            legend_Y = -50+(50*(i+1))
            ctx.move_to(560+x, legend_Y+10)
            ctx.show_text(motif_objects[i].ID)
            ctx.move_to(800+x, legend_Y)
            ctx.line_to(1000+x, legend_Y)
            ctx.set_source_rgb(colormap[i][0],colormap[i][1],colormap[i][2])
            ctx.stroke()
    #end legend

    counter = -1 #obligatory counter to make sure motif marks go to correct gene (used for motif.gene)
    for g, v in enumerate(gene_objects):  #for each gene found, map it
        y = 200 * g + 500  #set Y for current map
        seq = v.sequence #grab the sequence from the g'th gene
        seq = list(seq) #split the sequence
        for i,v in enumerate(seq): #converts each character to either (i)ntron or (e)xon.
            if seq[i].isupper():
                seq[i] = "e"
            else:
                seq[i] = "i"
        seq = "".join(seq)
        seq_lst = [] 
        for l in seq:
            seq_lst.append(l)  #append our i's and e's to a list, used for drawing
        prev = seq_lst[0]
        found_lst = []
        count_lst = []
        count = 0

        while len(seq_lst):
            if seq_lst[0] == prev:
                count += 1  #count current letter unless it has changed (from i to e)
                seq_lst.remove(seq_lst[0])  #remove from list when processed
            else:
                found_lst.append(prev)
                count_lst.append(count)
                count = 1
                prev = seq_lst[0]
                seq_lst.remove(seq_lst[0])
        found_lst.append(prev)  #for quiet tracking of found letters (should be [i, e, i])
        count_lst.append(count)
        #Gene Name
        ctx.set_source_rgb(1,1,1)
        ctx.move_to(100, y-60)
        ctx.show_text(gene_objects[g].ID) 
        #First intron
        ctx.move_to(100, y) 
        ctx.line_to((100+(2*count_lst[0])), y)
        #Exon
        ctx.rectangle((100+(2*count_lst[0])), (y-50), (2*count_lst[1]), 100)
        #Second intron
        ctx.move_to((100+(2*count_lst[0])+(2*count_lst[1])), y)
        ctx.line_to((100+(2*count_lst[0])+(2*count_lst[1])+(2*count_lst[2])), y)
        #Stroke
        ctx.set_source_rgb(1,1,1)
        ctx.set_line_width(5)
        ctx.stroke()
        #Draw motif ticks
        counter += 1
        for motif in motif_objects:
            if motif.gene == counter: #draw all motifs with
                for s in motif.start:
                    tick_size = motif.length * 2
                    ctx.rectangle((100+2*s), ((counter*200+500)-50), tick_size, 100)
                    ctx.set_source_rgb(colormap[motif.number][0],colormap[motif.number][1],colormap[motif.number][2])
                    ctx.set_line_width(1)
                    ctx.fill()
    canvas.write_to_png(ifile + ".png")
    return 0

class Gene:
    def __init__(self, ID, chrome, sequence, exon):
        self.ID = ID #name of gene
        self.chrome = chrome
        self.sequence = sequence
        self.exon = exon

class Motif:
    def __init__(self, ID, gene, start, length, number):
        self.ID = ID #motif name
        self.gene = gene #gene number
        self.start = start #list of start positions
        self.length = length #length of motif
        self.number = number #which motif (i.e. line number)

def find_motifs(gene_objects):
    global motif_count
    motif_count = 0
    with open(mfile, 'r') as motifs:
        for line in motifs:
            motif_lst.append(line.strip())
            motif_count += 1
    geneNUM = -1
    for g in gene_objects:
        seq = g.sequence.upper() #temporarily convert sequence to upper so we can find patterns
        geneNUM += 1
        for i, m in enumerate(motif_lst):  #m is current motif, i will eventually be motif.number
            m = m.upper() #because our motifs are random cases, thanks Leslie
            motifID = m
            motifGENE = geneNUM #number of gene, used for mapping
            motifLEN = len(m)  #length of motif, used for line thickness
            motifSTART = []    #all start points for the motif within the gene
            pattern = (m.replace("N", "[ATCG]").replace("U", "[A]")  #nucleic acid notation
                    .replace("Y", "[CT]").replace("R", "[AG]")
                    .replace("W", "[AT]").replace("S", "[CG]")
                    .replace("K", "[GT]").replace("M", "[AC]")
                    .replace("B", "[CGT]").replace("D", "[AGT]")
                    .replace("H", "[ACT]").replace("V", "[ACG]"))
            matches = re.finditer(pattern, seq)  #find matches in sequence
            for match in matches: #for each match, create a new start point
                motifSTART.append(match.start())
            this_motif = Motif(motifID, motifGENE, motifSTART, motifLEN, i) #create object for current motif
            motif_objects.append(this_motif) #add object to our list
            #we should now have a list of X objects, where X is # of motifs * # of genes
            #4 genes * 4 motifs = 16 objects

genes = []
gene_objects = []
motif_objects = []
motif_lst = []
genecount = -1
seq = ""
is_seq = 0      

with open(ifile, 'r') as sequences:
    for line in sequences:
        if line.startswith(">"):
            if is_seq == 1: #if pending sequence
                for i in seq:
                    if i.isupper():
                        exon += i
                this_gene = Gene(genes[genecount], chrome, seq, exon)
                gene_objects.append(this_gene)
                exon = ""
                genecount += 1
                line.strip()
                is_rev = re.search("reverse", line)
                this_gene = "".join(re.findall(r">[A-Za-z0-9]+\s", line)).replace(">","")
                genes.append(this_gene)
                chrome = "".join(re.findall(r"\s[a-z0-9]+:", line)).replace(":","").replace(" ","")
                seq = ""
            else: #if no pending sequence (aka is first sequence in file)
                is_rev = re.search("reverse", line)
                exon = ""
                genecount += 1
                line.strip()
                this_gene = "".join(re.findall(r">[A-Za-z0-9]+\s", line)).replace(">","")
                genes.append(this_gene)
                chrome = "".join(re.findall(r"\s[a-z0-9]+:", line)).replace(":","").replace(" ","")
                seq = ""
                is_seq = 1             
        else: #if seq line
            seq += line.strip()
    if genes: #process last gene from file
        for i in seq:
            if i.isupper():
                exon += i
        this_gene = Gene(genes[genecount], chrome, seq, exon)
        gene_objects.append(this_gene) 

print("Processing sequences...")
print("Finding motifs...")
print("Reticulating splines...")
find_motifs(gene_objects)
print("Bending spoons...")
print("Filtering morale...")
print("Swapping time and space...")
print("Generating sequence maps...")
generate_map(gene_objects, motif_objects, motif_count)
print("Done!")