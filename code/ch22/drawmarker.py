import os
import re

from pymongo import MongoClient
from Bio.Graphics import BasicChromosome
from reportlab.lib import colors
from reportlab.lib.units import cm

CONNECTION_STRING = os.getenv('MONGODB_CS', 'localhost:27017')

def sortmarkers(crms,end):
    """ Sort markers into chromosomes """
    i = 0
    crms_o = [[] for r in range(len(end))]
    crms_fo = [[] for r in range(len(end))]
    for crm in crms:
        for marker in crm:
            # add the marker start position at each chromosome.
            crms_fo[i].append(marker[1])
        crms_fo[i].sort() # Sort the marker positions.
        i += 1
    i = 0
    for order in crms_fo:
        # Using the marker order set in crms_fo, fill crms_o
        # with all the marker information
        for pos in order:
            for mark in crms[i]:
                if pos==mark[1]:
                    crms_o[i].append(mark)
        i += 1
    return crms_o

def getchromo(crms_o, end):
    """ From an ordered list of markers, generate chromosomes.
    """
    chromo = [[] for r in range(len(end))]
    i = 0
    for crm_o in crms_o:
        j = 0
        if len(crm_o)>1:
            for mark in crm_o:
                if mark==crm_o[0]: #first marker
                    chromo[i].append(('', None, mark[1]))
                    chromo[i].append((mark[0], colors.red,
                                      mark[2]-mark[1]))
                    ant = mark[2]
                elif mark==crm_o[-1]: #last marker
                    chromo[i].append(('', None, mark[1]-ant))
                    chromo[i].append((mark[0], colors.red,
                                      mark[2]-mark[1]))
                    chromo[i].append(('', None, end[i]-mark[2]))
                else:
                    chromo[i].append(('', None, mark[1]-ant))
                    chromo[i].append((mark[0], colors.red,
                                      mark[2]-mark[1]))
                    ant=mark[2]
        elif len(crm_o)==1: # For chromosomes with one marker
            chromo[i].append(('', None, crm_o[0][1]))
            chromo[i].append((crm_o[0][0], colors.red,
                              crm_o[0][2]-crm_o[0][1]))
            chromo[i].append(('', None, end[i]-crm_o[0][2]))
        else:
            # For chromosomes without markers
            # Add 3% of each chromosome.
            chromo[i].append(('', None, int(0.03*end[i])))
            chromo[i].append(('', None, end[i]))
            chromo[i].append(('', None, int(0.03*end[i])))
        i += 1
        j += 1
    return chromo

def addends(chromo):
    """ Adds a 3% of blank region at both ends for better
        graphic output.
    """
    size = 0
    for x in chromo:
        size += x[2]
    # get 3% of size of each chromosome:
    endsize = int(float(size)*.03)
    # add this size to both ends in chromo:
    chromo.insert(0,('', None, endsize))
    chromo.append(('', None, endsize))
    return chromo

def load_chrom(chr_name):
    """ Generate a chromosome with information
    """
    cur_chromosome = BasicChromosome.Chromosome(chr_name[0])
    chr_segment_info = chr_name[1]

    for seg_info_num in range(len(chr_segment_info)):
        label, color, scale = chr_segment_info[seg_info_num]
        # make the top and bottom telomeres
        if seg_info_num == 0:
            cur_segment = BasicChromosome.TelomereSegment()
        elif seg_info_num == len(chr_segment_info) - 1:
            cur_segment = BasicChromosome.TelomereSegment(1)
        # otherwise, they are just regular segments
        else:
            cur_segment = BasicChromosome.ChromosomeSegment()
        cur_segment.label = label
        cur_segment.label_size = 12
        cur_segment.fill_color = color
        cur_segment.scale = scale
        cur_chromosome.add(cur_segment)

    cur_chromosome.scale_num = max(END) + (max(END)*.04)
    return cur_chromosome

def dblookup(atgids):
    """ Code to retrieve all marker data from name using MongoDB
    """
    client = MongoClient(CONNECTION_STRING)
    db = client.pr4
    collection = db.markers_map4
    markers = []
    for marker in atgids:
        mrk = collection.find_one({'marker_id': marker})
        if mrk:
            markers.append((marker, (mrk['chromosome'],
                            mrk['start'], mrk['end'])))
        else:
            print('Marker {0} is not in the DB'.format(marker))
    return markers

# Size of each chromosome:
END = (30427563, 19696817, 23467989, 18581571, 26986107)
gids = []
rx_rid = re.compile('^AT[1-5]G\d{5}$')
print('''Enter AT ID or press 'enter' to stop entering IDs.
Valid IDs:
AT2G28000
AT3G03020

Also you can enter DBDEMO to use predefined set of markers
fetched from a MongoDB database. Enter NODBDEMO to use a
predefined set of markers without database access.''')
while True:
    rid = input('Enter Gene ID: ')
    if not rid:
        break
    if rid=='DBDEMO':
        gids = ['AT3G14890','AT1G66160','AT3G55260','AT5G59570',
                'AT4G32551','AT1G01430','AT4G26000','AT2G28000',
                'AT5G21090','AT5G10470']
        break
    elif rid=='NODBDEMO':
        samplemarkers=[('AT3G14890', (3, 5008749, 5013275)),
                       ('AT1G66160', (1, 24640827, 24642411)),
                       ('AT3G55260', (3, 20500225, 20504056)),
                       ('AT1G10960', (1, 3664385, 3665040)),
                       ('AT5G23350', (5, 7857646, 7859280)),
                       ('AT5G15250', (5, 4950414, 4952780)),
                       ('AT1G55700', (1, 20825263, 20827306)),
                       ('AT5G21090', (5, 7164583, 7167257)),
                       ('AT5G10470', (5, 3289228, 3297249)),
                       ('AT2G28000', (2, 11933524, 11936523)),
                       ('AT3G03020', (3, 680920, 682009)),
                       ('AT4G26000', (4, 13197255, 13199845)),
                       ('AT4G32551', (4, 15707516, 15713587))]
        break
    if rx_rid.match(rid):
        gids.append(rid)
    else:
        print("Bad format, please enter it again")

if rid!='NODBDEMO':
    samplemarkers = dblookup(gids)

crms = [[] for r in range(len(END))]
for x in samplemarkers:
    crms[int(x[1][0])-1].append((x[0], x[1][1], x[1][2]))

crms_o = sortmarkers(crms, END)
chromo = getchromo(crms_o, END)
all_chr_info = [('Chr I', chromo[0]), ('Chr II', chromo[1]),
                ('Chr III', chromo[2]), ('Chr IV', chromo[3]),
                ('Chr V', chromo[4])]

organism = BasicChromosome.Organism()
organism.page_size = (29.7*cm, 21*cm) #A4 landscape
for chr_info in all_chr_info:
    newcrom = (chr_info[0], addends(chr_info[1]))
    organism.add(load_chrom(newcrom))

organism.draw('at.pdf','Arabidopsis thaliana')
