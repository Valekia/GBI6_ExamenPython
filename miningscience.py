
#Angiely Camacho 

import re
import pandas as pd
import matplotlib.pyplot as plt
import csv as csv
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank
    
def download_pubmed (keyword):
    """
    Este script permite buscar articulos en pubmed por medio de  palabras claves
    """
     
    Entrez.email = 'gualapuro.moises@gmail.com'
    busq = Entrez.read(Entrez.esearch(db="pubmed", 
                            term=keyword,
                            usehistory="y"))
    webenv = busq["WebEnv"]
    query_key = busq["QueryKey"]
    handle = Entrez.efetch(db="pubmed",
                           rettype="medline", 
                           retmode="text", 
                           retstart=0,
                           retmax=543, webenv=webenv, query_key=query_key)
    data = handle.read()
    dataexp = re.sub(r'\n\s{6}','', data)
    return dataexp
    


def science_plots(tipo):
    """creación de un pie con los países"""
    #if tipo == "AD":
    with open(tipo) as f:
        my_text = f.read()
    my_text = re.sub(r'\n\s{6}', ' ', my_text)  
    zipcodes = re.findall(r'[A-Z]{2}\s(\d{5}), USA', my_text)
    unique_zipcodes = list(set(zipcodes))
    zip_coordinates = {}
    with open('zip_coordinates.txt') as f:
        csvr = csv.DictReader(f)
        for row in csvr:
         zip_coordinates[row['ZIP']] = [float(row['LAT']), float(row['LNG'])]
    zip_code = []
    zip_long = []
    zip_lat = []
    zip_count = []
    for z in unique_zipcodes:
    # if we can find the coordinates
        if z in zip_coordinates.keys():
            zip_code.append(z)
            zip_lat.append(zip_coordinates[z][0])
            zip_long.append(zip_coordinates[z][1])
            zip_count.append(zipcodes.count(z))
    import matplotlib.pyplot as plt
    #%matplotlib inline
    fig = plot.figure()
    plt.scatter(zip_long, zip_lat, s = zip_count, c= zip_count)
    plt.colorbar()
# only continental us without Alaska
    plt.xlim(-125,-65)
    plt.ylim(23, 50)
# add a few cities for reference (optional)
    ard = dict(arrowstyle="->")
    plt.annotate('Houston', xy = (-95.36327, 29.76328), 
                   xytext = (-95.36327, 29.76328), arrowprops = ard)
    plt.annotate('San Diego', xy = (-117.16472, 32.71571), 
                   xytext = (-117.16472, 32.71571), arrowprops= ard)
    plt.annotate('San Francisco', xy = (-122.41942, 37.77493), 
                   xytext = (-122.41942, 37.77493), arrowprops= ard)
    plt.annotate('Manhattan', xy = (-73.96625, 40.78343), 
                   xytext = (-73.96625, 40.78343), arrowprops= ard)
    plt.annotate('Jacksonville', xy = (-81.65565, 30.33218), 
                   xytext = (-81.65565, 30.33218), arrowprops= ard)
    plt.annotate('Miami', xy = (-80.21, 25.7753), 
                   xytext = (-80.21, 25.7753), arrowprops= ard)
    params = plt.gcf()
    plSize = params.get_size_inches()
    params.set_size_inches( (plSize[0] * 3, plSize[1] * 3) )
    return plt.show()