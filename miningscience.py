
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
    
#Angiely Camacho

def science_plots (file):
    """
    Codigo para seleccionar los paises de mayor a menor incidencia  
    """
    
    email = re.sub(r'\s[\w._%+-]+@[\w.-]+\.[a-zA-Z]{1,4}','',file)
    puntos = re.sub(r'\..\d.\,',',',email)
    numb = re.sub(r'\..\d.','',puntos)
    x=numb[1:].split('PMID-')
    
    Countries_A=[]
    for PMID in x:
        q=PMID.split('\n')
        for fila in q:
            w=fila.split(' ')
            if w[0] == 'AD':
                e=fila.split(',')
                Countries_A.append(e[-1])
    
    a=0
    Countries_B =[0]*len(Countries_A)
    for lis in Countries_A:
        bytes(lis,encoding="utf8")
        if lis != '':
            w=lis
            if w[0] == ' ':
                w = re.sub (r'^\s','',w)
            if w[-1] == '.':
                w = re.sub (r'\.$','',w)
            w = re.sub (r'\.$','',w)
            w = re.sub (r'\s$','',w)
        Countries_B[a]=w
        a=a+1
        
    Contries_all=['Andorra','United Arab Emirates ','Afghanistan','Antigua and Barbuda','Anguilla','Albania','Armenia','Netherlands Antilles','Angola','Antarctica','Argentina','American Samoa','Austria','Australia','Aruba','Azerbaijan','Bosnia and Herzegovina','Barbados','Bangladesh','Belgium','Burkina Faso','Bulgaria','Bahrain', 'Burundi','Benin','Bermuda','Brunei','Bolivia', 'Brazil','Bahamas','Bhutan','Bouvet Island','Botswana','Belarus','Belize','Canada','Cocos [Keeling] Islands','Congo [DRC]','Central African Republic','Congo Republic', 'Switzerland',"Côte d'Ivoire",'Cook Islands','Chile','Cameroon','China','Colombia','Costa Rica','Cuba', 'Cape Verde','Christmas Island','Cyprus','Czech Republic','Germany','Djibouti','Denmark','Dominica','Dominican Republic','Algeria','Ecuador' ,'Estonia','Egypt','Western Sahara','Eritrea','Spain','Ethiopia','Finland','Fiji','Falkland Islands [Islas Malvinas]','Micronesia','Faroe Islands','France','Gabon', 'United Kingdom','Grenada','Georgia','French Guiana','Guernsey','Ghana','Gibraltar','Greenland','Gambia', 'Guinea','Guadeloupe','Equatorial Guinea','Greece','South Georgia and the South Sandwich Islands','Guatemala','Guam','Guinea-Bissau','Guyana','Gaza Strip','Hong Kong','Heard Island and McDonald Islands','Honduras','Croatia', 'Haiti','Hungary','Indonesia','Ireland' ,'Israel','Isle of Man','India','British Indian Ocean Territory','Iraq', 'Iran','Iceland','Italy','Jersey','Jamaica','Jordan', 'Japan','Kenya','Kyrgyzstan','Cambodia','Kiribati','Comoros','Saint Kitts and Nevis','North Korea','South Korea','Kuwait','Cayman Islands','Kazakhstan','Laos','Lebanon','Saint Lucia','Liechtenstein','Sri Lanka','Liberia','Lesotho','Lithuania','Luxembourg','Latvia' ,'Libya','Morocco','Monaco','Moldova','Montenegro','Madagascar','Marshall Islands','Macedonia [FYROM]','Mali','Myanmar [Burma]','Mongolia' ,'Macau','Northern Mariana Islands','Martinique','Mauritania','Montserrat','Malta','Mauritius','Maldives','Malawi','Mexico','Malaysia' ,'Mozambique','Namibia','New Caledonia','Niger','Norfolk Island','Nigeria','Nicaragua','The Netherlands','Norway','Nepal','Nauru', 'Niue','New Zealand','Oman','Panama','Peru','French Polynesia', 'Papua New Guinea','Philippines','Pakistan','Poland','Saint Pierre and Miquelon' ,'Pitcairn Islands','Puerto Rico','Palestinian Territories','Portugal','Palau','Paraguay','Qatar','Réunion','Romania', 'Serbia','Russia' ,'Rwanda','Saudi Arabia','Solomon Islands','Seychelles','Sudan','Sweden','Singapore','Saint Helena','Slovenia', 'Svalbard and Jan Mayen','Slovakia','Sierra Leone','San Marino','Senegal','Somalia','Suriname','São Tomé and Príncipe','El Salvador','Syria', 'Swaziland' ,'Turks and Caicos Islands','Chad','French Southern Territories','Togo','Thailand','Tajikistan','Tokelau','Timor-Leste','Turkmenistan' ,'Tunisia','Tonga','Turkey','Trinidad and Tobago','Tuvalu','Taiwan','Tanzania','Ukraine','Uganda','U.S. Minor Outlying Islands','United States of America','Uruguay','Uzbekistan','Vatican City','Saint Vincent and the Grenadines','Venezuela', 'British Virgin Islands','U.S. Virgin Islands','Vietnam','Vanuatu','Wallis and Futuna','Samoa','Kosovo','Yemen','Mayotte','South Africa','Zambia','Zimbabwe']
    Countries_C=Countries_B
    h=Contries_all
    f=len(h)
    CountriesCount=[0]*f
    k=0
    for elem in h:
        d=0
        for comp in Countries_C:
            if elem == str(comp):
                d=d+1
        CountriesCount[k]=d
        k=k+1
    
    Countries_D=[]
    Counter=[]
    o=0
    for line in CountriesCount:
        if str(line) != '0': 
            Counter.append(line) 
            m=Contries_all[o]
            Countries_D.append(m) 
        o=o+1

    Table_A = pd.DataFrame({'Country' : Countries_D,
                           'num_auth' : Counter})
    Order=Table_A.sort_values(by=['num_auth'], ascending=[False])
    Taken = Order.iloc[0:5]
    suma=Taken['num_auth'].sum()
    iu= Taken.iloc[:,0]
    su=pd.Series(iu)
    i = Taken.iloc[:,1]
    sa = pd.Series(i)
    s = sa.tolist()
    
    prom = []
    for number in s:
        xa=(number/suma)*100
        prom.append(xa)
    Table_B = pd.DataFrame({'Country' : su,
                           'Porcent' : prom})
    
    fig1, ax1 = plt.subplots()
    ax1.pie(prom, labels=su, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  
    plt.show()

    return (Table_B)