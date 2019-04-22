import os
import datetime
import pandas as pd
from collections import OrderedDict
import requests
import json
import sys


#convert vcf file to string (meta-information) and dataframe (variant information)
def VCFtoDF(fileName):
    processInfo = ''
    df = pd.DataFrame()
    #read vcf line by line. meta-information started with '##' to string. Variant information to df 
    try:
        vcf = open(fileName,'r')
        vcfLines = []
        for line in vcf:
            vcfLines.append(line)
        vcf.close()

        #Get date to record the time for file annotation
        date = datetime.datetime.now().strftime('%Y-%m-%d')
        metaInfo = '##fileAnnotationDate='+date+'\n'

        #Add some other self-defined info to metaInfo
        defineInfo ='''##INFO=<ID=compTwoSampl,Type=Integer,Description="Compare two samples. The value is 1, if all numbers are same; otherwise 0">
##INFO=<ID=Gene Symbol,Type=String,Description="Gene Symbol of gene where variants are located">
##INFO=<ID=CCDS ID,Type=String,Description="CCDS id of canonical transcript of the gene">
##INFO=<ID=Motif Name,Type=String,Description="Motif name if variant is located in a motif">
##INFO=<ID=Mutation type,Type=String,Description="Mutation type caused by variant">
##INFO=<ID=Amino acid change,Type=String,Description="Amino acid change caused by variant in canonical transcript if there is any">
##INFO=<ID=Clincal Significance,Type=String,Description="Classification of variant from ClinVar">
##INFO=<ID=PubMed,Type=String,Description="PubMed ID about variant">
##INFO=<ID=ENSP ID,Type=String,Description="ENSP id of gene">
##INFO=<ID=CDC position,Type=String,Description="variant position on canonical transcript of the gene">
##INFO=<ID=Intron,Type=String,Description="In which intron the variant located, defined by canonical transcript">
##INFO=<ID=Exon,Type=String,Description="In which exon the variant located and total exon numbers, defined by canonical transcript">
##INFO=<ID=STRAND,Type=String,Description="Corresponding DNA strand for gene">
##INFO=<ID=SIFT,Type=String,Description="SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids.">
##INFO=<ID=PolyPhen,Type=String,Description="Predicting possible impact of an amino acid substitution on the structure and function of a human protein">
##INFO=<ID=GTof1st,Type=String,Description="GT of first sample">
##INFO=<ID=DPof1st,Type=Integer,Description="DP of first sample">
##INFO=<ID=ReadSupport1st,Type=Integer,Description="Reads number support annotated variant in first sample">
##INFO=<ID=ratioAO1st,Type=Floating point number,Description="Percentage of alternate allele observation count for annotated variant in first sample">
##INFO=<ID=ratioRO1st,Type=Floating point number,Description="Percentage of reference allele observation count for annotated variant in first sample">
##INFO=<ID=GTof2nd,Type=String,Description="GT of second sample">
##INFO=<ID=DPof2nd,Type=Integer,Description="DP of second sample">
##INFO=<ID=ReadSupport2nd,Type=Integer,Description="Reads number support annotated variant in second sample">
##INFO=<ID=ratioAO2nd,Type=Floating point number,Description="Percentage of alternate allele observation count for annotated variant in second sample">
##INFO=<ID=ratioRO2nd,Type=Floating point number,Description="Percentage of reference allele observation count for annotated variant in second sample">
##INFO=<ID=vepAnnotationMostImpactTranscripts,Type=String,Description="Prediction of most severe impact on any transcript from gene">
##INFO=<ID=vepAnnotationCanonicalTranscript,Type=String,Description="Prediction of impact on canonical transcript from gene">'''
            
        metaInfo = metaInfo + defineInfo
        
        n=0
        for line in vcfLines:
            if(str(line)[0]=='#' and str(line)[1]=='#'):
                metaInfo = metaInfo + str(line)
                n += 1
            else:
                break

        #header for df
        headers = vcfLines[n][:-1].split('\t')

        #data for each column
        dataInfo = [[] for i in range(len(headers))]

        for j in range (n+1,len(vcfLines)):
            lineInfo = vcfLines[j][:-1].split('\t')
            for i in range (0, len(headers)):
                dataInfo[i].append(lineInfo[i])

        #Create dictionary and then to df
        vcfData = OrderedDict()
        for i in range (0, len(headers)):
            vcfData[headers[i]] = dataInfo[i]

        df = pd.DataFrame(vcfData)
        processInfo = 'Success! File convert to dataframe.'

        return([processInfo, metaInfo, df])

    except Exception as e:
        processInfo = 'Error! ' + str(e)
        metaInfo = ''
        return ([processInfo, metaInfo, df])
        



#Query ExAC using variant position, return following items:
#allele_freq, vep_annotation(first one, since it is the most severe impact amongst the transcripts according to guide),
# vep_annotation(canonical, which are we use more often), rsid

def ExacQuery(urlString):
    #Use ExAC API for annotation
    #return as dictionary including allele_freq, rsid, vep_annotation_first, vep_annotation-canonical
    #other relevant information (I think) include: gene(symbol, ccds-id, ensp-id),consequence, protein level change if it is
    # any, clinial significance, pubmed, variant location (exon or intron), function prediction(SIFT, PolyPhen). These information
    # are extracted from canonical vep.
    exacInfo = {}
    exacResponse = requests.get(urlString)
    exacjson = exacResponse.json()
    #if there is allele_freq, return this number; if not, return '.'
    if ("allele_freq" in exacjson):
        allele_frequency = exacjson["allele_freq"]
    else:
        allele_frequency = '.'
    #if there is rsid, return it; sometimes, it can also be found in vep_annotation; else return '.'
    if ("rsid" in exacjson):
        rsID = exacjson["rsid"]
        if(rsID == '.' and ("vep_annotations" in exacjson)):
            rsID = exacjson["vep_annotations"][0]["Existing_variation"]
    else:
        rsID = '.'

    # if there is vep-annotatio, return first one, most severe impact on function prediction 
    GENESYM = '.'
    CCDS = '.'
    MOTIF_NAME = '.'
    Conquecence = '.'
    PROTEINCHANGE = '.'
    CLIN_SIG = '.'
    PUBMED = '.'
    ENSP = '.'
    CDC_postion = '.'
    INTRON = '.'
    EXON = '.'
    STRAND = '.'
    SIFT = '.'
    PolyPhen = '.'

    vepAnnotationFirst = OrderedDict()
    vepAnnotationCanonical = OrderedDict()
    

    if ("vep_annotations" in exacjson):
        try:
            vepAnnotationFirst = exacjson["vep_annotations"][0]
            vepAnnotationCanonical = {}
            for vep in exacjson["vep_annotations"]:
                if (vep["CANONICAL"] == "YES"):
                    vepAnnotationCanonical = vep
                    GENESYM = vep["SYMBOL"]
                    CCDS = vep["CCDS"]
                    MOTIF_NAME = vep["MOTIF_NAME"]
                    Conquecence = vep["Consequence"]
                    PROTEINCHANGE = vep["HGVSp"].split(':')[-1]
                    CLIN_SIG = vep["CLIN_SIG"]
                    PUBMED = vep["PUBMED"]
                    ENSP = vep["ENSP"]
                    CDC_postion = vep["CDS_position"]
                    INTRON = vep["INTRON"] + '.'
                    EXON = vep["EXON"] + '.'
                    STRAND = vep["STRAND"]
                    SIFT = vep["SIFT"]
                    PolyPhen = vep["PolyPhen"]
        except Exception as e:
            print()
            print(str(e))
            print(urlString)
                
    else:
        vepAnnotationFirst = {}
        vepAnnotationCanonical = {}

    
    return {'allele_frequency':allele_frequency,
            'rsID': rsID,
            'vepAnnotationFirst':vepAnnotationFirst,
            'vepAnnotationCanonical':vepAnnotationCanonical,
            'GENESYM':GENESYM,
            'CCDS':CCDS,
            'MOTIF_NAME' : MOTIF_NAME, 
            'Conquecence' : Conquecence,
            'PROTEINCHANGE' : PROTEINCHANGE,
            'CLIN_SIG' : CLIN_SIG,
            'PUBMED' : PUBMED,
            'ENSP' : ENSP,
            'CDC_postion' : CDC_postion,
            'INTRON' : INTRON,
            'EXON' : EXON,
            'STRAND' : STRAND,
            'SIFT' : SIFT,
            'PolyPhen' : PolyPhen}


#Iterate VEP from ExAC (format change)
def IteVEP(vep):
    strFormat = ''
    for key, value in vep.items():
        if(value == ''):
            value = '.'
        strFormat = strFormat + key + ' : ' + value + '||'
    return strFormat


#Add all annotation to vcf dataframe

def Annotation(vcfDF):
    #Information from ExAc query
    alleleFrequency = []
    rsID = []
    vepAnnotationFirst = []
    vepAnnotationCanonical = []
    GENESYM = []
    CCDS = []
    MOTIF_NAME = []
    Conquecence = []
    PROTEINCHANGE = []
    CLIN_SIG = []
    PUBMED = []
    ENSP = []
    CDC_postion = []
    INTRON = []
    EXON = []
    STRAND = []
    SIFT = []
    PolyPhen = []

    #compare result data of firstSample vs. secondSample. If they are identical, set value to 1, else set value to 0
    compTwoSampl=[]

    #variant type
    varinatType = []    

    #Genotype of normal tissue
    GTof1st = []
    #Depth of normal tissue
    DPof1st = []
    #Support variant reads of normal tissue
    ReadSupport1st = []
    #ratio of alternative observation of normal tissue
    ratioAO1st=[]
    #ratio of reference observation of normal tissue
    ratioRO1st=[]

    #Genotype of vaf tissue
    GTof2nd = []
    #Depth of vaf tissue
    DPof2nd = []
    #Support variant reads of vaf tissue
    ReadSupport2nd = []
    #ratio of alternative observation of vaf tissue
    ratioAO2nd=[]
    #ratio of reference observation of vaf tissue
    ratioRO2nd=[]

    n=0

    #iterating every row
    for index, row in vcfDF.iterrows():

        #find most deleterious variant if there are  multiple possibilities
        #usually deletion and insertion cause more deleterious impact.
        #assume the more deletion or insertion will be more deleterious
        #find largest length difference between all multiple possibilities and reference
        #if same, choose first one with largest length diffenrence
    
        #if there is only one type of variant, set index to 0

        variant = row['ALT']
        index = 0
        if (variant.find(',') == -1):
            index = 0

        #find largest length difference between all multiple possibilities and reference
        #set index of variant with largest length difference (it will be use for query, find supporting reads)
        else:
            #calculate length differences
            diff = []
            for v in variant.split(','):
                diff.append(abs(len(v)-len(row['REF'])))
            max_value = max(diff)
            index = diff.index(max_value)
                


        #compare result data of normal vs. vaf5. Any difference will result in value to 0
        if(row['normal'] == row['vaf5']):
            compTwoSampl.append(1)
        else:
            compTwoSampl.append(0)


        #extract variant type from vcf file
        varinatType.append(row['INFO'].split(';')[-1][5:])

        #extract depth information
        gt=row['FORMAT'].split(':').index('GT')
        dp=row['FORMAT'].split(':').index('DP')
        ao=row['FORMAT'].split(':').index('AO')
        ro=row['FORMAT'].split(':').index('RO')

        #find genotype and depth at site of variant
        GTof1st.append(row['normal'].split(':')[gt]+'.')
        GTof2nd.append(row['vaf5'].split(':')[gt]+'.')
        DPof1st.append(row['normal'].split(':')[dp])
        DPof2nd.append(row['vaf5'].split(':')[dp])

        #calculate supporting read number and ratio of reference observation and alternative oberservation

        AOread1st = 0
        for i in row['normal'].split(':')[ao].split(','):
            AOread1st = AOread1st + int(i)

        AOread2nd = 0
        for i in row['vaf5'].split(':')[ao].split(','):
            AOread2nd = AOread2nd + int(i)
            
        ReadSupport1st.append(AOread1st+int(row['normal'].split(':')[ro]))
        ratioAO1st.append(int(row['normal'].split(':')[ao].split(',')[index])/(AOread1st+int(row['normal'].split(':')[ro]))*100)
        ratioRO1st.append(int(row['normal'].split(':')[ro])/(AOread1st + int(row['normal'].split(':')[ro]))*100)


        ReadSupport2nd.append(AOread2nd+int(row['vaf5'].split(':')[ro]))
        ratioAO2nd.append(int(row['vaf5'].split(':')[ao].split(',')[index])/(AOread2nd+int(row['vaf5'].split(':')[ro]))*100)
        ratioRO2nd.append(int(row['vaf5'].split(':')[ro])/(AOread2nd+int(row['vaf5'].split(':')[ro]))*100)


        #setup ExAC query using most deleterious variant        
        urlString = 'http://exac.hms.harvard.edu/rest/variant/variant/'

        chrom = row['#CHROM']
        pos = row['POS']
        ref = row['REF']
        alt = row['ALT'].split(',')[index]
        urlString = urlString+chrom+'-'+pos+'-'+ref+'-'+alt
        exacInfo = ExacQuery(urlString)

        
        alleleFrequency.append(exacInfo["allele_frequency"])
        rsID.append(exacInfo["rsID"])
        vepAnnotationFirst.append(IteVEP(exacInfo["vepAnnotationFirst"]))
        vepAnnotationCanonical.append(IteVEP(exacInfo["vepAnnotationCanonical"]))

        GENESYM.append(exacInfo["GENESYM"])
        CCDS.append(exacInfo["CCDS"])
        MOTIF_NAME.append(exacInfo["MOTIF_NAME"])
        Conquecence.append(exacInfo["Conquecence"])
        PROTEINCHANGE.append(exacInfo["PROTEINCHANGE"])
        CLIN_SIG.append(exacInfo["CLIN_SIG"])
        PUBMED.append(exacInfo["PUBMED"])
        ENSP.append(exacInfo["ENSP"])
        CDC_postion.append(exacInfo["CDC_postion"])
        INTRON.append(exacInfo["INTRON"])
        EXON.append(exacInfo["EXON"])
        STRAND.append(exacInfo["STRAND"])
        SIFT.append(exacInfo["SIFT"])
        PolyPhen.append(exacInfo["PolyPhen"])

        #show progress information, how many variants are queried.
        n+=1
        proInfor = str(n) + ' variants have been queried.'
        sys.stdout.write('\r' + proInfor)        

        
    vcfDF['compTwoSampl'] = compTwoSampl
    vcfDF['GTof1st'] = GTof1st
    vcfDF['DPof1st'] = DPof1st
    vcfDF['ReadSupport1st'] = ReadSupport1st
    vcfDF['ratioAO1st (%)'] = ratioAO1st
    vcfDF['ratioRO1st (%)'] = ratioRO1st
    vcfDF['GTof2nd'] = GTof2nd
    vcfDF['DPof2nd'] = DPof2nd
    vcfDF['ReadSupport2nd'] = ReadSupport2nd
    vcfDF['ratioAO2nd (%)'] = ratioAO2nd
    vcfDF['ratioRO2nd (%)'] = ratioRO2nd


    vcfDF['alleleFrequency'] = alleleFrequency
    vcfDF['rsID'] = rsID
    vcfDF['vepAnnotationMostImpactTranscripts'] = vepAnnotationFirst
    vcfDF['vepAnnotationCanonicalTranscript'] = vepAnnotationCanonical

    vcfDF['Gene Symbol'] = GENESYM
    vcfDF['CCDS ID'] = CCDS
    vcfDF['Motif Name'] = MOTIF_NAME
    vcfDF['Mutation type'] = Conquecence
    vcfDF['Amino acid change'] = PROTEINCHANGE
    vcfDF['Clincal Significance'] = CLIN_SIG
    vcfDF['PubMed'] = PUBMED
    vcfDF['ENSP ID'] = ENSP
    vcfDF['CDC position'] = CDC_postion
    vcfDF['Intron'] = INTRON
    vcfDF['Exon'] = EXON
    vcfDF['STRAND'] = STRAND
    vcfDF['SIFT'] = SIFT
    vcfDF['PolyPhen'] = PolyPhen
    vcfDF['VariantType'] = varinatType


def AnnoFile(fileName):
    fileName = fileName
    processInfo = VCFtoDF(fileName)[0]
    metaInfo = VCFtoDF(fileName)[1]
    vcfDF = VCFtoDF(fileName)[2]

    print(processInfo)
    if (processInfo == 'Success! File convert to dataframe.'):
        Annotation(vcfDF)
    else:
        print('Please check your file and try again!')


    vcfDF = vcfDF[['compTwoSampl','#CHROM','POS','rsID','REF','ALT','QUAL','FILTER','VariantType','Gene Symbol','CCDS ID','Motif Name',
                   'Mutation type','Amino acid change','Clincal Significance','PubMed','ENSP ID','CDC position',
                   'Intron','Exon','STRAND','SIFT','PolyPhen','GTof1st','DPof1st','ReadSupport1st','ratioAO1st (%)',
                   'ratioRO1st (%)','GTof2nd','DPof2nd', 'ReadSupport2nd','ratioAO2nd (%)','ratioRO2nd (%)',
                   'vepAnnotationMostImpactTranscripts','vepAnnotationCanonicalTranscript','INFO','FORMAT','normal','vaf5']]

    try:
        #save dataframe to temporary csv file.
        #creat new file to add metaInformation, then delete temporary file
        csvFileNamet = './CSVreport/' + fileName[15:-4] + 'AnnoCSVt.csv'
        csvFileNamef = './CSVreport/' + fileName[15:-4] + 'AnnoCSV.csv'
        vcfDF.to_csv(csvFileNamet, index=False)
        
        f = open(csvFileNamet,'r')
        newf = open(csvFileNamef,'w')
        lines = f.readlines() 
        newf.write(metaInfo) 
        for line in lines: 
            newf.write(line)
        f.close()
        newf.close()

        os.remove(csvFileNamet)

        print()
        print('Annotation to csv file completed!')
        print('File saved to '+ csvFileNamef + '!')

    except Exception as e:
        print('Save to csv failed! ')
        print(e)

#checking input file name for annotation
def main():    
    AnnotationProcess = True
            
    while (AnnotationProcess):
        try:
            vcfFileName = input("please input file name without file extension (use 'n' to quit): ")
            if(vcfFileName == 'n' ):
                print('Program ended.')
                AnnotationProcess = False
            else:                         
                vcfFileName = './OriginalFile/' + vcfFileName + '.vcf'
                exists = os.path.isfile(vcfFileName)

                if (exists):
                    AnnoFile(vcfFileName)
                    AnnotationProcess = False
                    
                else:
                    print('Please check file name and check file is in OriginalFile folder!')
                

          
            
        except Exception as e:
            processInfo = 'Error! ' + str(e)




main()








    



    
    



    
