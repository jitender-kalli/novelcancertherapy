import numpy as np
import pandas as pd
import os


#splits CSV into multiple CSVs based on row length
def split(filehandler, delimiter=',', row_limit=499,
          output_name_template='output_%s.csv', output_path='.', keep_headers=True):
    import csv
    reader = csv.reader(filehandler, delimiter=delimiter)
    current_piece = 1
    current_out_path = os.path.join(
        output_path,
        output_name_template % current_piece
    )
    current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
    current_limit = row_limit
    if keep_headers:
        headers =next(reader)
        current_out_writer.writerow(headers)
    for i, row in enumerate(reader):
        if i + 1 > current_limit:
            current_piece += 1
            current_limit = row_limit * current_piece
            current_out_path = os.path.join(
                output_path,
                output_name_template % current_piece
            )
            current_out_writer = csv.writer(open(current_out_path, 'w'), delimiter=delimiter)
            if keep_headers:
                current_out_writer.writerow(headers)
        current_out_writer.writerow(row)



#counts number of sequences in fasta file
def countfasta(filename):
    f = open(filename, "r");
    count=0;
    for x in f:
        for y in x:
            #print (y)
            if y==">":
                count+=1
    return(count);

#counts number of sequences in IEDB output file
def getiedbcount(iedbfilename):
    df = pd.read_csv(iedbfilename, sep='\t')
    df=df.sort_values(by=['seq_num'])
    last=df["seq_num"].iloc[-1];
    return (last);

#Lists all tissues from the csv file (of RawFile type), modify as required
def gettissuelist(rawfilename):
    df=pd.read_csv(rawfilename,sep = '\t',engine = 'python');
    listx=df["Tissue"].unique();
    return (listx);

#processes fasta file and makes all changes so that file is ready for IEDB processing.
def processfasta():
    df = pd.read_csv("output.csv", sep='\t', header=None)
    df.columns = ["Gene name", "Sequence"]
    df.to_csv('output1.csv') ;
    listy=df["Gene name"].tolist()

    listx= df["Sequence"].tolist()

    #list numbers contains all instances of X
    listnumbersx=[];
    listnumbersu=[]
    for i in range(len(listx)) :
      sequence=str(listx[i]);
      for aminoacid in sequence:
        if aminoacid == "X":
            listnumbersx.append(i);
            break;
        if aminoacid == "U":
            listnumbersu.append(i);
            break;


    f = open("modified.txt", "w")
    f.write("Modified X")

    for genenumbers in listnumbersx:
      #print(listy[genenumbers]);
      f.write(str(listy[genenumbers])+"\n");
      
    f.close()

    f=open("modified.txt", "a+")
    f.write("Modified U")
    for genenumbers in listnumbersu:
      #print(listy[genenumbers]);
      f.write(str(listy[genenumbers])+"\n");
      
    f.close()

    #numbers=[];
    #for i in range(len(listx)) :
      #sequence=str(listx[i]);
      #if sequence == 'Sequence unavailable':
        #numbers.append(i)

    #print (numbers);
    #for i in numbers:
      #df=df.drop(df.index[i])

    df['Sequence']= df['Sequence'].replace('[X]+', 'Q', regex=True)
    df['Sequence']= df['Sequence'].replace('[U]+', 'Q', regex=True)
    #df['Sequence']= df['Sequence'].replace('Sequence\sunavailable', np.nan, regex=True)
    df['Sequence']= df['Sequence'].replace('Sequence\sunavailable', np.nan, regex=True)
    #df.dropna(inplace= True)

    df['Sequence']= df['Sequence'].replace('[*]+', '', regex=True)
    df.dropna(inplace= True)
    df=df[df['Sequence'].map(len) > 9]
    df.to_csv('filex.fasta', sep="\n" ,header=False, index=False) ;

#converts fasta to csv
def fasta2csv(fastafilename):
    count=0;
    header="";
    seq="";
    f=open(fastafilename,"r");
    f1 = open("op.csv", "a");
    for line in f:
        line=line.replace("\n", "")
        if str(line).startswith(">"):
            if count == 0:
                header=line+"\t";
               # print(header);
                #f1.write(line+"\t"); 
                count+=1;
                
            else:
                print(header+seq);
                f1.write(header+seq);
                seq="";
                header="\n"+line+"\t";
               # print(header);
                #f1.write(" "+line+"\t");
                count+=1;
                
        else :
            seq+=line;
            #f1.write(line);
           # print (line);
    f1.write(header+seq);

#makes new files with 500 genes each from the tissue list provided    
def filtergenes(tissuelist):

    df =pd.read_csv('normal_tissue.tsv',sep='\t' ,engine='python')
    df=df[["Gene","Tissue"]]
    #for i in tissuelist:    
        #print (i,i in df['Tissue'].values)
    
    df=df[df['Tissue'].isin(tissuelist)]
    df= df['Gene']
    df.drop_duplicates(keep='first',inplace=True) 
    df.to_csv('filtered_normal_tissue.csv', sep="\n" ,header=False, index=False) ;
    split(open('filtered_normal_tissue.csv', 'r'));
    
    

    
    return("new files and filtered_normal_tissue.tsv created ")


tissuelist=["breast","bronchus","bone marrow",'cerebellum','cerebral cortex','cervix, uterine','colon','duodenum','epididymis','esophagus','fallopian tube','gallbladder','heart muscle','kidney','liver','lung','ovary','pancreas','placenta','prostate','skeletal muscle','spleen','stomach 1','stomach 2','thyroid gland','urinary bladder','pituitary gland']
print(filtergenes(tissuelist))








