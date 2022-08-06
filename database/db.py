#!/usr/bin/python
import sys
# get sys package for file arguments etc
import mysql.connector
from mysql.connector import Error
import pandas as pd
from pymysql.converters import escape_string
try:
    con = mysql.connector.connect(host='localhost',
    database='',user='',password='') # input database name, user name, and password
    if con.is_connected():
        print(11)
        cursor = con.cursor() # buffered=True
        cursor.execute("select database();")
        record = cursor.fetchone()
        print("You're connected to database: ", record)
except Error as e:
    print("Error while connecting to MySQL", e)

project = pd.read_csv('/datastore/home/s2172876/project/recount3_bulkRNA_project_description.csv',
            usecols=['project','organism', 'project_home', 'n_samples', 'study_title', 'study_abstract'])
project = project[['project','organism', 'project_home', 'n_samples', 'study_title', 'study_abstract']]
project['n_samples'] = project['n_samples'].astype(str)
project.fillna('', inplace=True)
project['project_home'] = project['project_home'].str[13:]
sample = pd.read_csv('/datastore/home/s2172876/project/sample_annotation.csv')
sample.fillna('', inplace=True)
summary = pd.read_csv('/datastore/home/s2172876/project/DEG_summary_fix.txt',sep='\t')
group = pd.read_csv('/datastore/home/s2172876/project/group.txt',sep = '\t')
group.fillna('', inplace=True)
TCGA = pd.read_csv('/datastore/home/s2172876/project/TCGA_annotation.csv')
TCGA.fillna('', inplace=True)
TCGA = TCGA.applymap(lambda x: x.replace('\'', ''))
GTEx = pd.read_csv('/datastore/home/s2172876/project/gtex_annotation.csv',sep=' ')
GTEx.fillna('', inplace=True)
motif = pd.read_csv('/datastore/home/s2172876/project/motif_summary_agg.txt',sep = '\t')
motif = motif.sort_values(by=['count'], ascending=False)
mirna = pd.read_csv('/datastore/home/s2172876/project/mirna_info.txt',sep='\t')
mirna = mirna[(mirna['organism'] == 'Homo sapiens') | (mirna['organism'] == 'Mus musculus')]
mirna = mirna[['mir','id','organism','seq_aa','seed','motif']]

cursor.execute("CREATE TABLE project (\
`projectID`  VARCHAR(80) NOT NULL,\
`organism` VARCHAR(80) NOT NULL,\
`project_home` VARCHAR(80),\
`n_samples`  VARCHAR(80),\
`study_title`VARCHAR(10000),\
`study_abstract` VARCHAR(1000) )")
for i in range(0,project.shape[0]):
    sql = "INSERT INTO project (projectID, organism, project_home,\
            n_samples,study_title,study_abstract) \
            VALUES ('%s','%s','%s','%s','%s','%s')" % \
                (project.iloc[i,0],project.iloc[i,1],\
                    project.iloc[i,2],project.iloc[i,3],\
                        escape_string(project.iloc[i,4]),escape_string(project.iloc[i,5]))
    cursor.execute(sql)
    print(i)

con.commit()


cursor.execute("CREATE TABLE sample (\
`sampleID`  VARCHAR(80) NOT NULL,\
`projectID`  VARCHAR(80) NOT NULL,\
`tissue` VARCHAR(80) ,\
`cellline` VARCHAR(80),\
`celltype`  VARCHAR(80),\
`samplecondition` VARCHAR(1000),\
`treatment` VARCHAR(1000), \
`treatment_duration` VARCHAR(1000),\
`genotype` VARCHAR(100),\
`subtype` VARCHAR(1000))")

for i in range(0,sample.shape[0]):
    sql = "INSERT INTO sample (sampleID, projectID, tissue,cellline,celltype,samplecondition,\
        treatment,treatment_duration,genotype,subtype) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
        % (sample.iloc[i,1],sample.iloc[i,0],escape_string(sample.iloc[i,2]),escape_string(sample.iloc[i,3]),escape_string(sample.iloc[i,4]),escape_string(sample.iloc[i,5]),\
            escape_string(sample.iloc[i,6]),escape_string(sample.iloc[i,7]),escape_string(sample.iloc[i,8]),escape_string(sample.iloc[i,9]))
    cursor.execute(sql)
    print(i)

con.commit()
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir('/datastore/home/s2172876/project/DEG/') if isfile(join('/datastore/home/s2172876/project/DEG/', f))]
DEG_Result = pd.read_csv('/datastore/home/s2172876/project/DEG/'+ onlyfiles[0],sep = '\t',)
DEG_Result = DEG_Result.rename_axis("Gene").reset_index()
DEG_Result = DEG_Result[(abs(DEG_Result['log2FoldChange']) > 1) & (DEG_Result['padj']<0.05)]
DEG_Result['ExpID'] = str((onlyfiles[0][0:-3]))
for i in range(1,len(onlyfiles)):
    name = str((onlyfiles[i][0:-3]))
    globals()[name] = pd.read_csv('/datastore/home/s2172876/project/DEG/'+ onlyfiles[i],sep = '\t',)
    globals()[name] = globals()[name].rename_axis("Gene").reset_index()
    globals()[name] = globals()[name][(abs(globals()[name]['log2FoldChange']) > 1) & ( globals()[name]['padj']<0.05)]
    globals()[name]['ExpID'] = name
    DEG_Result = pd.concat([DEG_Result, globals()[name]])
    print(i,globals()[name].shape[0])


cursor.execute("CREATE TABLE DEGResult (\
`Gene`  VARCHAR(80) NOT NULL,\
`log2FoldChange`  DOUBLE,\
`padj` DOUBLE ,\
`ExpID` VARCHAR(80))")
# DEG_Result = pd.read_csv('/datastore/home/s2172876/project/allDEG.csv')
# DEG_Result = DEG_Result[['Gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat',
#        'pvalue', 'padj', 'ExpID']]
for i in range(30310504,DEG_Result.shape[0]):
    sql = "INSERT INTO DEGResult (Gene, log2FoldChange, padj,ExpID) VALUES ('%s',%f,%f,'%s')" \
        % (DEG_Result.iloc[i,0],DEG_Result.iloc[i,2],DEG_Result.iloc[i,6],DEG_Result.iloc[i,7])
    cursor.execute(sql)
    con.commit()
    print(i)

cursor.execute("CREATE TABLE DEGSummary (\
`ExpID`  VARCHAR(80) NOT NULL,\
`ProjectID`  VARCHAR(80) NOT NULL,\
`Comparison` VARCHAR(1000) ,\
`DEG` DOUBLE,\
`Up`  DOUBLE,\
`Down` DOUBLE,\
`non_DEG` DOUBLE, \
`Group1` VARCHAR(100000),\
`Group2` VARCHAR(10000))")

for i in range(0,summary.shape[0]):
    sql = "INSERT INTO DEGSummary (ExpID, ProjectID, Comparison, DEG, Up, Down, non_DEG, Group1, Group2)\
         VALUES ('%s','%s','%s', %d, %d, %d, %d,'%s','%s')" \
        % (summary.iloc[i,0],summary.iloc[i,1],escape_string(summary.iloc[i,2]),summary.iloc[i,3],\
        summary.iloc[i,4],summary.iloc[i,5],summary.iloc[i,6],summary.iloc[i,7],summary.iloc[i,8])
    cursor.execute(sql)
    print(i)

con.commit()

cursor.execute("CREATE TABLE samplegroup (\
`ProjectID`  VARCHAR(80) NOT NULL,\
`tissue`  VARCHAR(80) NOT NULL,\
`cellline` VARCHAR(100) ,\
`celltype` VARCHAR(100),\
`samplecondition`  VARCHAR(1000),\
`treatment` VARCHAR(1000),\
`treatment_duration` VARCHAR(1000), \
`genotype` VARCHAR(1000),\
`subtype` VARCHAR(1000),\
`external_id` VARCHAR(100000))")

for i in range(0,group.shape[0]):
    sql = "INSERT INTO samplegroup (ProjectID, tissue, cellline, celltype, samplecondition,\
         treatment, treatment_duration, genotype, subtype, external_id)\
         VALUES ('%s','%s','%s', '%s', '%s', '%s', '%s','%s','%s','%s')" \
        % (group.iloc[i,0],escape_string(group.iloc[i,1]),escape_string(group.iloc[i,2]),escape_string(group.iloc[i,3]),\
        escape_string(group.iloc[i,4]),escape_string(group.iloc[i,5]),escape_string(group.iloc[i,6]),escape_string(group.iloc[i,7]),\
            escape_string(group.iloc[i,8]),group.iloc[i,9])
    cursor.execute(sql)
    print(i)
con.commit()


cursor.execute("CREATE TABLE TCGAannotation (`study`  VARCHAR(80) NOT NULL,\
`external_id`  VARCHAR(80) NOT NULL,\
`project_name` VARCHAR(100) ,\
`sample_type` VARCHAR(100),\
`tumor_stage`  VARCHAR(1000),\
`diagnoses_vital_status` VARCHAR(1000),\
`case_sample_type` VARCHAR(1000),\
`tumor_status` VARCHAR(1000),\
`case_vital_status` VARCHAR(1000),\
`pathologic_n` VARCHAR(1000),\
`clinical_m` VARCHAR(1000),\
`clinical_stage` VARCHAR(1000),\
`pathologic_stage` VARCHAR(1000),\
`pathologic_t` VARCHAR(1000),\
`histological_type` VARCHAR(1000),\
`submitter_id` VARCHAR(1000))")

for i in range(0,TCGA.shape[0]):
    sql = "INSERT INTO TCGAannotation (study, external_id, project_name, sample_type, tumor_stage,diagnoses_vital_status,\
         case_sample_type, tumor_status, case_vital_status, pathologic_n, clinical_m, clinical_stage, pathologic_stage,\
         pathologic_t,histological_type,submitter_id)\
         VALUES ('%s','%s','%s', '%s', '%s', '%s', '%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" \
        % (TCGA.iloc[i,0],TCGA.iloc[i,1],TCGA.iloc[i,2],TCGA.iloc[i,3],\
        TCGA.iloc[i,4],TCGA.iloc[i,5],TCGA.iloc[i,6],TCGA.iloc[i,7],\
        TCGA.iloc[i,8],TCGA.iloc[i,9],TCGA.iloc[i,10],TCGA.iloc[i,11],\
        TCGA.iloc[i,12],TCGA.iloc[i,13],TCGA.iloc[i,14],TCGA.iloc[i,15])
    cursor.execute(sql)
    print(i)
con.commit()
cursor.execute("DELETE FROM TCGAannotation where study='';")
con.commit()

cursor.execute("CREATE TABLE GTEx (`external_id`  VARCHAR(80) NOT NULL,\
`run_acc`  VARCHAR(80) NOT NULL,\
`study` VARCHAR(100) ,\
`SUBJID` VARCHAR(100),\
`SEX`  VARCHAR(1000),\
`AGE` VARCHAR(1000),\
`SAMPID` VARCHAR(1000))")
for i in range(0,GTEx.shape[0]):
    sql = "INSERT INTO GTEx (external_id, run_acc, study, SUBJID, SEX, AGE, SAMPID)\
         VALUES ('%s','%s','%s', '%s', '%s', '%s', '%s')" \
        % (GTEx.iloc[i,0],GTEx.iloc[i,1],GTEx.iloc[i,2],GTEx.iloc[i,3],\
        GTEx.iloc[i,4],GTEx.iloc[i,5],GTEx.iloc[i,6])
    cursor.execute(sql)
    print(i)
con.commit()

cursor.execute("CREATE TABLE motif (`motif`  VARCHAR(80) NOT NULL,\
`type`  VARCHAR(100) NOT NULL,\
`organism`  VARCHAR(100) NOT NULL,\
`ExpID`  VARCHAR(100000) NOT NULL,\
`count` DOUBLE)")

for i in range(0,motif.shape[0]):
    sql = "INSERT INTO motif (motif, type, organism, ExpID, count)\
         VALUES ('%s','%s','%s','%s',%d)" \
        % (motif.iloc[i,0],motif.iloc[i,1],motif.iloc[i,2],motif.iloc[i,3],motif.iloc[i,4])
    cursor.execute(sql)
    print(i)
con.commit()

cursor.execute("CREATE TABLE mirna (`mir`  VARCHAR(80) NOT NULL,\
`id`  VARCHAR(80) NOT NULL,\
`organism`  VARCHAR(80),\
`seq_aa`  VARCHAR(500),\
`seed`  VARCHAR(80),\
`motif`  VARCHAR(80))")

for i in range(0,mirna.shape[0]):
    sql = "INSERT INTO mirna (mir, id, organism,seq_aa,seed,motif)\
         VALUES ('%s','%s','%s','%s','%s','%s')" \
        % (mirna.iloc[i,0],escape_string(mirna.iloc[i,1]),mirna.iloc[i,2],mirna.iloc[i,3],mirna.iloc[i,4],mirna.iloc[i,5])
    cursor.execute(sql)
    print(i)
con.commit()

# from pickle import TRUE
# from sqlalchemy import create_engine
# import mysql.connector
# url = 'mysql+pymysql://s2172876:123456@localhost/s2172876'
# engine = create_engine(url, echo=True)
# # connection = engine.connect()
# connection = engine.raw_connection()
# cursor = connection.cursor()
# project.to_sql('project',con = con)
# print(cur)
