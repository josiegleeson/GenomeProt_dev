import re
import sys
import os
#from  FileHandelingUtilities import  *
import operator
from collections import Counter


#############GTF FILE PARSER###########################################
# Script By: Hitesh Kore

#Usage: python GTFFileParser.py <GTF file>
#######################################################################

class GTFParser():
    def __init__(self,line):
      self.rec=line.strip()
      if len(self.rec.split("\t"))>2:
        self.chr=self.rec.split("\t")[0]
        self.feature = self.rec.split("\t")[2]
        self.start = int(self.rec.split("\t")[3])
        self.end = int(self.rec.split("\t")[4])
        self.strand = self.rec.split("\t")[6]
        self.meta=self.rec.split("\t")[8].strip().split(";")
        self.ensgene=",".join([i.replace("gene_id","").replace('"',"") for i in self.meta if "gene_id" in i]).strip()
        self.transcriptid=",".join([i.replace("transcript_id","").replace('"',"") for i in self.meta if "transcript_id" in i]).strip()
        self.supportevi=",".join([i.replace("transcript_support_level","").replace('"',"") for i in self.meta if "transcript_support_level" in i]).strip()
        self.tag=",".join([i.replace("tag","").replace('"',"").strip() for i in self.meta if "tag" in i]).strip()
        self.havanagene =",".join([i.replace("havana_gene","").replace('"',"").strip() for i in self.meta if "havana_gene" in i]).strip()
        self.havanatranscript = ",".join([i.replace("havana_transcript","").replace('"',"").strip() for i in self.meta if "havana_transcript" in i]).strip()
        self.genename = ",".join([i.replace("gene_name","").replace('"',"").strip() for i in self.meta if "gene_name" in i]).strip()
        self.genetype = ",".join([i.replace("gene_type","").replace('"',"").strip() for i in self.meta if "gene_type" in i]).strip()
        self.trascript_name = ",".join([i.replace("transcript_name","").replace('"',"").strip() for i in self.meta if "transcript_name" in i]).strip()
        self.ccds_id = ",".join([i.replace("ccdsid","").replace('"',"").strip() for i in self.meta if "ccdsid" in i]).strip()
        self.exn_no = ",".join([i.replace("exon_number","").replace('"',"").strip() for i in self.meta if "exon_number" in i]).strip()
        self.exn_id = ",".join([i.replace("exon_id","").replace('"',"").strip() for i in self.meta if "exon_id" in i]).strip()
        self.transcript_type = ",".join([i.replace("transcript_type","").replace('"',"").strip() for i in self.meta if "transcript_type" in i]).strip()




    def getCoordinates(self):
       return self.chr+":"+str(self.start)+"-"+str(self.end)
    def featureExists(self,feature):
        fg=False
        if feature.strip() in self.feature:
            fg=True
        return fg

    def getFeatureLength(self,end,start):
        return int(end)-int(start)

    def fieldsCheck(self):
        fg=False
        if len(self.rec.split("\t"))==9:
            fg=True
        return fg
    def getTranscriptTypeCount(self,trid,trtype,trmap):
        trmap.setdefault(trtype,[]).append(trid)

class RefSeqGTFParser():
    def __init__(self,line):
        self.rec=line.strip()
        self.chr=self.rec.split("\t")[0]
        self.feature = self.rec.split("\t")[2]
        self.start = int(self.rec.split("\t")[3])
        self.end = int(self.rec.split("\t")[4])
        self.strand = self.rec.split("\t")[6]
        self.meta=self.rec.split("\t")[8].strip().split(";")
        self.gene=",".join([i.replace("gene_id","").replace('"',"") for i in self.meta if "gene_id" in i]).strip()
        self.transcriptid=",".join([i.replace("transcript_id","").replace('"',"") for i in self.meta if "transcript_id" in i]).strip()
        self.transcript_type = ",".join([i.replace("transcript_biotype","").replace('"',"").strip() for i in self.meta if "transcript_biotype" in i]).strip()
        self.genetype = ",".join([i.replace("gene_biotype","").replace('"',"").strip() for i in self.meta if "gene_biotype" in i]).strip()

    def getCoordinates(self):
       return self.chr+":"+str(self.start)+"-"+str(self.end)
    def featureExists(self,feature):
        fg=False
        if feature.strip() in self.feature:
            fg=True
        return fg

    def getFeatureLength(self,end,start):
        return int(end)-int(start)

    def fieldsCheck(self):
        fg=False
        if len(self.rec.split("\t"))==9:
            fg=True
        return fg
    def getTranscriptTypeCount(self,trid,trtype,trmap):
        trmap.setdefault(trtype,[]).append(trid)


