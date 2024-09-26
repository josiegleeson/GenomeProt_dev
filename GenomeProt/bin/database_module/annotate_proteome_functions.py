from pycdhit import cd_hit, read_clstr
import peptides as pep



class clusterSequences():
    def SeqClust(self,fasta,out):
      res = cd_hit(i=fasta,o=out,c=1.0,d=0,sc=1,n=3,G=0,aS=1,aL=0,A=0)
      return res
class refDb():
    def getDbAnnotations(self,db,openprot_map,ref_prot_map,uniprot_map): #Reference protein sequence database
      fh=open(db)
      for i in fh:
        if "#" not in i.strip():
          id=i.strip().split("\t")[0].split("|")[0].replace(">",'')
          seq=i.strip().split("\t")[1].strip()
          uniprot_header=i.strip().split("\t")[2].strip()
          
          openprot_map[seq]=id
          if '>IP_' not in i.strip()  and '>II_' not in i.strip():
            ref_prot_map[seq]=id
          if uniprot_header !="-":
            uniprot_accession=i.strip().split("\t")[2].split("|")[0].replace(">","")+"|"+i.strip().split("\t")[2].split("|")[1] #trEMBL or Reviewed| Accession
            uniprot_gene_description=i.strip().split("\t")[3].strip()
            uniprot_map[seq]=uniprot_accession+"|"+uniprot_gene_description
      fh.close()

            
    def seqRefAnnotations(self,orfmap,seq):
        if seq in orfmap.keys():
            return orfmap[seq]

class Annotations():
    def UTRAnnotations(self,tr,utrmap,orf,st,trcds): #transcript,utrcord, orfcoord,strand,transcript coordinates
        if ":" in orf and "-" in orf:
          orfBT=""
          ochr=orf.split(':')[0]
          utr_anno={} #annotates 5' and 3' UTR coordinates: k: UTR coordinates v:3UTR/5UTR 
          utr_orf_pos=[]
          #segreagate 3UTR and 5UTR coordinates and determine postion of ORFs in UTR exons
          if st=="+":
            ostart=int(orf.split(":")[1].split('-')[0])
            oend=int(orf.split(":")[1].split('-')[1])
            first_cds=int(trcds[0].split(':')[1].split("-")[0]) #start position of first CDs of transcript
            last_cds=int(trcds[-1].split(":")[1].split("-")[1]) # end position of last CDS of transcript
            for utr in utrmap:  #Identifying 5' UTR and 3'UTR exons among UTR cooordinates of transcript
              uchr=utr.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utr.split(':')[1].split('-')[0])
                u_end=int(utr.split(':')[1].split('-')[1])
                #marking 5UTR and 3UTR 
                if u_start < first_cds and u_end < first_cds:#if UTR start and end is lesser than first cds start then it is 5UTR
                  utr_anno[utr]='5UTR'
                elif u_start > last_cds and u_end>last_cds: #if UTR start and end is greater than first cds then it is 5UTR
                  utr_anno[utr]='3UTR'
            
            for utrs in utrmap:
              uchr=utrs.split(':')[0].strip()
              u_start=int(utrs.split(':')[1].split('-')[0])
              u_end=int(utrs.split(':')[1].split('-')[1])
              if uchr==ochr:
                if ostart >=u_start and ostart < u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("*"+utr_anno[utrs]) # * indicates if ORF starts in 5' or 3' UTR
                elif ostart >u_start and ostart <= u_end: 
                  utr_orf_pos.append("*"+utr_anno[utrs]) 
                  
                  
                if oend >= u_start and oend< u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR
                elif oend > u_start and oend<= u_end: #if ORF start is greater than UTR start and less than than UTR end
                  utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR
          
          if st=="-":
            ostart=int(orf.split(":")[1].split('-')[1])
            oend=int(orf.split(":")[1].split('-')[0])
            first_cds=int(trcds[0].split(':')[1].split("-")[1])
            last_cds=int(trcds[-1].split(":")[1].split("-")[0])
            for utr in utrmap:
              uchr=utr.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utr.split(':')[1].split('-')[1])
                u_end=int(utr.split(':')[1].split('-')[0])
                if u_start > first_cds and u_end > first_cds :
                  utr_anno[utr]='5UTR'
                elif u_start < last_cds and u_end<last_cds:
                  utr_anno[utr]='3UTR'
                

            for utrs in utrmap:
              uchr=utrs.split(':')[0].strip()
              if uchr==ochr:
                u_start=int(utrs.split(':')[1].split('-')[1])
                u_end=int(utrs.split(':')[1].split('-')[0])
                #print(ostart,oend,"***",u_start,u_end)
                #start
                if ostart <= u_start and ostart > u_end: #if ORF start is lesser than UTR start and greater than than UTR end
                  utr_orf_pos.append("*"+utr_anno[utrs]) # * indicates if ORF starts in 5' or 3' UTR
                elif ostart < u_start and ostart >= u_end: 
                  utr_orf_pos.append("*"+utr_anno[utrs])
                  #print(ostart,oend,"<==>",u_start,u_end)
                #end
                if oend<= u_start and oend > u_end: #if ORF end is lesser than UTR start and greater than than UTR end
                    utr_orf_pos.append("**"+utr_anno[utrs]) # ** indicates if ORF ends in 5' or 3' UTR
                    #print(ostart,oend,"<==>",u_start,u_end)
                elif oend< u_start and oend >= u_end:
                  utr_orf_pos.append("**"+utr_anno[utrs])
                
          ######Determining ORF type based on utr_orf_pos
          #print(tr,utr_orf_pos,orf)
          #print(tr, utr_anno.items(),"cds",trcds)

          if len(utr_orf_pos)==1:
            if utr_orf_pos[0]=='*5UTR':
              orfBT="5UTR:CDS"
            elif utr_orf_pos[0]=='*3UTR':
              orfBT="3UTR"
            elif utr_orf_pos[0]=='**3UTR':
              orfBT="CDS:3UTR" #Alt-CDS:3UTR
            elif utr_orf_pos[0]=='**5UTR': #handeling ambiguity due to ORFs locations predicted by ORFik are off by certain nucleiotide differences
              orfBT="5UTR"
            
              
              
          elif len(utr_orf_pos)==2:
            if utr_orf_pos[0]=='*3UTR' and utr_orf_pos[1]=='**3UTR':
              orfBT="3UTR"
            elif utr_orf_pos[0]=='*5UTR' and utr_orf_pos[1]=='**5UTR':
              orfBT="5UTR"
            elif utr_orf_pos[0]=='*5UTR' and utr_orf_pos[1]=='**3UTR':
              orfBT="5UTR:3UTR"
            elif utr_orf_pos[0]=='**3UTR' and utr_orf_pos[1]=='*5UTR':
              orfBT="5UTR:3UTR"
          return orfBT
        
    def isIntergenic(self,orf_cord,gene_cord_map):
      #ORF
      if ":" in orf_cord and "-" in orf_cord:
        ochr=orf_cord.split(":")[0]
        ostart=int(orf_cord.split(":")[1].split("-")[0])
        oend=int(orf_cord.split(":")[1].split("-")[1])
        orf_overlap=""
        for gnpost in gene_cord_map.values():
          #gene
          genchr=gnpost.split(":")[0]
          gnstart=int(gnpost.split(":")[1].split('-')[0])
          gnend=int(gnpost.split(":")[1].split('-')[1])
          if ochr==genchr:
            if ostart >= gnstart and oend <= gnend:
              orf_overlap="gene_overlap"
              break
        return orf_overlap
      
class SequenceProperties():
  def calculateMolWt(self,seq):
    peptide = pep.Peptide(seq)
    mol_weight=round(float(peptide.molecular_weight())/1000,2) #Molecular weight in KDa
    return mol_weight
  def calculateIsoElectricPoint(self,seq):
    peptide = pep.Peptide(seq)
    Iso_ele_point=round(peptide.isoelectric_point(),2) #Isoelectric point
    return Iso_ele_point
  def calculateHydrophobicity(self,seq):
    peptide = pep.Peptide(seq)
    hydrophobicity=round(peptide.hydrophobicity(),2) #hydrophobicity
    return hydrophobicity
  def calculateAliphatic_index(self,seq):
    peptide = pep.Peptide(seq)
    aliphatic_index=round(peptide.aliphatic_index(),2) #hydrophobicity
    return aliphatic_index
