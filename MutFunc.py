#!/usr/bin/env python3

import os
import sys
import subprocess
from os.path import basename

sep=' '
dir1='1_unmapped'
dir2='2_sort'
dir3='3_dedup'
dir4='4_realign'
dir5='5_recal'
dira='annotate'
dirc='cnv'
dirf='fusion'
dirs='snv_indel'
dirt='tmb'
dirm='msi'
split_chr={
  4:[
    ['chr2','chrM','chrY','chr18'],
    ['chr1','chr7','chr5','chr16','chr8','chr14'],
    ['chr19','chr6','chr3','chr4','chr10','chr22','chr13'],
    ['chr17','chr12','chr9','chr11','chr20','chrX','chr15','chr21']],
  12:[
    ['chr1','chrY'],
    ['chr2','chrM'],
    ['chr19','chr15'],
    ['chr17','chr13'],
    ['chr12','chr21'],
    ['chr6','chr14'],
    ['chr7','chr22'],
    ['chr5','chrX'],
    ['chr3','chr20','chr18'],
    ['chr9','chr10'],
    ['chr11','chr8'],
    ['chr4','chr16']]
}


def MSI(n_bam,t_bam,nMark,tMark,mantis=None,msibeddir=None,ref=None,mf=None,outdir='.',name='Sample'):
  '[I] MSI analysis'  # TODO: Tumor only
  wdir=os.path.join(outdir,'variants','somatic',dirm)

  beds=os.listdir(msibeddir)
  mf.write('\n\n.MAKEFLOW CATEGORY MSI')
  for i in beds:
    log=os.path.join(wdir,i+'_msi.log')
    out=os.path.join(wdir,i+'_msi.xls')
    done_log=os.path.join(wdir,i+'_msi_done.log')
    done_cmd='echo "%s" > %s' % (name+'_'+i+'_msi_done',done_log)
    cmd=sep.join(['python',mantis,'--bedfile',msibeddir+'/'+i,'--genome',ref,'-n',n_bam,'-t',t_bam,'-o',wdir+'/'+i+'_msi.xls 2>',log,'&&',done_cmd])
    oi=sep.join([out,done_log,':',n_bam,t_bam])
    mf.write('\n%s\n\t%s\n' % (oi,cmd))
  return([out,done_log])


def TMB(vcf,mark,tmb=None,mf=None,outdir='.',name='Sample'):
  '[H] TMB analysis'
  wdir=os.path.join(outdir,'variants','somatic',dirt)

  log=os.path.join(wdir,'tmb.log')
  out=os.path.join(wdir,'tmb.xls')
  done_log=os.path.join(wdir,'tmb_done.log')
  done_cmd='echo "%s" > %s' % (name+'_tmb_done',done_log)
  cmd=sep.join(['perl',tmb,vcf,wdir+'/tmb 2>',log,'&&',done_cmd])
  oi=sep.join([out,done_log,':',vcf])
  mf.write('\n\n.MAKEFLOW CATEGORY TMB\n')
  mf.write('%s\n\t%s\n' % (oi,cmd))
  return([out,done_log])


def CNV(n_bam,t_bam,nMark,tMark,contra=None,filter=None,lr=None,rd=None,bed=None,sex=None,mf=None,outdir='.',name='Sample'):
  '[G] Detect CNVs'  # TODO: Tumor only
  wdir=os.path.join(outdir,'variants','somatic',dirc)
  mf.write('\n\n.MAKEFLOW CATEGORY CNV_CONTRA\n')

  rn_out=os.path.join(wdir,'n_rd_result')
  tmp=os.path.join(wdir,'n_rd')
  err=os.path.join(wdir,'n_rd.err')
  done_log=os.path.join(wdir,'n_rd_done.log')
  done_cmd='echo "%s" > %s' % (name+'_n_rd_done',done_log)
  r_cmd=sep.join(['perl',rd,n_bam,rn_out,tmp,bed,'2>',err]) # &> -> 2> (for ubuntu)
  cmd=sep.join([r_cmd,'&&',done_cmd])
  oi=sep.join([rn_out,done_log,':',n_bam])
  mf.write('%s\n\t%s\n' % (oi,cmd))

  rt_out=os.path.join(wdir,'t_rd_result')
  tmp=os.path.join(wdir,'t_rd')
  err=os.path.join(wdir,'t_rd.err')
  done_log=os.path.join(wdir,'t_rd_done.log')
  done_cmd='echo "%s" > %s' % (name+'_t_rd_done',done_log)
  r_cmd=sep.join(['perl',rd,t_bam,rt_out,tmp,bed,'2>',err])
  cmd=sep.join([r_cmd,'&&',done_cmd])
  oi=sep.join([rt_out,done_log,':',t_bam])
  mf.write('\n%s\n\t%s\n' % (oi,cmd))

  c_out=os.path.join(wdir,'table',name+'.CNATable.10rd.10bases.20bins.txt')
  log=os.path.join(wdir,'contra.log')
  done_log=os.path.join(wdir,'contra_done.log')
  done_cmd='echo "%s" > %s' % (name+'_contra_done',done_log)
  cmd=sep.join([contra,'--target',bed,'--test',t_bam,'--control',n_bam,'--outFolder',wdir,'--sampleName',name,'--largeDeletion --minControlRdForCall 100 --minTestRdForCall 100 --minAvgForCall 20 --minExon 50 --nomultimapped --removeDups --plot','2>',log,'&&',done_cmd])
  oi=sep.join([c_out,done_log,':',n_bam,t_bam])
  mf.write('\n%s\n\t%s\n' % (oi,cmd))

  l_out=os.path.join(wdir,'ERBB2.csv')
  f1_out=os.path.join(wdir,'all_filter.txt')
  f1_err=os.path.join(wdir,'all_filter.err')
  f2_out=os.path.join(wdir,'mycg_filter.txt')
  f2_err=os.path.join(wdir,'mycg_filter.err')
  l1_out=os.path.join(wdir,'cnv_num_all.csv')
  l1_err=os.path.join(wdir,'cnv_num_all.err')
  l2_out=os.path.join(wdir,'cnv_num_mycg.csv')
  l2_err=os.path.join(wdir,'cnv_num_mycg.err')
  f1_cmd=sep.join(['perl',filter,c_out,'all 1',rn_out,rt_out,'>',f1_out,'2>',f1_err])
  f2_cmd=sep.join(['perl',filter,c_out,'mycg 1',rn_out,rt_out,'>',f2_out,'2>',f2_err])
  l1_cmd=sep.join(['perl',lr,f1_out,sex,l1_out,l_out,'2>',l1_err])
  l2_cmd=sep.join(['perl',lr,f2_out,sex,l2_out,'2>',l2_err])
  #l1_cmd=sep.join(['perl',lr,f1_out,'>',l1_out,'2>',l1_err])
  #l2_cmd=sep.join(['perl',lr,f2_out,'>',l2_out,'2>',l2_err])
  #g_out=os.path.join(wdir,'ERBB2.csv')
  #g_cmd=sep.join(['grep ERBB2',l1_out,'>',g_out]) # always failed with code 1 with makeflow
  done_log=os.path.join(wdir,'cnv_done.log')
  done_cmd='echo "%s" > %s' % (name+'_cnv_done',done_log)
  cmd=sep.join([f1_cmd,'&&',f2_cmd,'&&',l1_cmd,'&&',l2_cmd,'&&',done_cmd])
  oi=sep.join([l_out,done_log,':',rn_out,rt_out,c_out])
  mf.write('\n%s\n\t%s\n' % (oi,cmd))
  return([l_out,done_log])


def Fusion(n_bam,t_bam,nMark,tMark,factera=None,compfusion=None,hg2bit=None,bed=None,ref=None,bedexon=None,mf=None,outdir='.',name='Sample'):
  '[F] Detect fusions'  # TODO: Tumor only
  wdir=os.path.join(outdir,'variants','somatic',dirf)
  ndir=os.path.join(wdir,'normal')
  tdir=os.path.join(wdir,'tumor')
  assert not os.makedirs(ndir,exist_ok=True),'Cannot create %s' % ndir
  assert not os.makedirs(tdir,exist_ok=True),'Cannot create %s' % tdir

  log=os.path.join(ndir,'factera.log')
  nout=os.path.join(ndir,'merge_sorted.aln.factera.fusions.txt')
  done_log=os.path.join(ndir,'factera_done.log')
  done_cmd='echo "%s" > %s' % (name+'_factera_done',done_log)
  cmd=sep.join(['perl',factera,'-o',ndir,n_bam,bedexon,hg2bit,bed,'>',log,'&&',done_cmd]) # &> -> > (for ubuntu)
  oi=sep.join([nout,done_log,':',n_bam])
  mf.write('\n\n.MAKEFLOW CATEGORY FusionFactera\n')
  mf.write('%s\n\t%s\n' % (oi,cmd))

  log=os.path.join(tdir,'factera.log')
  tout=os.path.join(tdir,'merge_sorted.aln.factera.fusions.txt')
  done_log=os.path.join(tdir,'factera_done.log')
  done_cmd='echo "%s" > %s' % (name+'_factera_done',done_log)
  cmd=sep.join(['perl',factera,'-o',tdir,t_bam,bedexon,hg2bit,bed,'>',log,'&&',done_cmd])
  oi=sep.join([tout,done_log,':',t_bam])
  mf.write('\n%s\n\t%s\n' % (oi,cmd))

  log=os.path.join(wdir,'compFactera.log')
  done_log=os.path.join(wdir,'compFactera_done.log')
  done_cmd='echo "%s" > %s' % (name+'_compFactera_done',done_log)
  cmd=sep.join(['python',compfusion,nout,tout,'5',name,wdir+'/','>',log,'&&',done_cmd])
  oi=sep.join([done_log,':',nout,tout])
  mf.write('\n%s\n\t%s\n' % (oi,cmd))
  return([done_log])


def CleanFq(normalFastq,tumorFastq,trimmomatic,mf=None,outdir='.',name='Sample'):
  '[1]'
  n_pe_fq = sep.join(normalFastq)
  n_clean_fq_1 = os.path.join(outdir,'normal',dir1,"clean_1.fq.gz")
  n_clean_fq_1_U = os.path.join(outdir,'normal',dir1,"clean_1U_1.fq.gz")
  n_clean_fq_2 = os.path.join(outdir,'normal',dir1,"clean_2.fq.gz")
  n_clean_fq_2_U = os.path.join(outdir,'normal',dir1,"clean_2U_2.fq.gz")
  n_clean_outs = sep.join([n_clean_fq_1,n_clean_fq_1_U,n_clean_fq_2,n_clean_fq_2_U])

  t_pe_fq = sep.join(tumorFastq)
  t_clean_fq_1 = os.path.join(outdir,'tumor',dir1,"clean_1.fq.gz")
  t_clean_fq_1_U = os.path.join(outdir,'tumor',dir1,"clean_1U_1.fq.gz")
  t_clean_fq_2 = os.path.join(outdir,'tumor',dir1,"clean_2.fq.gz")
  t_clean_fq_2_U = os.path.join(outdir,'tumor',dir1,"clean_2U_1.fq.gz")
  t_clean_outs = sep.join([t_clean_fq_1,t_clean_fq_1_U,t_clean_fq_2,t_clean_fq_2_U])

  adapters_path = os.path.join(os.path.dirname(trimmomatic),"adapters/TruSeq3-PE-2.fa")

  n_trim_log = os.path.join(outdir,'normal',dir1,"trim.log")
  t_trim_log = os.path.join(outdir,'tumor',dir1,"trim.log")
  n_done_log = os.path.join(outdir,'normal',dir1,"success_trim.log")
  t_done_log = os.path.join(outdir,'tumor',dir1,"success_trim.log")

  n_done_cmd = "> %s && echo \"%s\" > %s" % (n_trim_log,name+"_n_clean", n_done_log)
  t_done_cmd = "> %s && echo \"%s\" > %s" % (t_trim_log,name+"_t_clean", t_done_log)
  n_cmd = sep.join(["java -jar",trimmomatic,"PE -threads 4",n_pe_fq,n_clean_outs,"ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % adapters_path,n_done_cmd])
  t_cmd = sep.join(["java -jar",trimmomatic,"PE -threads 4",t_pe_fq,t_clean_outs,"ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % adapters_path,t_done_cmd]) # -phred33 removed, as Trimmomatic can automatically determine the encoding since version 0.32

  Un=sep.join([n_clean_fq_1_U,n_clean_fq_2_U])
  Ut=sep.join([t_clean_fq_1_U,t_clean_fq_2_U])
  U=sep.join([Un,Ut])

  n_oi=sep.join([n_clean_fq_1,n_clean_fq_2,Un,n_done_log,':',n_pe_fq])
  t_oi=sep.join([t_clean_fq_1,t_clean_fq_2,Ut,t_done_log,':',t_pe_fq])
  mf.write("\n\n.MAKEFLOW CATEGORY FastqClean\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_cmd))
  return([sep.join([n_clean_fq_1,n_clean_fq_2]),sep.join([t_clean_fq_1,t_clean_fq_2]),n_done_log,t_done_log,U])


def FQ2Bam(normalFastq=None,tumorFastq=None,picard=None,qcbase=None,qcwarn1=None,test=None,mf=None,outdir='.',name='Sample'):
  '[A] Do FastqToSam for GATK4'
  if normalFastq==None and tumorFastq==None:
    print("Error: Both normal and tumor FASTQ files do not exist!", file=sys.stderr)
    sys.exit(1)
  if normalFastq!=None:
    n_pe_fq = sep.join(normalFastq)
  else:
    n_pe_fq = None
  if tumorFastq!=None:
    t_pe_fq = sep.join(tumorFastq)
  else:
    t_pe_fq = None
  ndir=os.path.join(outdir,'normal',dir1)
  tdir=os.path.join(outdir,'tumor',dir1)
  n_ubam = os.path.join(ndir,'unmap.bam')
  t_ubam = os.path.join(tdir,'unmap.bam')
  ntmpdir=os.path.join(ndir,'tmp')
  ttmpdir=os.path.join(tdir,'tmp')
  if normalFastq!=None:
    assert not os.makedirs(ntmpdir,exist_ok=True) ,'Cannot create %s' % ntmpdir
  if tumorFastq!=None:
    assert not os.makedirs(ttmpdir,exist_ok=True) ,'Cannot create %s' % ttmpdir

  n_log = os.path.join(ndir,'fq2sam.log')
  t_log = os.path.join(tdir,'fq2sam.log')
  n_done_log = os.path.join(ndir,'fq2sam_done.log')
  t_done_log = os.path.join(tdir,'fq2sam_done.log')

  n_done_cmd = '2> %s && echo "%s" > %s' % (n_log,name+'_n_fq2sam', n_done_log) # &> -> 2> (for ubuntu)
  t_done_cmd = '2> %s && echo "%s" > %s' % (t_log,name+'_t_fq2sam', t_done_log)

  mf.write('\n\n.MAKEFLOW CATEGORY FQ2Bam\n')

  if normalFastq!=None:
    n_cmd = sep.join(['java -jar',picard,'FastqToSam TMP_DIR='+ntmpdir,'F1='+normalFastq[0],'F2='+normalFastq[1],'O='+n_ubam,'SM='+name+'_normal RG=normal PL=illumina LB=normal_lib PU=PE',n_done_cmd])
    n_oi=sep.join([n_ubam,n_done_log,':',n_pe_fq])
    mf.write('\n%s\n\t%s\n' % (n_oi,n_cmd))
  if tumorFastq!=None:
    t_cmd = sep.join(['java -jar',picard,'FastqToSam TMP_DIR='+ttmpdir,'F1='+tumorFastq[0],'F2='+tumorFastq[1],'O='+t_ubam,'SM='+name+'_tumor RG=tumor PL=illumina LB=tumor_lib PU=PE',t_done_cmd])
    t_oi=sep.join([t_ubam,t_done_log,':',t_pe_fq])
    mf.write('\n%s\n\t%s\n' % (t_oi,t_cmd))
  # QC warn 1
  if tumorFastq!=None and normalFastq!=None and 'GS519' in name:
    qcdir=os.path.join(outdir,'qc')
    log=os.path.join(qcdir,'qcwarn1.err')
    done_log = os.path.join(qcdir,'qcwarn1_done.log')
    done_cmd = 'echo "%s" > %s' % (name+'_qcwarn1',done_log)
    cmd=sep.join([qcwarn1,name,normalFastq[0],normalFastq[1],tumorFastq[0],tumorFastq[1],qcbase,'2>',log,'&&',done_cmd])
    oi=sep.join([done_log,':',normalFastq[0],normalFastq[1],tumorFastq[0],tumorFastq[1]])
    if not test:
      mf.write('\n%s\n\t%s\n' % (oi,cmd))

  if tumorFastq!=None and normalFastq!=None:
    return([n_ubam,t_ubam,n_done_log,t_done_log])
  elif tumorFastq==None:
    return([n_ubam,None,n_done_log,None])
  else:
    return([None,t_ubam,None,t_done_log])


def bwa(normalFastq=None,tumorFastq=None,mark=None,picard=None,cores=None,bwa='bwa',ref=None,mf=None,outdir='.',name='sample'):
  '[B][C] Do bwa mem or aln without -R option (read group) for GATK4'
  if normalFastq==None and tumorFastq==None:
    print("Error: Both normal and tumor FASTQ files do not exist!", file=sys.stderr)
    sys.exit(1)
  if normalFastq!=None:
    n_pe_fq = sep.join(normalFastq)
  else:
    n_pe_fq = None
  if tumorFastq!=None:
    t_pe_fq = sep.join(tumorFastq)
  else:
    t_pe_fq = None
  if mark=='aln':
    return bwaAlnMap(n_pe_fq,t_pe_fq,['gatk4'],picard=picard,cores=cores,ref=ref,bwa=bwa,mf=mf,outdir=outdir,name=name)
  else:
    n_bam_out = os.path.join(outdir,'normal',dir2,'map.bam')
    t_bam_out = os.path.join(outdir,'tumor',dir2,'map.bam')

    n_err_log = n_bam_out + '.err'
    t_err_log = t_bam_out + '.err'
    n_done_log = os.path.join(outdir,'normal',dir2,'map_done.log')
    t_done_log = os.path.join(outdir,"tumor",dir2,'map_done.log')

    n_done_cmd = 'echo "%s" > %s' % (name+'_normal_mapping_done',n_done_log)
    t_done_cmd = 'echo "%s" > %s' % (name+'_tumor_mapping_done',t_done_log)

    mf.write('\n\n.MAKEFLOW CATEGORY BWA_MEM\n')
    if normalFastq!=None:
      n_mf_cmd=sep.join([bwa,'mem -t %d' % cores,'-M',ref,n_pe_fq,'2>',n_err_log,'| samtools view -bS >',n_bam_out,'&&',n_done_cmd])
      n_mf_oi=sep.join([n_done_log,n_bam_out,':',n_pe_fq])
      mf.write('\n%s\n\t%s\n' % (n_mf_oi,n_mf_cmd))
    if tumorFastq!=None:
      t_mf_cmd=sep.join([bwa,'mem -t %d' % cores,'-M',ref,t_pe_fq,'2>',t_err_log,'| samtools view -bS >',t_bam_out,'&&',t_done_cmd])  # -F 4
      t_mf_oi=sep.join([t_done_log,t_bam_out,':',t_pe_fq])
      mf.write('\n%s\n\t%s\n' % (t_mf_oi,t_mf_cmd))

    if normalFastq!=None and tumorFastq!=None:
      return([n_bam_out,t_bam_out,n_done_log,t_done_log])
    elif tumorFastq==None:
      return([n_bam_out,None,n_done_log,None])
    else:
      return([None,t_bam_out,None,t_done_log])


#def MergeBam(umapBams,mapBams,mark,picard,cores,ref=None,qcsite=None,mf=None,outdir='.',name='sample'):
def MergeBam(umapBams,mapBams,mark,picard,cores,ref=None,mf=None,outdir='.',name='sample'):
  '[D] Do MergeBamAlignment for GATK4'
  ndir=os.path.join(outdir,'normal',dir2)
  tdir=os.path.join(outdir,'tumor',dir2)
  #qcdir = os.path.join(outdir,'qc')
  #nqcdir=os.path.join(qcdir,'normal','sites')
  #tqcdir=os.path.join(qcdir,'tumor','sites')
  if umapBams[0]!=None:
    #assert not os.makedirs(nqcdir,exist_ok=True) ,'Cannot create %s' % nqcdir
    bn=basename(mapBams[0])
  if umapBams[1]!=None:
    #assert not os.makedirs(tqcdir,exist_ok=True) ,'Cannot create %s' % tqcdir
    bn=basename(mapBams[1])
  if bn=='map.bam':
    n_bam_out = os.path.join(ndir,'merge_sorted.bam')
    t_bam_out = os.path.join(tdir,'merge_sorted.bam')
    ntmpdir=os.path.join(ndir,'tmp')
    ttmpdir=os.path.join(tdir,'tmp')
    n_done_log = os.path.join(ndir,'merge_sort_done.log')
    t_done_log = os.path.join(tdir,'merge_sort_done.log')
  elif bn=='aln.bam':
    n_bam_out = os.path.join(ndir,'merge_sorted.aln.bam')
    t_bam_out = os.path.join(tdir,'merge_sorted.aln.bam')
    ntmpdir=os.path.join(ndir,'tmp.aln')
    ttmpdir=os.path.join(tdir,'tmp.aln')
    n_done_log = os.path.join(ndir,'merge_sort_done.aln.log')
    t_done_log = os.path.join(tdir,'merge_sort_done.aln.log')
  else:
    print("The basename of bam file '%s' is neither 'map.bam' nor 'aln.bam'." % bn, file=sys.stderr)
    sys.exit(1)
  if umapBams[0]!=None:
    assert not os.makedirs(ntmpdir,exist_ok=True) ,'Cannot create %s' % ntmpdir
  if umapBams[1]!=None:
    assert not os.makedirs(ttmpdir,exist_ok=True) ,'Cannot create %s' % ttmpdir

  n_idx_log = n_bam_out + '.idx.log'
  t_idx_log = t_bam_out + '.idx.log'
  n_err_log = n_bam_out + '.err'
  t_err_log = t_bam_out + '.err'
  #nsum=os.path.join(nqcdir,'summary.xls')
  #tsum=os.path.join(tqcdir,'summary.xls')

  #if bn=='map.bam':
  #  n_sqc_cmd = sep.join(['sh',qcsite,n_bam_out,nqcdir,nsum])
  #  t_sqc_cmd = sep.join(['sh',qcsite,t_bam_out,tqcdir,tsum])
  n_idx_cmd = 'samtools index -@ %s %s 2> %s' % (cores,n_bam_out,n_idx_log) # &> -> 2> (for ubuntu)
  t_idx_cmd = 'samtools index -@ %s %s 2> %s' % (cores,t_bam_out,t_idx_log)
  n_done_cmd = 'echo "%s" > %s' % (name+'_normal_merge_done',n_done_log)
  t_done_cmd = 'echo "%s" > %s' % (name+'_tumor_merge_done',t_done_log)
  if bn=='map.bam':
    if umapBams[0]!=None:
      #n_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ntmpdir,'UNMAPPED='+umapBams[0],'ALIGNED='+mapBams[0],'O='+n_bam_out,'R='+ref,'SO=coordinate 2>',n_err_log,'&&',n_idx_cmd,'&&',n_sqc_cmd,'&&',n_done_cmd]) # &> -> 2> (for ubuntu)
      n_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ntmpdir,'UNMAPPED='+umapBams[0],'ALIGNED='+mapBams[0],'O='+n_bam_out,'R='+ref,'SO=coordinate 2>',n_err_log,'&&',n_idx_cmd,'&&',n_done_cmd]) # &> -> 2> (for ubuntu)
    else:
      n_mf_cmd=None
    if umapBams[1]!=None:
      #t_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ttmpdir,'UNMAPPED='+umapBams[1],'ALIGNED='+mapBams[1],'O='+t_bam_out,'R='+ref,'SO=coordinate 2>',t_err_log,'&&',t_idx_cmd,'&&',t_sqc_cmd,'&&',t_done_cmd])
      t_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ttmpdir,'UNMAPPED='+umapBams[1],'ALIGNED='+mapBams[1],'O='+t_bam_out,'R='+ref,'SO=coordinate 2>',t_err_log,'&&',t_idx_cmd,'&&',t_done_cmd])
    else:
      t_mf_cmd=None
  elif bn=='aln.bam':
    if umapBams[0]!=None:
      n_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ntmpdir,'UNMAPPED='+umapBams[0],'ALIGNED='+mapBams[0],'O='+n_bam_out,'R='+ref,'SO=coordinate 2>',n_err_log,'&&',n_idx_cmd,'&&',n_done_cmd])
    if umapBams[1]!=None:
      t_mf_cmd=sep.join(['java -jar',picard,'MergeBamAlignment VALIDATION_STRINGENCY=LENIENT TMP_DIR='+ttmpdir,'UNMAPPED='+umapBams[1],'ALIGNED='+mapBams[1],'O='+t_bam_out,'R='+ref,'SO=coordinate 2>',t_err_log,'&&',t_idx_cmd,'&&',t_done_cmd])
  else:
    print("The basename of bam file '%s' is neither 'map.bam' nor 'aln.bam'." % bn, file=sys.stderr)
    sys.exit(1)

  mf.write('\n\n.MAKEFLOW CATEGORY MergeBam\n')
  if umapBams[0]!=None:
    n_mf_oi=sep.join([n_done_log,n_bam_out,':',umapBams[0],mapBams[0],mark[0],mark[2]])
    mf.write('\n%s\n\t%s\n' % (n_mf_oi,n_mf_cmd))
  if umapBams[1]!=None:
    t_mf_oi=sep.join([t_done_log,t_bam_out,':',umapBams[1],mapBams[1],mark[1],mark[3]])
    mf.write('\n%s\n\t%s\n' % (t_mf_oi,t_mf_cmd))

  if umapBams[0]!=None and umapBams[1]!=None:
    return([n_bam_out,t_bam_out,n_done_log,t_done_log])
  elif umapBams[1]==None:
    return([n_bam_out,None,n_done_log,None])
  else:
    return([None,t_bam_out,None,t_done_log])


def CleanP(nfq,tfq,fastp,mf=None,outdir='.',name='Sample'):
  '[f]'
  n_pe_fq = sep.join(nfq)
  n_clean_fq_1 = os.path.join(outdir,'normal',dir1,"cleanP_1.fq.gz")
  n_clean_fq_2 = os.path.join(outdir,'normal',dir1,"cleanP_2.fq.gz")
  n_clean_outs = sep.join([n_clean_fq_1,n_clean_fq_2])
  t_pe_fq = sep.join(tfq)
  t_clean_fq_1 = os.path.join(outdir,'tumor',dir1,"cleanP_1.fq.gz")
  t_clean_fq_2 = os.path.join(outdir,'tumor',dir1,"cleanP_2.fq.gz")
  t_clean_outs = sep.join([t_clean_fq_1,t_clean_fq_2])

  n_done_log = os.path.join(outdir,'normal',dir1,"success_cleanP.log")
  t_done_log = os.path.join(outdir,'tumor',dir1,"success_cleanP.log")

  n_done_cmd='&& echo "%s" > %s' % (name+'_n_cleanP', n_done_log)
  t_done_cmd='&& echo "%s" > %s' % (name+'_t_cleanP', t_done_log)
  n_cmd=sep.join([fastp,'-c -i',nfq[0],'-I',nfq[1],'-o',n_clean_fq_1,'-O',n_clean_fq_2,n_done_cmd])
  t_cmd=sep.join([fastp,'-c -i',tfq[0],'-I',tfq[1],'-o',t_clean_fq_1,'-O',t_clean_fq_2,t_done_cmd])

  n_oi=sep.join([n_clean_fq_1,n_clean_fq_2,n_done_log,':',n_pe_fq])
  t_oi=sep.join([t_clean_fq_1,t_clean_fq_2,t_done_log,':',t_pe_fq])
  mf.write("\n\n.MAKEFLOW CATEGORY FastqCleanP\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_cmd))
  return([n_clean_outs,t_clean_outs,n_done_log,t_done_log])


def FqQualControl(nFq,tFq,nMark=None,tMark=None,mf=None,outdir='.',name='sample'):
  '[2] (also called by [J]) nFq and tFq both have 2 FASTQ filenames joined by a space'
  qcdir = os.path.join(outdir,'qc')
  ndir=os.path.join(qcdir,'normal')
  tdir=os.path.join(qcdir,'tumor')
  if nFq!=None:
    assert not os.makedirs(ndir,exist_ok=True),'Cannot create %s' % ndir
  if tFq!=None:
    assert not os.makedirs(tdir,exist_ok=True),'Cannot create %s' % tdir

  n_run_log = os.path.join(ndir,'fastqc.log')
  t_run_log = os.path.join(tdir,'fastqc.log')
  n_done_log = os.path.join(ndir,'fastqc_done.log')
  t_done_log = os.path.join(tdir,'fastqc_done.log')

  n_done_cmd = 'echo "%s" > %s' % (name+'_n_fastqc',n_done_log)
  t_done_cmd = 'echo "%s" > %s' % (name+'_t_fastqc',t_done_log)
  mf.write("\n\n.MAKEFLOW CATEGORY FqQualControl\n")
  if nFq!=None:
    n_cmd=sep.join(['fastqc -t 2 -o',ndir,nFq,'>',n_run_log,'&&',n_done_cmd]) # &> -> > (for ubuntu)
    n_oi=sep.join([n_done_log,':',nFq])
    mf.write('\n%s\n\t%s\n' % (n_oi,n_cmd))
  if tFq!=None:
    t_cmd=sep.join(['fastqc -t 2 -o',tdir,tFq,'>',t_run_log,'&&',t_done_cmd])
    t_oi=sep.join([t_done_log,':',tFq])
    mf.write('\n%s\n\t%s\n' % (t_oi,t_cmd))
  if nFq!=None and tFq!=None:
    return([nFq,tFq,n_done_log,t_done_log])
  elif nFq==None:
    return([None,tFq,None,t_done_log])
  else:
    return([nFq,None,n_done_log,None])


def QualControl(nFqs,tFqs,nBam,tBam,nMark,tMark,qcprog=None,qcsite=None,gatk3=None,ref=None,intvl=None,qcwarn2=None,test=None,bed=None,mf=None,outdir='.',name='sample'):
  '[J] nFqs and tFqs both have a list of 2 FASTQ filenames'
  if nFqs==None and tFqs==None:
    print("Error: Both normal and tumor FASTQ files do not exist!", file=sys.stderr)
    sys.exit(1)
  if nFqs!=None:
    nFq = sep.join(nFqs)
  else:
    nFq = None
  if tFqs!=None:
    tFq = sep.join(tFqs)
  else:
    tFq = None
  qcdir = os.path.join(outdir,'qc')
  ndir=os.path.join(qcdir,'normal')
  tdir=os.path.join(qcdir,'tumor')
  if nFqs!=None:
    assert not os.makedirs(ndir,exist_ok=True),'Cannot create %s' % ndir
  if tFqs!=None:
    assert not os.makedirs(tdir,exist_ok=True),'Cannot create %s' % tdir

  fastqc_done_log = os.path.join(tdir,'fastqc_done.log')
  if not os.path.isfile(fastqc_done_log):  # only for 2nd run
    FqQualControl(nFq,tFq,mf=mf,outdir=outdir,name=name)

  n_run_log = os.path.join(ndir,'qc.log')
  t_run_log = os.path.join(tdir,'qc.log')
  n_done_log = os.path.join(ndir,'qc_done.log')
  t_done_log = os.path.join(tdir,'qc_done.log')
  nout = os.path.join(ndir,name+'.QC.12indexes.xls')
  tout = os.path.join(tdir,name+'.QC.12indexes.xls')

  n_done_cmd = 'echo "%s" > %s' % (name+'_n_qc',n_done_log)
  t_done_cmd = 'echo "%s" > %s' % (name+'_t_qc',t_done_log)
  mf.write("\n\n.MAKEFLOW CATEGORY QualControl\n")
  if nFqs!=None:
    n_cmd=sep.join(['cd',ndir,'&& perl',qcprog,nFq,nBam,bed,name,'>',n_run_log,'&&',n_done_cmd]) # &> -> > (for ubuntu)
    n_oi=sep.join([nout,n_done_log,':',nFq,nBam])
    mf.write('\n%s\n\t%s\n' % (n_oi,n_cmd))
  if tFqs!=None:
    t_cmd=sep.join(['cd',tdir,'&& perl',qcprog,tFq,tBam,bed,name,'>',t_run_log,'&&',t_done_cmd])
    t_oi=sep.join([tout,t_done_log,':',tFq,tBam])
    mf.write('\n%s\n\t%s\n' % (t_oi,t_cmd))

  if nFqs!=None and tFqs!=None:
    qcdir=os.path.join(outdir,'qc','sites')
    assert not os.makedirs(qcdir,exist_ok=True) ,'Cannot create %s' % qcdir
    sum=os.path.join(qcdir,'summary.xls')
    log=os.path.join(qcdir,'qcsites.err')
    sqc_cmd=sep.join(['python3',qcsite,nBam,tBam,qcdir,sum,gatk3,ref,intvl,'2>',log])
    #sqc_cmd=sep.join(['sh',qcsite,nBam,tBam,qcdir,sum,'2>',log])
    if 'GS519' in name:
      log=os.path.join(qcdir,'qcwarn2.err')
      qcwarn2_cmd=sep.join(['ruby',qcwarn2,name,nout,tout,sum,'2>',log])
    oi=sep.join([sum,':',nBam,tBam,nout,tout])
    if 'GS519' in name:
      if test:
        mf.write('\n%s\n\t%s\n' % (oi,sep.join([sqc_cmd])))
      else:
        mf.write('\n%s\n\t%s\n' % (oi,sep.join([sqc_cmd,'&&',qcwarn2_cmd])))
    else:
      mf.write('\n%s\n\t%s\n' % (oi,sep.join([sqc_cmd])))

  if nFqs!=None and tFqs!=None:
    return([nout,tout,n_done_log,t_done_log])
  elif nFqs==None:
    return([None,tout,None,t_done_log])
  else:
    return([nout,None,n_done_log,None])


def QCrun(nrfq,trfq,mf=None,outdir='.',name='Sample'):
  '[q] Make makeflow script for long-time QC runs'
  wdir=os.path.join(outdir,'qc')
  mf.write('\n\n.MAKEFLOW CATEGORY QCrun\n')

  n=-1
  fqs=[nrfq[0],nrfq[1],trfq[0],trfq[1]]
  tis=['n1','n2','t1','t2']
  for i in range(0,4):
    stat=fqs[i]
    stat=stat.replace('.gz','.fqStat.txt')
    out=os.path.join(wdir,'qcout_%s.txt' % tis[i])
    done_log=os.path.join(wdir,'qcdone_%s.log' % tis[i])
    if not os.path.isfile(stat):
      n=n+1
      done_cmd='echo \"%s\" > %s' % ('%s get base and read number done' % fqs[i],done_log)
      cmd='zcat %s | awk \'NR%%4==2{n++;l+=length($0)}END{print "ReadNum\\t"n"\\nBaseNum\\t"l}\' > %s && %s' % (fqs[i],out,done_cmd)
      oi=sep.join([out,done_log,':',fqs[i]])
      mf.write('\n') if n else ''
      mf.write('%s\n\t%s\n' % (oi,cmd))
  return([out,done_log])


def QCstat(nrfq,trfq,mf=None,outdir='.',name='Sample'):
  'Do short-time QC runs and stats'
  wdir=os.path.join(outdir,'qc')
  qcout=os.path.join(wdir,'qcout.txt')
  done_log=os.path.join(wdir,'qc_done.log')

  nr=nb=0
  fqs=[nrfq[0],nrfq[1],trfq[0],trfq[1]]
  tis=['n1','n2','t1','t2']
  for i in range(0,4):
    stat=fqs[i]
    stat=stat.replace('.gz','.fqStat.txt')
    if not os.path.isfile(stat):
      stat=os.path.join(wdir,'qcout_%s.txt' % tis[i])
      if not os.path.isfile(stat):
        print("Stat file for '%s' NOT exist." % fqs[i], file=sys.stderr)
        sys.exit(1)
    nr=nr+int(subprocess.run('grep ReadNum %s | cut -f2' % stat,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.strip())
    nb=nb+int(subprocess.run('grep BaseNum %s | cut -f2' % stat,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.strip())


def bwaMapping(nCleanFq,tCleanFq,nMark,tMark,cores,bwa='bwa',ref=None,mf=None,outdir='.',name='sample',PL='illumina',LB='library'):
  '[3]'
  n_RG = r'"@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\\tPU:PE"' % ('normal',name+'_normal',PL,LB)
  t_RG = r'"@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\\tPU:PE"' % ("tumor",name+'_tumor',PL,LB)

  n_bam_out = os.path.join(outdir,'normal',dir2,'map.bam')
  t_bam_out = os.path.join(outdir,'tumor',dir2,'map.bam')

  n_err_log = n_bam_out + '.err'
  t_err_log = t_bam_out + '.err'
  n_done_log = os.path.join(outdir,'normal',dir2,"map_done.log")
  t_done_log = os.path.join(outdir,"tumor",dir2,"map_done.log")

  n_done_cmd = "echo \"%s\" > %s" % (name+"_normal_mapping_done",n_done_log)
  t_done_cmd = "echo \"%s\" > %s" % (name+"_tumor_mapping_done",t_done_log)
  n_mf_cmd = "%s mem -t %d -M -R %s %s %s 2> %s | samtools view -bS > %s && %s" % (bwa,cores,n_RG,ref,nCleanFq,n_err_log,n_bam_out,n_done_cmd)
  t_mf_cmd = "%s mem -t %d -M -R %s %s %s 2> %s | samtools view -bS > %s && %s" % (bwa,cores,t_RG,ref,tCleanFq,t_err_log,t_bam_out,t_done_cmd)  # -F 4

  n_mf_oi=sep.join([n_done_log,n_bam_out,':',nCleanFq,nMark])
  t_mf_oi=sep.join([t_done_log,t_bam_out,':',tCleanFq,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY BWA-MEM\n")
  mf.write("%s\n\t%s\n\n" % (n_mf_oi,n_mf_cmd))
  mf.write("%s\n\t%s\n" % (t_mf_oi,t_mf_cmd))
  return([n_bam_out,t_bam_out,n_done_log,t_done_log])


def bwaAlnMap(nFq,tFq,mark,picard,cores,bwa='bwa',ref=None,mf=None,outdir='.',name='sample',PL='ILLUMINA',LB='library'):
  '[b], also called by [C]'
  if nFq==None and tFq==None:
    print("Error: Both normal and tumor FASTQ files do not exist!", file=sys.stderr)
    sys.exit(1)
  n_sai1 = os.path.join(outdir,'normal',dir2,'1.sai')
  t_sai1 = os.path.join(outdir,'tumor',dir2,'1.sai')
  n_sai2 = os.path.join(outdir,'normal',dir2,'2.sai')
  t_sai2 = os.path.join(outdir,'tumor',dir2,'2.sai')
  n_bam_tmp = os.path.join(outdir,'normal',dir2,'temp.bam')
  t_bam_tmp = os.path.join(outdir,'tumor',dir2,'temp.bam')
  n_bam_out = os.path.join(outdir,'normal',dir2,'aln.bam')
  t_bam_out = os.path.join(outdir,'tumor',dir2,'aln.bam')

  n_err_log1 = n_bam_out + '1.err'
  t_err_log1 = t_bam_out + '1.err'
  n_err_log2 = n_bam_out + '2.err'
  t_err_log2 = t_bam_out + '2.err'
  n_err_log = n_bam_out + '.err'
  t_err_log = t_bam_out + '.err'
  n_adrg_log = os.path.join(outdir,'normal',dir2,'picardAddRG.log')
  t_adrg_log = os.path.join(outdir,'tumor',dir2,'picardAddRG.log')
  n_done_log = os.path.join(outdir,'normal',dir2,'aln_done.log')
  t_done_log = os.path.join(outdir,'tumor',dir2,'aln_done.log')

  if nFq!=None:
    nCFq = nFq.split()
    n_done_cmd = "echo \"%s\" > %s" % (name+"_normal_aln_mapping_done",n_done_log)
    n_mf_cmd1 = "%s aln -t %d %s %s > %s 2> %s" % (bwa,cores,ref,nCFq[0],n_sai1,n_err_log1)
    n_mf_cmd2 = "%s aln -t %d %s %s > %s 2> %s" % (bwa,cores,ref,nCFq[1],n_sai2,n_err_log2)
    if mark.__len__()==2:
      n_mf_cmd3 = "%s sampe %s %s %s %s 2> %s | samtools view -bS > %s && rm %s %s" % (bwa,ref,n_sai1,n_sai2,nFq,n_err_log,n_bam_tmp,n_sai1,n_sai2)
      n_mf_cmd4 = "java -jar %s AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=%s O=%s RGID=normal RGLB=%s RGPL=%s RGPU=PE RGSM=%s 2> %s && %s" % (picard,n_bam_tmp,n_bam_out,LB,PL,name+'_normal',n_adrg_log,n_done_cmd) # &> -> 2> (for ubuntu)
    else:  # gatk4
      n_mf_cmd3 = "%s sampe %s %s %s %s 2> %s | samtools view -bS > %s && rm %s %s && %s" % (bwa,ref,n_sai1,n_sai2,nFq,n_err_log,n_bam_out,n_sai1,n_sai2,n_done_cmd)

  if tFq!=None:
    tCFq = tFq.split()
    t_done_cmd = "echo \"%s\" > %s" % (name+"_tumor_aln_mapping_done",t_done_log)
    t_mf_cmd1 = "%s aln -t %d %s %s > %s 2> %s" % (bwa,cores,ref,tCFq[0],t_sai1,t_err_log1)
    t_mf_cmd2 = "%s aln -t %d %s %s > %s 2> %s" % (bwa,cores,ref,tCFq[1],t_sai2,t_err_log2)
    if mark.__len__()==2:
      t_mf_cmd3 = "%s sampe %s %s %s %s 2> %s | samtools view -bS > %s && rm %s %s" % (bwa,ref,t_sai1,t_sai2,tFq,t_err_log,t_bam_tmp,t_sai1,t_sai2)  # -F 4
      t_mf_cmd4 = "java -jar %s AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=%s O=%s RGID=tumor RGLB=%s RGPL=%s RGPU=PE RGSM=%s 2> %s && %s" % (picard,t_bam_tmp,t_bam_out,LB,PL,name+'_tumor',t_adrg_log,t_done_cmd)
    else:  # gatk4
      t_mf_cmd3 = "%s sampe %s %s %s %s 2> %s | samtools view -bS > %s && rm %s %s && %s" % (bwa,ref,t_sai1,t_sai2,tFq,t_err_log,t_bam_out,t_sai1,t_sai2,t_done_cmd)

  mf.write("\n\n.MAKEFLOW CATEGORY BWA_ALN\n")
  if mark.__len__()==2:  # mark: n_done_log, t_done_log
    if nFq!=None:
      n_mf_oi=sep.join([n_sai1,':',nCFq[0],mark[0]])
    if tFq!=None:
      t_mf_oi=sep.join([t_sai1,':',tCFq[0],mark[1]])
  else:
    if nFq!=None:
      n_mf_oi=sep.join([n_sai1,':',nCFq[0]])
    if tFq!=None:
      t_mf_oi=sep.join([t_sai1,':',tCFq[0]])
  if nFq!=None:
    mf.write("%s\n\t%s\n\n" % (n_mf_oi,n_mf_cmd1))
    n_mf_oi=sep.join([n_sai2,':',nCFq[1]])
    mf.write("%s\n\t%s\n\n" % (n_mf_oi,n_mf_cmd2))
  if tFq!=None:
    mf.write("%s\n\t%s\n\n" % (t_mf_oi,t_mf_cmd1))
    t_mf_oi=sep.join([t_sai2,':',tCFq[1]])
    mf.write('%s\n\t%s\n\n' % (t_mf_oi,t_mf_cmd2))
  if mark.__len__()==2:
    if nFq!=None:
      n_mf_oi=sep.join([n_bam_tmp,':',nFq,n_sai1,n_sai2])
      mf.write('%s\n\t%s\n\n' % (n_mf_oi,n_mf_cmd3))
      n_mf_oi=sep.join([n_done_log,n_bam_out,':',n_bam_tmp])
      mf.write('%s\n\t%s\n\n' % (n_mf_oi,n_mf_cmd4))
    if tFq!=None:
      t_mf_oi=sep.join([t_bam_tmp,':',tFq,t_sai1,t_sai2])
      mf.write('%s\n\t%s\n\n' % (t_mf_oi,t_mf_cmd3))
      t_mf_oi=sep.join([t_done_log,t_bam_out,':',t_bam_tmp])
      mf.write('%s\n\t%s\n\n' % (t_mf_oi,t_mf_cmd4))
  else:
    if nFq!=None:
      n_mf_oi=sep.join([n_done_log,n_bam_out,':',nFq,n_sai1,n_sai2])
      mf.write("%s\n\t%s\n\n" % (n_mf_oi,n_mf_cmd3))
    if tFq!=None:
      t_mf_oi=sep.join([t_done_log,t_bam_out,':',tFq,t_sai1,t_sai2])
      mf.write("%s\n\t%s\n\n" % (t_mf_oi,t_mf_cmd3))

  if nFq!=None and tFq!=None:
    return([n_bam_out,t_bam_out,n_done_log,t_done_log])
  elif nFq==None:
    return([None,t_bam_out,None,t_done_log])
  else:
    return([n_bam_out,None,n_done_log,None])


def Sort(n_bam,t_bam,nMark,tMark,samtools="samtools",mf=None,outdir='.',name='Sample'):
  '[4] Sort bam and output bam file, then samtools index bam 2> index.err.log'
  ndir=os.path.join(outdir,'normal',dir2)
  tdir=os.path.join(outdir,'tumor',dir2)
  bn=basename(n_bam)
  if bn=='map.bam':
    n_bam_out = os.path.join(ndir,'sorted.bam')
    t_bam_out = os.path.join(tdir,'sorted.bam')
    n_done_log = os.path.join(ndir,'Sort_done.log')
    t_done_log = os.path.join(tdir,'Sort_done.log')
  elif bn=='aln.bam':
    n_bam_out = os.path.join(ndir,'sorted.aln.bam')
    t_bam_out = os.path.join(tdir,'sorted.aln.bam')
    n_done_log = os.path.join(ndir,'Sort_done.aln.log')
    t_done_log = os.path.join(tdir,'Sort_done.aln.log')
  else:
    print("The basename of bam file '%s' is neither 'map.bam' nor 'aln.bam'." % bn, file=sys.stderr)
    sys.exit(1)

  n_sort_log = n_bam_out + '.sort.log'
  t_sort_log = t_bam_out + '.sort.log'
  n_idx_log = n_bam_out + '.idx.log'
  t_idx_log = t_bam_out + '.idx.log'

  n_done_cmd = 'echo "%s" > %s' % (name+'_n_sortIndex_done',n_done_log)
  t_done_cmd = 'echo "%s" > %s' % (name+'_t_sortIndex_done',t_done_log)
  n_sort_cmd = '%s sort --threads 5 %s -o %s 2> %s' % (samtools,n_bam,n_bam_out,n_sort_log) # &> -> 2> (for ubuntu)
  t_sort_cmd = '%s sort --threads 5 %s -o %s 2> %s' % (samtools,t_bam,t_bam_out,t_sort_log)
  n_idx_cmd = '%s index %s 2> %s' % (samtools,n_bam_out,n_idx_log)
  t_idx_cmd = '%s index %s 2> %s' % (samtools,t_bam_out,t_idx_log)
  n_cmd = '%s && %s && %s' % (n_sort_cmd,n_idx_cmd,n_done_cmd)
  t_cmd = '%s && %s && %s' % (t_sort_cmd,t_idx_cmd,t_done_cmd)

  n_oi=sep.join([n_done_log,n_bam_out,':',n_bam,nMark])
  t_oi=sep.join([t_done_log,t_bam_out,':',t_bam,tMark])
  mf.write('\n\n.MAKEFLOW CATEGORY bam2sortBamAndIndex\n')
  mf.write('%s\n\t%s\n\n' % (n_oi,n_cmd))
  mf.write('%s\n\t%s\n' % (t_oi,t_cmd))
  return([n_bam_out,t_bam_out,n_done_log,t_done_log])


def withDupBamStats():
  " bam stats with duplicates"
  #TODO


def MarkDuplicates(n_bam=None,nMark=None,t_bam=None,tMark=None,picard=None,mf=None,outdir='.',name='Sample'):
  '[5] Also for GATK4'
  ndir=os.path.join(outdir,'normal',dir3)
  tdir=os.path.join(outdir,'tumor',dir3)
  if n_bam!=None:
    bn=basename(n_bam)
  if t_bam!=None:
    bn=basename(t_bam)
  if 'aln' not in bn:
    n_bam_out = os.path.join(ndir,'sortedDedup.bam')
    t_bam_out = os.path.join(tdir,'sortedDedup.bam')
    n_dedup_metrics = os.path.join(ndir,'dedup_metrics.txt')
    t_dedup_metrics = os.path.join(tdir,'dedup_metrics.txt')
    ntmpdir=os.path.join(ndir,'tmp')
    ttmpdir=os.path.join(tdir,'tmp')
    n_dedup_log = os.path.join(ndir,'dedup.log')
    t_dedup_log = os.path.join(tdir,'dedup.log')
    n_done_log = os.path.join(ndir,'dedup_done.log')
    t_done_log = os.path.join(tdir,'dedup_done.log')
  else:
    n_bam_out = os.path.join(ndir,'sortedDedup.aln.bam')
    t_bam_out = os.path.join(tdir,'sortedDedup.aln.bam')
    n_dedup_metrics = os.path.join(ndir,"dedup_metrics.aln.txt")
    t_dedup_metrics = os.path.join(tdir,"dedup_metrics.aln.txt")
    ntmpdir=os.path.join(ndir,'tmp.aln')
    ttmpdir=os.path.join(tdir,'tmp.aln')
    n_dedup_log = os.path.join(ndir,'dedup.aln.log')
    t_dedup_log = os.path.join(tdir,'dedup.aln.log')
    n_done_log = os.path.join(ndir,'dedup_done.aln.log')
    t_done_log = os.path.join(tdir,'dedup_done.aln.log')

  n_done_cmd = 'echo "%s" > %s' % (name+'_n_dedup_done',n_done_log)
  t_done_cmd = 'echo "%s" > %s' % (name+'_t_dedup_done',t_done_log)

  mf.write('\n\n.MAKEFLOW CATEGORY markDuplicates\n')
  if n_bam!=None:
    n_cmd=sep.join(['java -jar',picard,'MarkDuplicates I='+n_bam,'O='+n_bam_out,'M='+n_dedup_metrics,'ASO=coordinate CREATE_INDEX=true TMP_DIR='+ntmpdir,'2>',n_dedup_log,'&&',n_done_cmd]) # &> -> 2> (for ubuntu)
    n_oi=sep.join([n_done_log,n_bam_out,':',n_bam,nMark])
    mf.write('\n%s\n\t%s\n' % (n_oi,n_cmd))
  if t_bam!=None:
    t_cmd=sep.join(['java -jar',picard,'MarkDuplicates I='+t_bam,'O='+t_bam_out,'M='+t_dedup_metrics,'ASO=coordinate CREATE_INDEX=true TMP_DIR='+ttmpdir,'2>',t_dedup_log,'&&',t_done_cmd])
    t_oi=sep.join([t_done_log,t_bam_out,':',t_bam,tMark])
    mf.write('\n%s\n\t%s\n' % (t_oi,t_cmd))
  if n_bam!=None and t_bam!=None:
    return([n_bam_out,t_bam_out,n_done_log,t_done_log])
  elif n_bam==None:
    return([None,t_bam_out,None,t_done_log])
  else:
    return([n_bam_out,None,n_done_log,None])


def Dedup(n_bam,t_bam,nMark,tMark,dedup=None,mf=None,outdir='.',name='Sample'):
  '[g] OpenGene pipeline'
  n_bam_out = os.path.join(outdir,'normal',dir3,"sort.dedup.bam")
  t_bam_out = os.path.join(outdir,'tumor',dir3,"sort.dedup.bam")
  n_done_log = os.path.join(outdir,'normal',dir3,"dedupOp_done.log")
  t_done_log = os.path.join(outdir,"tumor",dir3,"dedupOp_done.log")

  n_done_cmd = "echo \"%s\" > %s" % (name+"_n_dedup_done",n_done_log)
  t_done_cmd = "echo \"%s\" > %s" % (name+"_t_dedup_done",t_done_log)
  n_idx_cmd="samtools index %s" % n_bam_out
  t_idx_cmd="samtools index %s" % t_bam_out
  n_cmd = "python %s -1 %s -o %s && %s && %s" % (dedup,n_bam,n_bam_out,n_idx_cmd,n_done_cmd)
  t_cmd = "python %s -1 %s -o %s && %s && %s" % (dedup,t_bam,t_bam_out,t_idx_cmd,t_done_cmd)

  n_oi=sep.join([n_done_log,n_bam_out,':',n_bam,nMark])
  t_oi=sep.join([t_done_log,t_bam_out,':',t_bam,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY DedupOpengene\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_cmd))
  return([n_bam_out,t_bam_out,n_done_log,t_done_log])


def mpileup(n_bam,t_bam,nMark,tMark,samtools=None,ref=None,bed=None,mf=None,outdir='.',name='Sample'):
  '[h] OpenGene pipeline'
  n_out = os.path.join(outdir,'normal',dir3,"sort.dedup.mpileup")
  t_out = os.path.join(outdir,'tumor',dir3,"sort.dedup.mpileup")
  n_done_log = os.path.join(outdir,'normal',dir3,"mpileup_done.log")
  t_done_log = os.path.join(outdir,"tumor",dir3,"mpileup_done.log")
  n_err = os.path.join(outdir,'normal',dir3,"mpileup.err")
  t_err = os.path.join(outdir,"tumor",dir3,"mpileup.err")

  n_done_cmd = "echo \"%s\" > %s" % (name+"_n_dedup_done",n_done_log)
  t_done_cmd = "echo \"%s\" > %s" % (name+"_t_dedup_done",t_done_log)
  n_cmd = "%s mpileup -B -Q 20 -C 50 -q 20 -d 20000 -f %s -l %s %s > %s 2> %s && %s" % (samtools,ref,bed,n_bam,n_out,n_err,n_done_cmd)
  t_cmd = "%s mpileup -B -Q 20 -C 50 -q 20 -d 20000 -f %s -l %s %s > %s 2> %s && %s" % (samtools,ref,bed,t_bam,t_out,t_err,t_done_cmd)

  n_oi=sep.join([n_done_log,n_out,':',n_bam,nMark])
  t_oi=sep.join([t_done_log,t_out,':',t_bam,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY mpileup\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_cmd))
  return([n_out,t_out,n_done_log,t_done_log])


def varscan(n_mpileup,t_mpileup,nMark,tMark,varscan=None,mf=None,outdir='.',name='Sample'):
  '[i] OpenGene pipeline'
  wdir=os.path.join(outdir,'variants','somatic',dirs)
  n_snp_out = os.path.join(wdir,"vs.snp.normal.vcf")
  t_snp_out = os.path.join(wdir,"vs.snp.tumor.vcf")
  n_indel_out = os.path.join(wdir,"vs.indel.normal.vcf")
  t_indel_out = os.path.join(wdir,"vs.indel.tumor.vcf")
  n_snp_done_log = os.path.join(wdir,"varscan_done.snp.normal.log")
  t_snp_done_log = os.path.join(wdir,"varscan_done.snp.tumor.log")
  n_indel_done_log = os.path.join(wdir,"varscan_done.indel.normal.log")
  t_indel_done_log = os.path.join(wdir,"varscan_done.indel.tumor.log")
  n_snp_err = os.path.join(wdir,"varscan.normal.snp.err")
  t_snp_err = os.path.join(wdir,"varscan.tumor.snp.err")
  n_indel_err = os.path.join(wdir,"varscan.normal.indel.err")
  t_indel_err = os.path.join(wdir,"varscan.tumor.indel.err")

  n_done_cmd = "echo \"%s\" > %s" % (name+"_n_snp_varscan_done",n_snp_done_log)
  t_done_cmd = "echo \"%s\" > %s" % (name+"_t_snp_varscan_done",t_snp_done_log)
  n_snp_cmd = "java -jar %s mpileup2snp %s --min-coverage 4 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.001 --min-freq-for-hom 0.90 --output-vcf 1 > %s 2> %s" % (varscan,n_mpileup,n_snp_out,n_snp_err)
  t_snp_cmd = "java -jar %s mpileup2snp %s --min-coverage 4 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.001 --min-freq-for-hom 0.90 --output-vcf 1 > %s 2> %s" % (varscan,t_mpileup,t_snp_out,t_snp_err)
  n_cmd = "%s && %s" % (n_snp_cmd,n_done_cmd)
  t_cmd = "%s && %s" % (t_snp_cmd,t_done_cmd)

  n_oi=sep.join([n_snp_done_log,n_snp_out,':',n_mpileup,nMark])
  t_oi=sep.join([t_snp_done_log,t_snp_out,':',t_mpileup,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY varscan\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n\n" % (t_oi,t_cmd))

  n_done_cmd = "echo \"%s\" > %s" % (name+"_n_indel_varscan_done",n_indel_done_log)
  t_done_cmd = "echo \"%s\" > %s" % (name+"_t_indel_varscan_done",t_indel_done_log)
  n_indel_cmd = "java -jar %s mpileup2indel %s --min-coverage 4 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.001 --min-freq-for-hom 0.90 --output-vcf 1 > %s 2> %s" % (varscan,n_mpileup,n_indel_out,n_indel_err)
  t_indel_cmd = "java -jar %s mpileup2indel %s --min-coverage 4 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.001 --min-freq-for-hom 0.90 --output-vcf 1 > %s 2> %s" % (varscan,t_mpileup,t_indel_out,t_indel_err)
  n_cmd = "%s && %s" % (n_indel_cmd,n_done_cmd)
  t_cmd = "%s && %s" % (t_indel_cmd,t_done_cmd)

  n_oi=sep.join([n_indel_done_log,n_indel_out,':',n_mpileup,nMark])
  t_oi=sep.join([t_indel_done_log,t_indel_out,':',t_mpileup,tMark])
  mf.write("%s\n\t%s\n\n" % (n_oi,n_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_cmd))
  return([n_snp_out,t_snp_out,n_indel_out,t_indel_out,n_snp_done_log,t_snp_done_log,n_indel_done_log,t_indel_done_log])


def Mrbam(bam,bamLog,txt,txtLog,title,mrbamdir=None,filter=None,mf=None,outdir='.',name='Sample'):
  '[k] OpenGene pipeline'
  wdir=os.path.join(outdir,'variants','somatic',dira)
  out=os.path.join(wdir,'mbout.%s.txt' % title)
  outflt=os.path.join(wdir,'mbout.flt.%s.txt' % title)
  done_log=os.path.join(wdir,'mrbam_done.%s.log' % title)

  done_cmd='echo \"%s\" > %s' % (name+'_mrbam_done',done_log)
  cmd='cd %s && python3 -m MrBam.main -o %s -m 3 -q 25 --fast --cfdna %s --skip 1 %s && perl %s %s %s 2 0.01 && cd %s && %s' % (mrbamdir,out,bam,txt,filter,out,outflt,outdir,done_cmd)

  oi=sep.join([done_log,outflt,':',bam,bamLog,txt,txtLog])
  mf.write('\n\n.MAKEFLOW CATEGORY MrBam%s\n' % title)
  mf.write('%s\n\t%s\n' % (oi,cmd))
  return([outflt,done_log])


def Mutscan(fq,title,mutscan=None,mf=None,outdir='.',name='Sample'):
  '[l] OpenGene pipeline'
  wdir=os.path.join(outdir,'variants','somatic',dirs)
  out=os.path.join(wdir,'msout.%s.txt' % title)
  outhtm=os.path.join(wdir,'msout.%s.html' % title)
  done_log=os.path.join(wdir,'mutscan_done.%s.log' % title)

  done_cmd='echo \"%s\" > %s' % (name+'_mutscan_done',done_log)
  cmd='%s -1 %s -2 %s -t 4 -h %s -S 2 > %s && %s' % (mutscan,fq[0],fq[1],outhtm,out,done_cmd)

  oi=sep.join([done_log,out,outhtm,':',fq[0],fq[1]])
  mf.write('\n\n.MAKEFLOW CATEGORY Mutscan%s\n' % title)
  mf.write('%s\n\t%s\n' % (oi,cmd))
  return([out,done_log])


def Genefuse(fq,fa,csv,title,genefuse=None,mf=None,outdir='.',name='Sample'):
  '[m] OpenGene pipeline'
  wdir=os.path.join(outdir,'variants','somatic',dirf)
  out=os.path.join(wdir,'gfout.%s.txt' % title)
  outhtm=os.path.join(wdir,'gfout.%s.html' % title)
  done_log=os.path.join(wdir,'genefuse_done.%s.log' % title)

  done_cmd='echo \"%s\" > %s' % (name+'_genefuse_done',done_log)
  cmd='%s -r %s -f %s -1 %s -2 %s -u 2 -d 50 -t 4 -h %s > %s && %s' % (genefuse,fa,csv,fq[0],fq[1],outhtm,out,done_cmd)

  oi=sep.join([done_log,out,outhtm,':',fq[0],fq[1]])
  mf.write('\n\n.MAKEFLOW CATEGORY Genefuse%s\n' % title)
  mf.write('%s\n\t%s\n' % (oi,cmd))
  return([out,done_log])


def withoutDupBamStats():
  'bam stats with removing duplicates'
  #TODO

def contEst(n_bam,t_bam,nMark,tMark,gatk=None,hapmap=None,ref=None,bedFlank=None,outdir='.',mf=None,name='Sample'):
  '[a]'
  #TODO
  baseReport_tsv = os.path.join(outdir,'qc',"contEst_baseReport.tsv")
  contEstResult_tsv = os.path.join(outdir,'qc',"contEst_result.tsv")

  log_mark_success = os.path.join(outdir,'qc',"contEst_done.log")
  mark_success = "&& echo \"%s\" > %s" % (name+"_contEst_done",log_mark_success)

  contEst_mf_cmd = "java -jar %s -T ContEst -R %s -I:eval %s -I:genotype %s --popfile %s -L %s --population ALL -br %s -o %s %s" % (
    gatk,ref,t_bam,n_bam,hapmap,bedFlank,baseReport_tsv,contEstResult_tsv,mark_success)

  mf_input_cmd=sep.join([log_mark_success,baseReport_tsv,contEstResult_tsv,':',nMark,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY tumorNormalContEst\n")
  mf.write("%s\n\t%s\n" % (mf_input_cmd,contEst_mf_cmd))
  return([contEstResult_tsv])

def IndelRealigner(n_bam,t_bam,nMark,tMark,gatk="",millsAndIndels="",dbSNP="",bedFlank="",ref="",outdir='.',name='Sample',mf=None):
  '[6]'
  n_bam_out = os.path.join(outdir,'normal',dir4,"realign.bam")
  t_bam_out = os.path.join(outdir,'tumor',dir4,"realign.bam")
  n_idx_out = os.path.join(outdir,'normal',dir4,"realign.bai")
  t_idx_out = os.path.join(outdir,'tumor',dir4,"realign.bai")
  n_intval_out = os.path.join(outdir,'normal',dir4,"realign.intervals")
  t_intval_out = os.path.join(outdir,'tumor',dir4,"realign.intervals")

  n_intval_log = os.path.join(outdir,'normal',dir4,'interval.log')
  t_intval_log = os.path.join(outdir,'tumor',dir4,'interval.log')
  n_intval_done_log = os.path.join(outdir,'normal',dir4,'interval_done.log')
  t_intval_done_log = os.path.join(outdir,'tumor',dir4,'interval_done.log')

  n_intval_done_cmd = "echo \"%s\" > %s" % (name+"_n_interval_done",n_intval_done_log)
  t_intval_done_cmd = "echo \"%s\" > %s" % (name+"_t_interval_done",t_intval_done_log)
  n_intval_cmd = "java -jar %s -T RealignerTargetCreator -I %s -o %s -known %s -known %s -L %s -R %s -nt 12 -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment 2> %s && %s" % (gatk,n_bam,n_intval_out,millsAndIndels,dbSNP,bedFlank,ref,n_intval_log,n_intval_done_cmd) # &> -> 2> (for ubuntu)
  t_intval_cmd = "java -jar %s -T RealignerTargetCreator -I %s -o %s -known %s -known %s -L %s -R %s -nt 12 -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment 2> %s && %s" % (gatk,t_bam,t_intval_out,millsAndIndels,dbSNP,bedFlank,ref,t_intval_log,t_intval_done_cmd)

  n_oi = sep.join([n_intval_done_log,':',n_bam,nMark])
  t_oi = sep.join([t_intval_done_log,':',t_bam,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY IndelRealigner\n")
  mf.write("%s\n\t%s\n\n" % (n_oi,n_intval_cmd))
  mf.write("%s\n\t%s\n\n" % (t_oi,t_intval_cmd))

  n_realn_log = os.path.join(outdir,'normal',dir4,"realign.log")
  t_realn_log = os.path.join(outdir,"tumor",dir4,"realign.log")
  n_realn_done_log = os.path.join(outdir,'normal',dir4,"realign_done.log")
  t_realn_done_log = os.path.join(outdir,"tumor",dir4,"realign_done.log")

  n_realn_done_cmd = "echo \"%s\" > %s" % (name+"_n_realign_done",n_realn_done_log)
  t_realn_done_cmd = "echo \"%s\" > %s" % (name+"_t_realign_done",t_realn_done_log)
  n_index_cmd = "samtools index %s %s" % (n_bam_out,n_idx_out)
  t_index_cmd = "samtools index %s %s" % (t_bam_out,t_idx_out)
  n_realn_cmd = "java -jar %s -T IndelRealigner -I %s -o %s -targetIntervals %s -known %s -known %s -R %s -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment 2> %s && %s && %s" % (gatk,n_bam,n_bam_out,n_intval_out,millsAndIndels,dbSNP,ref,n_realn_log,n_index_cmd,n_realn_done_cmd) # &> -> 2> (for ubuntu)
  t_realn_cmd = "java -jar %s -T IndelRealigner -I %s -o %s -targetIntervals %s -known %s -known %s -R %s -allowPotentiallyMisencodedQuals -rf NotPrimaryAlignment 2> %s && %s && %s" % (gatk,t_bam,t_bam_out,t_intval_out,millsAndIndels,dbSNP,ref,t_realn_log,t_index_cmd,t_realn_done_cmd)

  n_oi = sep.join([n_realn_done_log,n_bam_out,':',n_bam,n_intval_done_log])
  t_oi = sep.join([t_realn_done_log,t_bam_out,':',t_bam,t_intval_done_log])
  mf.write("%s\n\t%s\n\n" % (n_oi,n_realn_cmd))
  mf.write("%s\n\t%s\n" % (t_oi,t_realn_cmd))
  return([n_bam_out,t_bam_out,n_realn_done_log,t_realn_done_log])


#def BQSR(n_bam=None,nMark=None,t_bam=None,tMark=None,outdir='.',gatk=None,qcsite=None,nqcout=None,tqcout=None,qcwarn2=None,ref=None,millsAndIndels=None,dbSNP=None,bedFlank=None,name='sample',mf=None):
def BQSR(n_bam=None,nMark=None,t_bam=None,tMark=None,outdir='.',gatk=None,ref=None,millsAndIndels=None,dbSNP=None,bedFlank=None,name='sample',mf=None):
  '[7] Also for GATK4'
  if n_bam!=None:
    bn=basename(n_bam)
  if t_bam!=None:
    bn=basename(t_bam)

  if 'aln' not in bn:
    n_bam_out = os.path.join(outdir,'normal',dir5,'recal.bam')
    t_bam_out = os.path.join(outdir,'tumor',dir5,'recal.bam')
    n_bai_out = os.path.join(outdir,'normal',dir5,'recal.bai')
    t_bai_out = os.path.join(outdir,'tumor',dir5,'recal.bai')
    n_tab_out = os.path.join(outdir,'normal',dir5,'recal.table')
    t_tab_out = os.path.join(outdir,'tumor',dir5,'recal.table')
    n_bqsr_log = os.path.join(outdir,'normal',dir5,'recal.log')
    t_bqsr_log = os.path.join(outdir,'tumor',dir5,'recal.log')
    n_bqsr_done_log = os.path.join(outdir,'normal',dir5,'bqsr_done.log')
    t_bqsr_done_log = os.path.join(outdir,'tumor',dir5,'bqsr_done.log')
  else:
    n_bam_out = os.path.join(outdir,'normal',dir5,'recal.aln.bam')
    t_bam_out = os.path.join(outdir,'tumor',dir5,'recal.aln.bam')
    n_tab_out = os.path.join(outdir,'normal',dir5,'recal.aln.table')
    t_tab_out = os.path.join(outdir,'tumor',dir5,'recal.aln.table')
    n_bqsr_log = os.path.join(outdir,'normal',dir5,'recal.aln.log')
    t_bqsr_log = os.path.join(outdir,'tumor',dir5,'recal.aln.log')
    n_bqsr_done_log = os.path.join(outdir,'normal',dir5,'bqsr_done.aln.log')
    t_bqsr_done_log = os.path.join(outdir,'tumor',dir5,'bqsr_done.aln.log')

  if 'GenomeAnalysisTK' in gatk:  # gatk3
    n_idx_out = os.path.join(outdir,'normal',dir5,'recal.bai')
    t_idx_out = os.path.join(outdir,'tumor',dir5,'recal.bai')
    n_printR_log = os.path.join(outdir,'normal',dir5,'printR.log')
    t_printR_log = os.path.join(outdir,'tumor',dir5,'printR.log')
    n_idx_log = os.path.join(outdir,'normal',dir5,'index.log')
    t_idx_log = os.path.join(outdir,'tumor',dir5,'index.log')
    n_bqsr_cmd = 'java -jar %s -T BaseRecalibrator -I %s -o %s -R %s -knownSites %s -knownSites %s -nct 8 -allowPotentiallyMisencodedQuals -L %s 2> %s' % (gatk,n_bam,n_tab_out,ref,millsAndIndels,dbSNP,bedFlank,n_bqsr_log) # &> -> 2> (for ubuntu)
    t_bqsr_cmd = 'java -jar %s -T BaseRecalibrator -I %s -o %s -R %s -knownSites %s -knownSites %s -nct 8 -allowPotentiallyMisencodedQuals -L %s 2> %s' % (gatk,t_bam,t_tab_out,ref,millsAndIndels,dbSNP,bedFlank,t_bqsr_log)
    n_printR_cmd = 'java -jar %s -T PrintReads -I %s -BQSR %s -o %s -R %s -nct 8 -allowPotentiallyMisencodedQuals 2> %s' % (gatk,n_bam,n_tab_out,n_bam_out,ref,n_printR_log)
    t_printR_cmd = 'java -jar %s -T PrintReads -I %s -BQSR %s -o %s -R %s -nct 8 -allowPotentiallyMisencodedQuals 2> %s' % (gatk,t_bam,t_tab_out,t_bam_out,ref,t_printR_log)
    n_index_cmd = 'samtools index %s %s 2> %s' % (n_bam_out,n_idx_out,n_idx_log)
    t_index_cmd = 'samtools index %s %s 2> %s' % (t_bam_out,t_idx_out,t_idx_log)
    n_bqsr_done_cmd = 'echo "%s" > %s' % (name+'_n_bqsr_done',n_bqsr_done_log)
    t_bqsr_done_cmd = 'echo "%s" > %s' % (name+'_t_bqsr_done',t_bqsr_done_log)
    n_cmd = '%s && %s && %s && %s' % (n_bqsr_cmd,n_printR_cmd,n_index_cmd,n_bqsr_done_cmd)
    t_cmd = '%s && %s && %s && %s' % (t_bqsr_cmd,t_printR_cmd,t_index_cmd,t_bqsr_done_cmd)
  else:  # gatk4
    if 'aln' not in bn:
      n_apply_log = os.path.join(outdir,'normal',dir5,'apply.log')
      t_apply_log = os.path.join(outdir,'tumor',dir5,'apply.log')
    else:
      n_apply_log = os.path.join(outdir,'normal',dir5,'apply.aln.log')
      t_apply_log = os.path.join(outdir,'tumor',dir5,'apply.aln.log')
    # (BQSR) For single sample
    if n_bam!=None:
      n_bqsr_cmd = sep.join([gatk,'BaseRecalibrator -I',n_bam,'-O',n_tab_out,'-R',ref,'--known-sites',dbSNP,'--known-sites',millsAndIndels,'-L',bedFlank,'2>',n_bqsr_log]) # &> -> 2> (for ubuntu)
      n_apply_cmd = sep.join([gatk,'ApplyBQSR -I',n_bam,'-O',n_bam_out,'-R',ref,'-L',bedFlank,'-bqsr',n_tab_out,'2>',n_apply_log])
      n_ln_cmd='ln -sf %s %s.bai' % (n_bai_out,n_bam_out) # -s -> -sf (for ubuntu)
      n_bqsr_done_cmd = 'echo "%s" > %s' % (name+'_n_bqsr_done',n_bqsr_done_log)
      n_cmd = '%s && %s && %s && %s' % (n_bqsr_cmd,n_apply_cmd,n_ln_cmd,n_bqsr_done_cmd)
    if t_bam!=None:
      t_bqsr_cmd = sep.join([gatk,'BaseRecalibrator -I',t_bam,'-O',t_tab_out,'-R',ref,'--known-sites',dbSNP,'--known-sites',millsAndIndels,'-L',bedFlank,'2>',t_bqsr_log])
      t_apply_cmd = sep.join([gatk,'ApplyBQSR -I',t_bam,'-O',t_bam_out,'-R',ref,'-L',bedFlank,'-bqsr',t_tab_out,'2>',t_apply_log])
      t_ln_cmd='ln -sf %s %s.bai' % (t_bai_out,t_bam_out)
      t_bqsr_done_cmd = 'echo "%s" > %s' % (name+'_t_bqsr_done',t_bqsr_done_log)
      t_cmd = '%s && %s && %s && %s' % (t_bqsr_cmd,t_apply_cmd,t_ln_cmd,t_bqsr_done_cmd)
    # (siteQC) For dual samples
    #if n_bam!=None and t_bam!=None:
    #  qcdir=os.path.join(outdir,'qc','sites')
    #  assert not os.makedirs(qcdir,exist_ok=True) ,'Cannot create %s' % qcdir
    #  sum=os.path.join(qcdir,'summary.xls')
    #  log=os.path.join(qcdir,'qcsites.err')
    #  sqc_cmd=sep.join(['sh',qcsite,n_bam_out,t_bam_out,qcdir,sum,'2>',log])
    #  if 'GS519' in name:
    #    log=os.path.join(qcdir,'qcwarn2.err')
    #    qcwarn2_cmd=sep.join(['ruby',qcwarn2,name,nqcout,tqcout,sum,'2>',log])

  mf.write('\n\n.MAKEFLOW CATEGORY BQSR\n')
  if n_bam!=None:
    n_oi=sep.join([n_bqsr_done_log,n_bam_out,':',n_bam,nMark])
    mf.write('\n%s\n\t%s\n' % (n_oi,n_cmd))
  if t_bam!=None:
    t_oi=sep.join([t_bqsr_done_log,t_bam_out,':',t_bam,tMark])
    mf.write('\n%s\n\t%s\n' % (t_oi,t_cmd))
  #if n_bam!=None and t_bam!=None:
  #  oi=sep.join([sum,':',n_bam_out,t_bam_out])
  #  if 'GS519' in name:
  #    mf.write('\n%s\n\t%s\n' % (oi,sep.join([sqc_cmd,'&&',qcwarn2_cmd])))
  #  else:
  #    mf.write('\n%s\n\t%s\n' % (oi,sep.join([sqc_cmd])))

  if n_bam!=None and t_bam!=None:
    return([n_bam_out,t_bam_out,n_bqsr_done_log,t_bqsr_done_log])
  if n_bam==None:
    return([None,t_bam_out,None,t_bqsr_done_log,None])
  else:
    return([n_bam_out,None,n_bqsr_done_log,None,None])


def Mutect2(n_bam=None,t_bam=None,nMark=None,tMark=None,out=None,mf=None,nsplit=1,gatk=None,samtools='samtools',outdir='.',ref=None,exac=None,cosmic=None,dbSNP=None,name='Sample',bedFlank=None):
  '[8] Also for GATK4'
  vardir=os.path.join(outdir,'variants','somatic',dirs)
  mutect2_log = os.path.join(vardir,'mutect2.log')
  if n_bam!=None:
    ndir=os.path.dirname(n_bam)
  tdir=os.path.dirname(t_bam)
  nsplit=int(nsplit)

  mf.write("\n\n.MAKEFLOW CATEGORY Mutect2\n")

  if 'GenomeAnalysisTK' in gatk:  # gatk3
    mutect2_vcf_out = out[0]
    mutect2_done_log = out[1]
    done_cmd = 'echo "%s" > %s' % (name+"_mutect2_done",mutect2_done_log)
    if n_bam==None:
      outin=sep.join([mutect2_done_log,mutect2_vcf_out,':',tMark])
      cmd='java -jar %s -T MuTect2 -nct 20 -R %s -I:tumor %s -L %s --dbsnp %s --cosmic %s -o %s 2> %s && %s' % (gatk,ref,t_bam,bedFlank,dbSNP,cosmic,mutect2_vcf_out,mutect2_log,done_cmd) # &> -> 2> (for ubuntu)
    else:
      outin=sep.join([mutect2_done_log,mutect2_vcf_out,':',nMark,tMark])
      cmd="java -jar %s -T MuTect2 -nct 20 -R %s -I:tumor %s -I:normal %s -L %s --dbsnp %s --cosmic %s -o %s 2> %s && %s" % (gatk,ref,t_bam,n_bam,bedFlank,dbSNP,cosmic,mutect2_vcf_out,mutect2_log,done_cmd)
    mf.write("%s\n\t%s\n" % (outin,cmd))
  else:  # gatk4
    mutect2_vcf_out = os.path.join(vardir,'raw.vcf')

    if nsplit==1:
      mutect2_done_log=os.path.join(vardir,'raw.log')
      done_cmd='echo "%s" > %s' % (name+'_mutect2_done',mutect2_done_log)
      if n_bam==None:
        outin=sep.join([mutect2_done_log,mutect2_vcf_out,':',tMark])
        cmd=sep.join([gatk,'Mutect2','-R',ref,'-I',t_bam,'-tumor',name+'_tumor -L',bedFlank,'-O',mutect2_vcf_out,'--germline-resource',dbSNP,'--af-of-alleles-not-in-resource 0.0000025 --max-reads-per-alignment-start 1000 --max-suspicious-reads-per-alignment-start 6 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>',mutect2_log,'&&',done_cmd])
        mf.write("%s\n\t%s\n" % (outin,cmd))
      else:
        outin=sep.join([mutect2_done_log,mutect2_vcf_out,':',nMark,tMark])
        cmd=sep.join([gatk,'Mutect2','-R',ref,'-I',n_bam,'-normal',name+'_normal -I',t_bam,'-tumor',name+'_tumor -L',bedFlank,'-O',mutect2_vcf_out,'--germline-resource',dbSNP,'--af-of-alleles-not-in-resource 0.0000025 --max-reads-per-alignment-start 1000 --max-suspicious-reads-per-alignment-start 6 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>',mutect2_log,'&&',done_cmd])
        mf.write("%s\n\t%s\n" % (outin,cmd))
    elif nsplit==4 or nsplit==12:
      if n_bam!=None:
        nsdir=ndir+'/split'
      tsdir=tdir+'/split'
      vsdir=vardir+'/split'
      if n_bam!=None:
        assert not os.makedirs(nsdir,exist_ok=True) ,'Cannot create %s' % nsdir
      assert not os.makedirs(tsdir,exist_ok=True) ,'Cannot create %s' % tsdir
      assert not os.makedirs(vsdir,exist_ok=True) ,'Cannot create %s' % vsdir
      chrs=split_chr[nsplit]
      nbams=[]
      tbams=[]
      vcfs=[]
      logs=[]
      for i in range(0,nsplit):
        chr=sep.join(chrs[i])

        if n_bam!=None:
          nbam=nsdir+'/'+str(i)+'.bam'
          log=nsdir+'/'+str(i)+'.err'
          nbams.append(nbam)
          ndone_log=nsdir+'/'+str(i)+'_extract_done.log'
          done_cmd='echo "%s" > %s' % (str(i)+'.'+name+'_n_extract_done',ndone_log)
          oi=sep.join([nbam,ndone_log,':',n_bam])
          cmd=sep.join([samtools,'view -b',n_bam,chr,'>',nbam,'2>',log,'&&',samtools,'index',nbam,'2>>',log,'&&',done_cmd])
          mf.write("\n%s\n\t%s\n" % (oi,cmd))

        tbam=tsdir+'/'+str(i)+'.bam'
        log=tsdir+'/'+str(i)+'.err'
        tbams.append(tbam)
        tdone_log=tsdir+'/'+str(i)+'_extract_done.log'
        done_cmd='echo "%s" > %s' % (str(i)+'.'+name+"_t_extract_done",tdone_log)
        oi=sep.join([tbam,tdone_log,':',t_bam])
        cmd=sep.join([samtools,'view -b',t_bam,chr,'>',tbam,'2>',log,'&&',samtools,'index',tbam,'2>>',log,'&&',done_cmd])
        mf.write("\n%s\n\t%s\n" % (oi,cmd))

        vcf=vsdir+'/'+str(i)+'.vcf'
        vcfs.append(vcf+'.gz')
        log=vsdir+'/'+str(i)+'.log'
        done_log=vsdir+'/'+str(i)+'_mutect2_done.log'
        logs.append(done_log)
        done_cmd='echo "%s" > %s' % (str(i)+'.'+name+'_mutect2_done',done_log)
        if n_bam==None:
          oi=sep.join([vcf+'.gz',done_log,':',tbam,tdone_log])
          cmd=sep.join([gatk,'Mutect2','-R',ref,'-I',tbam,'-tumor',name+'_tumor -L',bedFlank,'-O',vcf,'--germline-resource',dbSNP,'--af-of-alleles-not-in-resource 0.0000025 --max-reads-per-alignment-start 1000 --max-suspicious-reads-per-alignment-start 6 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>',log,'&& bgzip',vcf,'&& tabix -fp vcf',vcf+'.gz &&',done_cmd])
        else:
          oi=sep.join([vcf+'.gz',done_log,':',nbam,tbam,ndone_log,tdone_log])
          cmd=sep.join([gatk,'Mutect2','-R',ref,'-I',nbam,'-normal',name+'_normal -I',tbam,'-tumor',name+'_tumor -L',bedFlank,'-O',vcf,'--germline-resource',dbSNP,'--af-of-alleles-not-in-resource 0.0000025 --max-reads-per-alignment-start 1000 --max-suspicious-reads-per-alignment-start 6 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter 2>',log,'&& bgzip',vcf,'&& tabix -fp vcf',vcf+'.gz &&',done_cmd])
        mf.write("\n%s\n\t%s\n" % (oi,cmd))
    else:
      print("Error [mutect2]: The number of Mutect2 threads '%d' is not 1, 4 or 12." % nsplit, file=sys.stderr)
      sys.exit(1)

    if nsplit==1:
      ''
    elif nsplit==4 or nsplit==12:
      mvcf=sep.join(vcfs)
      mlog=sep.join(logs)
      log=vardir+'/merge.log'
      done_log=vardir+'/merge_done.log'
      done_cmd='echo "%s" > %s' % (name+'_merge_done',done_log)
      oi=sep.join([mutect2_vcf_out,done_log,':',mvcf,mlog])
      cmd=sep.join(['vcf-concat -s 100',mvcf,'>',mutect2_vcf_out,'2>',log,'&&',gatk,'IndexFeatureFile -F',mutect2_vcf_out,'&&',done_cmd])
      mf.write("\n%s\n\t%s\n" % (oi,cmd))
    else:
      print("Error [merge]: The number of Mutect2 threads '%d' is not 1, 4 or 12." % nsplit, file=sys.stderr)
      sys.exit(1)

    qci_vcf=os.path.join(vardir,'raw.QCI.vcf')
    done_log=os.path.join(vardir,'qci.log')
    done_cmd='echo "%s" > %s' % (name+'_qci_done',done_log)
    outin=sep.join([qci_vcf,done_log,':',mutect2_vcf_out])
    cmd=sep.join(['sed "2i ##reference=NCBIb37t"',mutect2_vcf_out,'>',qci_vcf,'&&',done_cmd])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

    if n_bam!=None:
      n_pile_tab=os.path.join(outdir,'normal',dir5,'pile.tab')
      n_pile_log=os.path.join(outdir,'normal',dir5,'pile.log')
      n_pile_cmd=sep.join([gatk,'GetPileupSummaries -I',n_bam,'-L',bedFlank,'-O',n_pile_tab,'-V',exac,'2>',n_pile_log]) # &> -> 2> (for ubuntu)
      outin=sep.join([n_pile_tab,':',nMark,n_bam])
      mf.write("\n%s\n\t%s\n" % (outin,n_pile_cmd))

    t_pile_tab=os.path.join(outdir,'tumor',dir5,'pile.tab')
    t_pile_log=os.path.join(outdir,'tumor',dir5,'pile.log')
    t_pile_cmd=sep.join([gatk,'GetPileupSummaries -I',t_bam,'-L',bedFlank,'-O',t_pile_tab,'-V',exac,'2>',t_pile_log])
    outin=sep.join([t_pile_tab,':',tMark,t_bam])
    mf.write("\n%s\n\t%s\n" % (outin,t_pile_cmd))

    contab=os.path.join(vardir,'contamin.tab')
    calcont_log=os.path.join(vardir,'contamin.log')
    if n_bam==None:
      calcont_cmd=sep.join([gatk,'CalculateContamination -I',t_pile_tab,'-O',contab,'2>',calcont_log])
      outin=sep.join([contab,':',t_pile_tab])
    else:
      calcont_cmd=sep.join([gatk,'CalculateContamination -I',t_pile_tab,'-matched',n_pile_tab,'-O',contab,'2>',calcont_log])
      outin=sep.join([contab,':',t_pile_tab,n_pile_tab])
    mf.write("\n%s\n\t%s\n" % (outin,calcont_cmd))

    filter_out=os.path.join(vardir,'filteredCon.vcf')
    log=os.path.join(vardir,'FilterMutect.log')
    cmd=sep.join([gatk,'FilterMutectCalls -V',mutect2_vcf_out,'-O',filter_out,'-L',bedFlank,'-contamination-table',contab,'2>',log])
    outin=sep.join([filter_out,':',contab,mutect2_vcf_out])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

    artifact_out=os.path.join(vardir,'artifact.pre_adapter_detail_metrics.txt')
    log=os.path.join(vardir,'getArtifact.log')
    cmd=sep.join([gatk,'CollectSequencingArtifactMetrics -R',ref,'-I',t_bam,'-O',vardir+'/artifact','--FILE_EXTENSION ".txt" 2>',log])
    outin=sep.join([artifact_out,':',t_bam,tMark])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

    foxg_out=os.path.join(vardir,'filteredCon2.vcf')
    log=os.path.join(vardir,'FilterOxoG.log')
    done_log=os.path.join(vardir,'FilterOxoG_done.log')
    done_cmd='echo "%s" > %s' % (name+'_FilterOxoG_done',done_log)
    cmd=sep.join([gatk,'FilterByOrientationBias -V',filter_out,'-O',foxg_out,'-P',artifact_out,'-L',bedFlank,"--artifact-modes 'G/T' 2>",log,'&&',done_cmd])  # -AM 'C/T' (for FFPE samples)
    outin=sep.join([done_log,foxg_out,':',filter_out,artifact_out])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

    out=os.path.join(vardir,'filteredCon2.split.vcf')
    log=os.path.join(vardir,'splitMultiAllele.log')
    done_log=os.path.join(vardir,'splitMultiAllele_done.log')
    done_cmd='echo "%s" > %s' % (name+'_splitMultiAllele_done',done_log)
    cmd=sep.join([gatk,'LeftAlignAndTrimVariants --split-multi-allelics -V',foxg_out,'-O',out,'-R',ref,'2>',log,'&&',done_cmd])
    outin=sep.join([done_log,out,':',foxg_out])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

    mutect2_vcf_out=out
    mutect2_done_log=done_log

  return([mutect2_vcf_out,mutect2_done_log])


def Germline(n_bam,nMark,t_bam=None,tMark=None,mf=None,gatk=None,oneKGp3=None,outdir='.',ref=None,name='Sample',bed=None):
  '[E] Also for GATK4'
  if 'GenomeAnalysisTK' in gatk:  # gatk3
    print('Germline SNP/INDEL calling by GATK3 has not been implemented yet.', file=sys.stderr)
    sys.exit(1)

  vardir=os.path.join(outdir,'variants','germline',dirs)
  mf.write("\n\n.MAKEFLOW CATEGORY Germline_SNP_INDEL\n")

  nvcf=os.path.join(vardir,'normal.g.vcf')
  log=os.path.join(vardir,'HaploCallNormal.log')
  done_log=os.path.join(vardir,'HaploCallNormal_done.log')
  done_cmd='echo "%s" > %s' % (name+'_HaplotypeCaller_Normal_done',done_log)
  cmd=sep.join([gatk,'HaplotypeCaller -R',ref,'-I',n_bam,'-L',bed,'-O',nvcf,'-ERC GVCF 2>',log,'&&',done_cmd]) # &> -> 2> (for ubuntu)
  outin=sep.join([nvcf,done_log,':',n_bam,nMark])
  mf.write("%s\n\t%s\n" % (outin,cmd))

  if t_bam!=None:
    tvcf=os.path.join(vardir,'tumor.g.vcf')
    log=os.path.join(vardir,'HaploCallTumor.log')
    done_log=os.path.join(vardir,'HaploCallTumor_done.log')
    done_cmd='echo "%s" > %s' % (name+'_HaplotypeCaller_Tumor_done',done_log)
    cmd=sep.join([gatk,'HaplotypeCaller -R',ref,'-I',t_bam,'-L',bed,'-O',tvcf,'-ERC GVCF 2>',log,'&&',done_cmd])
    outin=sep.join([tvcf,done_log,':',t_bam,tMark])
    mf.write("\n%s\n\t%s\n" % (outin,cmd))

  cvcf=os.path.join(vardir,'combined.g.vcf')
  log=os.path.join(vardir,'CombGVCFs.log')
  done_log=os.path.join(vardir,'CombGVCFs_done.log')
  done_cmd='echo "%s" > %s' % (name+'_CombineGVCFs_done',done_log)
  if t_bam!=None:
    cmd=sep.join([gatk,'CombineGVCFs -R',ref,'-V',nvcf,'-V',tvcf,'-L',bed,'-O',cvcf,'2>',log,'&&',done_cmd])
    outin=sep.join([cvcf,done_log,':',nvcf,tvcf])
  else:
    cmd=sep.join([gatk,'CombineGVCFs -R',ref,'-V',nvcf,'-L',bed,'-O',cvcf,'2>',log,'&&',done_cmd])
    outin=sep.join([cvcf,done_log,':',nvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  rvcf=os.path.join(vardir,'raw.vcf')
  qci_vcf=os.path.join(vardir,'raw.QCI.vcf')
  log=os.path.join(vardir,'GenoGVCFs.log')
  done_log=os.path.join(vardir,'GenoGVCFs_done.log')
  done_cmd='echo "%s" > %s' % (name+'_GenotypeGVCFs_done',done_log)
  cmd=sep.join([gatk,'GenotypeGVCFs -R',ref,'-V',cvcf,'-L',bed,'-O',rvcf,'2>',log,'&& sed "2i ##reference=NCBIb37t"',rvcf,'>',qci_vcf,'&&',done_cmd])
  outin=sep.join([rvcf,done_log,':',cvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  pvcf=os.path.join(vardir,'posterior.vcf')
  log=os.path.join(vardir,'CalGenoPost.log')
  done_log=os.path.join(vardir,'CalGenoPost_done.log')
  done_cmd='echo "%s" > %s' % (name+'_CalculateGenotypePosteriors_done',done_log)
  cmd=sep.join([gatk,'CalculateGenotypePosteriors -R',ref,'-V',rvcf,'-L',bed,'-O',pvcf,'-supporting',oneKGp3,'2>',log,'&&',done_cmd])
  outin=sep.join([pvcf,done_log,':',rvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  svcf=os.path.join(vardir,'posteriorSNP.vcf')
  log=os.path.join(vardir,'SelVarSNP.log')
  done_log=os.path.join(vardir,'SelVarSNP_done.log')
  done_cmd='echo "%s" > %s' % (name+'_SelectVariants_SNP_done',done_log)
  cmd=sep.join([gatk,'SelectVariants -R',ref,'-V',pvcf,'-L',bed,'-O',svcf,'--select-type-to-include SNP','2>',log,'&&',done_cmd])
  outin=sep.join([svcf,done_log,':',pvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  sfvcf=os.path.join(vardir,'filteredSNP.vcf')
  log=os.path.join(vardir,'VarFiltSNP.log')
  done_log=os.path.join(vardir,'VarFiltSNP_done.log')
  done_cmd='echo "%s" > %s' % (name+'_VariantFiltration_SNP_done',done_log)
  cmd=sep.join([gatk,'VariantFiltration -R',ref,'-V',svcf,'-L',bed,'-O',sfvcf,'--filter-name "SNPfilter" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --genotype-filter-name "lowGQ" --genotype-filter-expression "GQ < 20"','2>',log,'&&',done_cmd])
  outin=sep.join([sfvcf,done_log,':',svcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  ivcf=os.path.join(vardir,'posteriorINDEL.vcf')
  log=os.path.join(vardir,'SelVarINDEL.log')
  done_log=os.path.join(vardir,'SelVarINDEL_done.log')
  done_cmd='echo "%s" > %s' % (name+'_SelectVariants_INDEL_done',done_log)
  cmd=sep.join([gatk,'SelectVariants -R',ref,'-V',pvcf,'-L',bed,'-O',ivcf,'--select-type-to-include INDEL','2>',log,'&&',done_cmd])
  outin=sep.join([ivcf,done_log,':',pvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  ifvcf=os.path.join(vardir,'filteredINDEL.vcf')
  log=os.path.join(vardir,'VarFiltINDEL.log')
  done_log=os.path.join(vardir,'VarFiltINDEL_done.log')
  done_cmd='echo "%s" > %s' % (name+'_VariantFiltration_INDEL_done',done_log)
  cmd=sep.join([gatk,'VariantFiltration -R',ref,'-V',ivcf,'-L',bed,'-O',ifvcf,'--filter-name "INDELfilter" --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --genotype-filter-name "lowGQ" --genotype-filter-expression "GQ < 20"','2>',log,'&&',done_cmd])
  outin=sep.join([ifvcf,done_log,':',ivcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  mout=os.path.join(vardir,'filtered.vcf')
  log=os.path.join(vardir,'MergeVcfs.log')
  done_log=os.path.join(vardir,'MergeVcfs_done.log')
  done_cmd='echo "%s" > %s' % (name+'_MergeVcfs_done',done_log)
  cmd=sep.join([gatk,'MergeVcfs -I',sfvcf,'-I',ifvcf,'-O',mout,'2>',log,'&&',done_cmd])
  outin=sep.join([mout,done_log,':',sfvcf,ifvcf])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  out=os.path.join(vardir,'filtered.split.vcf')
  log=os.path.join(vardir,'splitMultiAllele.log')
  done_log=os.path.join(vardir,'splitMultiAllele_done.log')
  done_cmd='echo "%s" > %s' % (name+'_splitMultiAllele_done',done_log)
  cmd=sep.join([gatk,'LeftAlignAndTrimVariants --split-multi-allelics -V',mout,'-O',out,'-R',ref,'2>',log,'&&',done_cmd])
  outin=sep.join([done_log,out,':',mout])
  mf.write("\n%s\n\t%s\n" % (outin,cmd))

  return([out,done_log])


#def Annotate(vcf,passvcf,mark,passmark,xref,mf=None,annovar=None,outdir='.',anndb=None,name='Sample'):
def Annotate(vcf=None,mark=None,n_bam=None,t_bam=None,xref=None,title=None,vep=None,dir=None,ref=None,mf=None,annovar=None,outdir='.',anndb=None,name='Sample',fpfilter=None,cosmver=87):
  '[9][j] Annotate by Annovar and VEP. Also for GATK4'
  if 'somatic' in vcf:
    anndir=os.path.join(outdir,'variants','somatic',dira)
  elif 'germline' in vcf:
    anndir=os.path.join(outdir,'variants','germline',dira)
  else:
    print("Unknown type of vcf '%s'." % vcf, file=sys.stderr)
    sys.exit(1)

  #pcmd='-buildver hg19 -remove -operation g,r,r,f,f,f,f -nastring . -vcfinput -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,snp138,ljb26_all -out'
  pcmd='-xref %s -buildver hg19 -remove -nastring . -vcfinput -operation gx,r,r,f,f,f,f,f,f,f,f,f,f -protocol refGeneWithVer,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,snp142,snp138NonFlagged,dbnsfp33a,cosmic%s_coding,cosmic%s_noncoding,clinvar_20170905,clinvar20180429 -out' % (xref,cosmver,cosmver)  # intervar_20170202
  vcf_out=os.path.join(anndir,'annovar.hg19_multianno.vcf')
  run_log=os.path.join(anndir,'annovar.log')
  done_log=os.path.join(anndir,'annovar_done.log')
  txt_out=os.path.join(anndir,'annovar.hg19_multianno.txt')
  done='echo "%s" > %s' % (name+'_annovar_done',done_log)
  cmd=sep.join(['perl',annovar,vcf,anndb,pcmd,anndir+'/annovar','2>',run_log,'&&',done])
  outin=sep.join([vcf_out,txt_out,done_log,':',vcf,mark])
  mf.write('\n\n.MAKEFLOW CATEGORY AnnotationAnnovar_%s\n' % title)
  mf.write("%s\n\t%s\n" % (outin,cmd))

  mf.write('\n\n.MAKEFLOW CATEGORY AnnotationFpfilter_%s\n' % title)
  n_filter_vcf_out = os.path.join(anndir,'fpfilt.n.vcf')
  n_filter_xls_out = os.path.join(anndir,'fpfilt.n.xls')
  n_filter_vcf_log = os.path.join(anndir,'fpfilt.n.log')
  t_filter_vcf_out = os.path.join(anndir,'fpfilt.t.vcf')
  t_filter_xls_out = os.path.join(anndir,'fpfilt.t.xls')
  t_filter_vcf_log = os.path.join(anndir,'fpfilt.t.log')
  if n_bam!=None:
    n_filter_mf_cmd=sep.join(['perl',fpfilter,'--min-var-count 3 --min-var-freq 0.01 --vcf-file',vcf,'--bam-file',n_bam,'--sample',name+'_normal --reference',ref,'--output',n_filter_vcf_out,'2>',n_filter_vcf_log])
  #else:
  #  n_filter_mf_cmd=sep.join(['ln -sf',t_filter_vcf_out,n_filter_vcf_out,'&& ln -sf',t_filter_xls_out,n_filter_xls_out,'2>',n_filter_vcf_log])
    outin=sep.join([n_filter_vcf_out,':',vcf,mark])
    mf.write('%s\n\t%s\n' % (outin,n_filter_mf_cmd))
  if t_bam!=None:
    t_filter_mf_cmd=sep.join(['perl',fpfilter,'--min-var-count 3 --min-var-freq 0.01 --vcf-file',vcf,'--bam-file',t_bam,'--sample',name+'_tumor --reference',ref,'--output',t_filter_vcf_out,'2>',t_filter_vcf_log])
  #else:
  #  t_filter_mf_cmd=sep.join(['ln -sf',n_filter_vcf_out,t_filter_vcf_out,'&& ln -sf',n_filter_xls_out,t_filter_xls_out,'2>',t_filter_vcf_log])
    outin=sep.join([t_filter_vcf_out,':',vcf,mark])
    mf.write('\n%s\n\t%s\n' % (outin,t_filter_mf_cmd))

  cache_option='' if dir==None else '-dir '+dir+' '
  run_log=os.path.join(anndir,'vep.log')
  done_log=os.path.join(anndir,'vep_done.log')
  txt_out=os.path.join(anndir,'vepout.txt')
  done='echo "%s" > %s' % (name+"_vep_done",done_log)
  cmd=sep.join([vep,cache_option+'--no_escape','-i',vcf,'-o',txt_out,'--fasta',ref,'--verbose --no_progress --shift_hgvs 1 --everything --cache --offline --force_overwrite --use_given_ref --refseq --fork 3 --pick 2>',run_log,'&&',done])
  outin=sep.join([done_log,txt_out,':',vcf,mark])
  mf.write('\n\n.MAKEFLOW CATEGORY AnnotationVEP_%s\n' % title)
  mf.write("%s\n\t%s\n" % (outin,cmd))

  return([done_log,txt_out,vcf_out])

def callCNV(n_bam,t_bam,nMark,tMark,bed=None,mf=None,contra=None,outdir='.',name='Sample'):
  '[c]'
  odir=os.path.join(outdir,'variants','somatic',dirc)
  cnv_log=os.path.join(odir,'callcnv.log')
  done_log=os.path.join(odir,'cnvdone.log')
  cnv_out=os.path.join(odir,'table',name+'.CNATable.10rd.10bases.20bins.DetailsFILTERED.txt')

  cnv_mark_success_cmd="echo \"%s\" > %s" % (name+"_callCNV_done",done_log)
  cnv_mf_cmd="%s --target %s --test %s --control %s --outFolder %s --sampleName %s --nomultimapped --removeDups --plot 2> %s && %s" % (contra,bed,t_bam,n_bam,odir,name,cnv_log,cnv_mark_success_cmd)

  cnv_mf_oi=sep.join([done_log,cnv_out,':',n_bam,t_bam,nMark,tMark])
  mf.write("\n\n.MAKEFLOW CATEGORY callCNV\n")
  mf.write("%s\n\t%s\n" % (cnv_mf_oi,cnv_mf_cmd))
  return([cnv_out,done_log])

def TMBold(vcf,passvcf,mark,passmark,bed82,bedtools,bed=None,mf=None,outdir='.',name='Sample'):
  '[t]'
  tmbdir=os.path.join(outdir,'variants','somatic',dira)
  mf.write("\n\n.MAKEFLOW CATEGORY TumorMutationBurden\n")

  done_log=os.path.join(tmbdir,"tmb_ori_done.log")
  done_cmd="echo \"%s\" > %s" % (name+"_tmb_ori_done",done_log)
  outfile=os.path.join(tmbdir,"ori.tmb")
  vcfnew=os.path.join(tmbdir,"ori.82g.vcf")
  shf=os.path.join(tmbdir,"tmb_ori.sh")
  sh=open(shf,'w')
  sh.write("%s intersect -a %s -b %s > %s\nn=`grep -v -c '^#' %s`\nawk -v n=$n '{s+=($3-$2)}END{print n/s*1e6}' %s > %s\n" % (bedtools,vcf,bed82,vcfnew,vcfnew,bed,outfile))
  sh.close()
  cmd="sh %s && %s" % (shf,done_cmd)
  outin=sep.join([outfile,done_log,':',vcf,mark])
  mf.write("%s\n\t%s\n\n" % (outin,cmd))

  done_log2=os.path.join(tmbdir,"tmb_pass_done.log")
  done_cmd="echo \"%s\" > %s" % (name+"_tmb_psss_done",done_log2)
  outfile2=os.path.join(tmbdir,"pass.tmb")
  passvcfnew=os.path.join(tmbdir,"pass.82g.vcf")
  shf=os.path.join(tmbdir,"tmb_pass.sh")
  sh=open(shf,'w')
  sh.write("%s intersect -a %s -b %s > %s\nn=`grep -v -c '^#' %s`\nawk -v n=$n '{s+=($3-$2)}END{print n/s*1e6}' %s > %s\n" % (bedtools,passvcf,bed82,passvcfnew,passvcfnew,bed,outfile2))
  sh.close()
  cmd="sh %s && %s" % (shf,done_cmd)
  outin=sep.join([outfile2,done_log2,':',passvcf,passmark])
  mf.write("%s\n\t%s\n" % (outin,cmd))
  return([outfile,outfile2,done_log,done_log2])

def DelBigLess(anndone1,anndone2,fq_u,bam2_n,bam2_t,bam4_n,bam4_t,bqsr_mark_n,bqsr_mark_t,mf=None,outdir='.',name='Sample'):
#def DelBigLess(anndone1,anndone2,fq_nclean,fq_tclean,fq_u,bam2_n,bam2_t,bam3_n,bam3_t,bam4_n,bam4_t,bqsr_mark_n,bqsr_mark_t,mf=None,outdir='.',name='Sample'):
  '[y] Keep clean data, sorted bams, deduped bams'
  mf.write("\n\n.MAKEFLOW CATEGORY DeleteBigFilesLess\n")
  rmfiles=sep.join([anndone1,anndone2,fq_u,bam2_n,bam2_t,bam4_n,bam4_t])
  #rmfiles=sep.join([anndone1,anndone2,fq_nclean,fq_tclean,fq_u,bam2_n,bam2_t,bam3_n,bam3_t,bam4_n,bam4_t])
  donefile=os.path.join(outdir,"del_big_files_done.log")
  done="echo \"%s\" >> %s" % (name+"_del_big_files_done",donefile)
  cmd="ls -l %s > %s && rm -f %s && %s" % (rmfiles,donefile,rmfiles,done)
  outin=sep.join([donefile,':',rmfiles,bqsr_mark_n,bqsr_mark_t])
  mf.write("%s\n\t%s\n" % (outin,cmd))

def DelBig(anndone1,anndone2,fq_nclean,fq_tclean,fq_u,bam2_n,bam2_t,bam2_sort_n,bam2_sort_t,bam3_n,bam3_t,bam4_n,bam4_t,bqsr_mark_n,bqsr_mark_t,mf=None,outdir='.',name='Sample'):
  '[z]'
  mf.write("\n\n.MAKEFLOW CATEGORY DeleteBigFiles\n")
  rmfiles=sep.join([anndone1,anndone2,fq_nclean,fq_tclean,fq_u,bam2_n,bam2_t,bam2_sort_n,bam2_sort_t,bam3_n,bam3_t,bam4_n,bam4_t])
  donefile=os.path.join(outdir,"del_big_files_done.log")
  done="echo \"%s\" >> %s" % (name+"_del_big_files_done",donefile)
  cmd="ls -l %s > %s && rm -f %s && %s" % (rmfiles,donefile,rmfiles,done)
  outin=sep.join([donefile,':',rmfiles,bqsr_mark_n,bqsr_mark_t])
  mf.write("%s\n\t%s\n" % (outin,cmd))
