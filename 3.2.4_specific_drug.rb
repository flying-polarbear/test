#! /nfs/disk3/user/gaoyuhui/.rvm/rubies/ruby-2.5.1/bin/ruby 
#*****************************************************************************
#                Copyright (C) 2018, gaoyuhui. All Rights Reserved
#
#  FileName:   specific_drug_13.rb
#  Desc:       
#  Author:     gaoyuhui
#  Email:      gaoyuhui@genesmile.com
#  LastChange: 2018-11-22 17:47:27
#  History:    2018-10-24 Create by gaoyuhui
#              2018-11-07 Modify by gaoyuhui
#              2018-11-09 Modify by gaoyuhui
#              2018-11-15 Modify by gaoyuhui
#              2018-11-15 Modify by gaoyuhui
#              2018-11-19 Modify by gaoyuhui
#              2018-11-22 Modify by gaoyuhui
#*****************************************************************************
require 'roo'

somatic_file, sample_cancer = ARGV 

db_file = '/nfs/disk3/user/gaoyuhui/sop/drug_mutation_database.xlsx'
db = Roo::Excelx.new( db_file )
#p all_s = x.sheets
db.default_sheet = 'Immunotherapy_database'

select = {
  :cancer     => 'Cancer_type',
  :sub_cancer => 'Cancer_subtype',
  :gene       => 'Marker',  #MSH6
  :p_mut      => 'Value',
  :evdc       => 'Evidence_level',
  :drug       => 'Drug',
  :effe       => 'Effect',
  :clinic     => 'Clinical_evidence', #A B C
}

db_les, db_big = {}, {}
dup = []
### get info from database
db.each( select ) do |sele|
  p sele
  next if sele[:marker] == 'Marker'
  cancer = sele[:cancer]
  sub_cancer = sele[:sub_cancer].chomp 
  gene  = sele[:gene].chomp.split
  p_mut = sele[:p_mut].split
  drug  = sele[:drug].chomp.split.join ','
  evdc  = sele[:evdc].chomp
  effe  = sele[:effe]
  clinic  = sele[:clinic].chomp

  dup << gene if gene.size >= 2
  if gene.size == 1
    key_big = [ cancer,gene.join(','),p_mut ].join('|')
    key_les = [ sub_cancer,gene.join(','),p_mut ].join('|')
  else
    key_big = [ cancer,gene.join(',') ].join('|')
    key_les = [ sub_cancer,gene.join(',') ].join('|')
  end
  db_big[ key_big ] = [ drug, effe, clinic.gsub(/\n/,'NEWLINE'), evdc ]
  db_les[ key_les ] = [ drug, effe, clinic.gsub(/\n/,'NEWLINE'), evdc ] if sub_cancer != '-'
end
p db_big,db_les,dup

loss_std = %w(frameshift splice stop_gained start_lost)
### should the same as cancer or sub_cancer in database
somatic = Roo::Excelx.new( somatic_file )
#p rate.sheets
somatic.default_sheet = 'somatic_final'

s_sele = {
  :gene  => 'Gene.refGeneWithVer', #MLH1
  :loss_std   => 'Consequence <VEP>',
  :p_mut => 'HGVSp3 <VEP>', #p.Val384Asp
  :exon  => 'EXON <VEP>',
  :intron => 'INTRON <VEP>',
  :v_type => 'VARIANT_CLASS <VEP>',
}
fh = File.open('16_specific_drug_3.2.4.txt','w+')
fh.puts "药物\t药效\t受检者药物疗效临床解析"
#肺腺癌|EGFR,KRAS|p.E746_A750del|p.G12A"
write = []
genes_all = []
p_mut_h = {}
somatic.each( s_sele ) do |r|
  next if r[:gene] == 'Gene.refGeneWithVer' || r[:gene].nil?
  gene  = r[:gene]
  l_std = r[:loss_std]
  p_mut = r[:p_mut]
  exon = r[:exon]
  intron = r[:intron]
  v_type = r[:v_type]
  key = [ sample_cancer,gene,p_mut ].join '|'
  genes_all << gene
  p_mut_h[gene] = p_mut
  if db_les.key?(key)
    write << db_les[key] + [gene, p_mut]
  elsif db_big.key?(key)
    write << db_big[key] + [gene, p_mut]
  elsif loss_std.any?{|e| e == l_std }
    key = [ sample_cancer,gene,'loss' ].join '|'
    db_les.key?(key) ? write << db_les[key] + [gene, p_mut] :
    db_big.key?(key) ? write << db_big[key] + [gene, p_mut] : write << "-\t-\t-"
  elsif gene == 'EGFR' && [ exon.split('/')[0], v_type ] == [ '19', 'deletion' ]
    key = [ sample_cancer,gene,'EX19del' ].join '|'
    db_les.key?(key) ? write << db_les[key] + [gene, p_mut] :
    db_big.key?(key) ? write << db_big[key] + [gene, p_mut] : write << "-\t-\t-"
  elsif v_type == 'cnv'
    key = [ sample_cancer,gene,'copy number gain' ].join '|'
    db_les.key?(key) ? write << db_les[key] + [gene, p_mut] :
    db_big.key?(key) ? write << db_big[key] + [gene, p_mut] : write << "-\t-\t-"
  elsif gene == 'MET' && (exon.split('/')[0] == '14' || intron.split('/')[0] == '13' )
    key = [ sample_cancer,gene,'Exon14' ].join '|'
    db_les.key?(key) ? write << db_les[key] + [gene, p_mut] :
    db_big.key?(key) ? write << db_big[key] + [gene, p_mut] : write << "-\t-\t-"
  else
    write << "-\t-\t-"
  end
end

###  [[ drug, effe, clinic.gsub(/\n/,'NEWLINE'), evdc, gene, p_mut ]... ]
dup_gene = dup.find_all{|a| a.all?{|b| genes_all.include? b } }
###!!!!!!!!! dup_gene may ERROR from overlap gene in dup_gene & genes_all
#
fh1 = File.open('add_for_immun_effec_2','w+')
if dup_gene.empty?
  dup_gene.each do |ea|
    dup_key = [sample_cancer,ea.join(',')].join '|'
    if db_les.key?(dup_key)
      fh.puts(db_les[dup_key].take(3).join "\t")
      ea.map{|m| p_mut_h[m]}.join
      fh1.puts [ea.join(','), ea.map{|m| p_mut_h[m]}.join(','), '-', db_les[dup_key][0,2], db_les[dup_key][-1] ].join("\t")
    elsif db_big.key?(dup_key)
      fh.puts(db_big[dup_key].take(3).join "\t")
      fh1.puts [ea.join(','), ea.map{|m| p_mut_h[m]}.join(','), '-', db_big[dup_key][0,2], db_big[dup_key][-1] ].join("\t")
    end 
  end
end

w_final = write - ["-\t-\t-"]
w_final.empty? ? fh.puts( "-\t-\t-" ) : w_final.sort_by{|so| so[3]}.each do |e|
  fh.puts e.take(3).join("\t")
  fh1.puts (e[-2,2] + ['-'] + e[0,2] + [e[3]]).join("\t")
end
fh.close
fh1.close
