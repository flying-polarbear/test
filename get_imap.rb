#! /usr/bin/env ruby
#*****************************************************************************
#                Copyright (C) 2019, gaoyuhui. All Rights Reserved
#
#  FileName:   get_imap.rb
#  Desc:       
#  Author:     gaoyuhui
#  Email:      gaoyuhui@genesmile.com
#  LastChange: 2019-01-03 18:20:30
#  History:    2019-01-03 Create by gaoyuhui
#*****************************************************************************
require 'mail'
require 'json'
require 'optparse'

version = 'v1.0_2019_01_07'
OPTIONS = {}
optparse = OptionParser.new do |opts|
  opts.banner = "Usage: #{__FILE__} [options]"
  opts.on('-r', '--run', "run 519b"){|e| OPTIONS[:run] = TRUE }
  opts.on('-a', '--runa', "run blood"){|e| OPTIONS['blood'] = TRUE }
end.parse!
p OPTIONS

DIR = '/nfs/disk13/519_sample/novo_sample'

def work_dir oss
  File.join(DIR, File.basename(oss))
end
#work_dir = ->(oss){ File.join(DIR, File.basename(oss) }

def info info
  puts "[Info]: #{info}"
end
#info = ->(info){puts 'Info:' info}

info "Attention: main dir is #{DIR}"
info  "Starts ...\n\n"

t_gap = 60 * 30
#2019-01-02T12:56:17+08:00

def get_imap

  Mail.defaults do
    retriever_method :imap, :address    => "imap.exmail.qq.com",
                            :user_name  => 'autowork@genesmile.com',
                            :password   => 'Work1234'
                           # :password   => gpg -q -d /nfs/disk3/user/gaoyuhui/sop/autowork.ps`.chomp
  end
end

def recent(imap,t,t_gap)
  ###maybe it is empty arr
  imap.find(:what => :last, :keys => ['FROM','gaoyuhui','UNSEEN', 'SUBJECT', 'dl 519']).select{|s| t - Time.new(*s.date.to_s.scan(/\d+/).take(6)) < t_gap  } 
end

###[ [subject,oss1,oss2...],...   ]
def mail2link mails
#  mails.map{|m1| pp m1.body.to_s.split(/\r|\r\n|\n/),1 } 
  mails.map{|m1| [m1.subject] + m1.body.to_s.split.select{|s| s.match?('oss://')}.map{|m2| m2.match('oss://novo-medical-jinuosimei-nj/.*').to_s } }
end

def dl_oss oss
  `oss_auto_dl.rb LTAIuSRj69CrVT7V eCmoNs5VoKHKi2MazxqlwSOAVrR220 #{oss}`
end


#lambda dl_oss
# dl_oss = ->link{`oss_auto_dl.rb LTAIuSRj69CrVT7V eCmoNs5VoKHKi2MazxqlwSOAVrR220 #{link}`}
def md5chk(oss)
  #oss://novo-medical-jinuosimei-nj/20190112_HFF2JDMXX-P101SC18090028-01_Result
  #/nfs/disk22/sample/novo_sample/20190112_HFF2JDMXX-P101SC18090028-01_Result
  re_dl_arr = []
  md5_dir = work_dir(oss)
  if File.exist? File.join(md5_dir, 'MD5.txt')   ### MD5.txt must existed in md5_dir
    info "Checking MD5 of #{File.join(md5_dir, 'MD5.txt')}\n......"
    Dir.chdir(md5_dir) do 
      md5_res = `md5sum -c MD5.txt`
      info "md5 check result is : \n#{md5_res}"
      re_dl_arr = md5_res.split("\n").select{|s| s.end_with? 'FAILED'}
    end
  end
  re_dl_arr
end

def re_dl(oss,re_dl_arr)
  ###oss: oss://novo-medical-jinuosimei-nj/20190112_HFF2JDMXX-P101SC18090028-01_Result
  ###re_dl_arr=> ['Rawdata/WK519-190108-01-GS519-190104-02W-m_L1_1.fq.gz: FAILED',...]
  ### novo_single_dl.py 20190112_HFF2JDMXX-P101SC18090028-01_Result/MD5.txt
  dir1 = File.basename(oss)
  
  Dir.chdir(DIR) do 
    re_dl_arr.map{|m| File.join(dir1, m[/.*(?=:)/])}.each do |file| 
      info "re-downloading #{file}\n......"
      `novo_single_dl.py #{file}`
      info "#{file} re-downloaded"
    end
  end
end

def notice_done
  `/nfs/disk3/user/gaoyuhui/software/mail_for_notice.rb`
end

#lambda notice_done
# notice_done = ->{`/nfs/disk3/user/gaoyuhui/software/mail_for_notice.rb`}

def notice_dl oss
  subject = 'Start Downloading'
  content = "downloading #{oss}"
  to = 'gaoyuhui@genesmile.com'
  `/nfs/disk3/user/gaoyuhui/software/send_email.rb "#{subject}" "#{content}" "#{to}"`
end

def merge_fq(oss)
  #oss: oss://novo-medical-jinuosimei-nj/20190112_HFF2JDMXX-P101SC18090028-01_Result
  dir = File.join(work_dir(oss), 'Rawdata')
  info "make sh file for merge in dir #{dir}"
  `merge_fq.rb #{dir}`
  Dir.chdir "#{dir}/merge" do 
    info "merging in dir #{Dir.pwd}\n......"
    `sh merge.sh`
  end
end

def get_519_sample(oss, merge = '', opt = OPTIONS)
  pwd = File.join(work_dir(oss), 'Rawdata', merge)
  Dir.chdir(pwd) do
    info "Generate sample info file"
    `get_sample.rb #{Dir.pwd}`
    info "Submit to slurm......"
    `519b_run.rb -r #{Dir.pwd}/GS519_sample.txt` if opt[:run]
  end
end

def get_blood_sample(oss, merge = '')
  pwd = File.join(DIR, File.basename(oss), 'Rawdata', merge)
  Dir.chdir(pwd) do
    `get_blood_sample.rb #{Dir.pwd}`
  end
end

def dl_oss_link(oss)
  Dir.chdir(DIR) do
  info "downloading at dir posistion: #{DIR}"
    #oss: oss://novo-medical-jinuosimei-nj/20190112_HFF2JDMXX-P101SC18090028-01_Result
    unless Dir.exist? File.basename(oss)
      notice_dl oss
      dl_oss oss
      notice_done
    end
  end
end

def detail_dl oss, merge
  info "Starts downloading #{oss}\n......\n"
  dl_oss_link oss
  puts "Done\n\n"
  
  info "Check MD5 file ......"
  re_dl_arr = md5chk oss
  puts "Done\n\n"

  if re_dl_arr.empty?
    info "All files have passed the md5 chk\n\n"
  else
    info "The MD5 chk of #{re_dl_arr} didn't pass, MUST re-download"
    info "Re-downloading ......"
    re_dl(oss, re_dl_arr)
    info "All re-download done\n\n"
  end

  if merge.end_with? 'mg'
    info "merge fastq files\n......"
    merge_fq oss
    puts "Done\n\n"
  end

  info "Gets 519 sample info ......"
  case merge
  when '519'
    get_519_sample oss, merge
  when '519mg'
    get_519_sample oss, 'merge'
  when 'blood'
    get_blood_sample oss
  when 'bloodmg'
    get_blood_sample oss, 'merge'
  end
  puts "Done\n\n"
end

def analysis(mail_arr)
### mail_arr => [  [subject,oss1,oss2,oss3,...],... ]
  mail_arr.each do |d|
    case d.shift
    when 'dl 519'
      info "Mission is: dl 519\n\n"
      d.each{|oss|  detail_dl oss, '519'}
    when 'dl 519 mg'
      info "Mission is: dl 519 merge\n\n"
      d.each{|oss| detail_dl oss, '519mg'}
    when 'dl blood'
      info "Mission is: dl blood\n\n"
      d.each{|oss| detail_dl oss, 'blood'}
    when 'dl blood mg'
      info "Mission is: dl blood merge\n\n"
      d.each{|oss|  detail_dl oss, 'bloodmg'}
    end
  end
end

loop do
  t = Time.now
 
  info "Check mail at #{t}"
  imap = get_imap 
  mails = recent(imap, t, t_gap)
  mail_arr = mail2link(mails)

  unless mail_arr.empty?
    info "Processing #{mail_arr} from imap_mail"
    analysis(mail_arr)
  end
  info "Starts sleeping #{t_gap}s......"
  puts "......\n"
  STDOUT.flush
  sleep((1..6) === t.hour ? 60 * 60 * 2 : t_gap)
  #(1..6) === t.hour ? sleep(60 * 60 * 2) : sleep(t_gap)
  info "Continue check ......\n"
end
