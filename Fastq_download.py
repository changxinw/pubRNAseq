import os, sys
import urllib2
import re
# import commands
import argparse

#Set up the path to the fastq-dump here
fdmp = 'fastq-dump'
#Set up the path to salmon software
smn = '/data5/home/changxin/miniconda2/envs/salmon/bin/salmon'
#Set up the path to genome index
hg38_index = "/data5/home/changxin/genome_index/salmon_refseq_hg38"
mm10_index = "/data5/home/changxin/genome_index/salmon_refseq_mm10"

###This code still require the checking read length script


# check fastq file
def checkFastq(X):
    X = X.split(',')
    for q in X:
        if not os.path.exists(q):
            return False
        elif os.path.exists(q) and (os.path.getsize(q) == 0):
            return False
        else:
            pass
    return True  # return correct fastq flag


# download sra files and transfer the sra files to fastqfiles and cat the multiple fastq files
def catFastq(main_path, gsm, srr, lay_type):
    main_path = "%s/%s"%(main_path, gsm)
    cmd = 'mkdir %s \n'%main_path
    cmd = cmd + 'mkdir %s/fastq \n'%main_path
    cat_file1 = ''
    cat_file2 = ''
    srx_infor = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR'
    if lay_type == 'SINGLE':
        for i in range(len(srr)):
            srr[i] = srr[i][1:-12]
            ftp = srx_infor + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
            # download the sra files and transmit them into fastq files
            fsra = '%s/%s_%s.sra'%(main_path, gsm, i+1)
            cmd = cmd + 'wget %s -O %s \n'%(ftp, fsra) 
            cmd = cmd + '%s %s/%s_%s.sra -O %s/ \n'%(fdmp, main_path, gsm,i+1, main_path)
            cat_file1 = cat_file1 + '%s/%s_%s.fastq '%(main_path, gsm,i+1)
        cmd = cmd + 'cat %s> %s/fastq/%s.fastq \n'%(cat_file1,main_path,gsm)
        cmd = cmd + 'rm %s/%s_* \n'%(main_path, gsm)
        # cmd = cmd + 'rmdir %s/ \n'%(main_path)
    elif lay_type == 'PAIRED':
        for i in range(len(srr)):
            srr[i] = srr[i][1:-12]
            ftp = srx_infor + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
            fsra = '%s/%s_%s.sra'%(main_path, gsm, i+1)
            cmd = cmd + 'wget %s -O %s \n'%(ftp, fsra) 
            cmd = cmd + '%s --split-files %s/%s_%s.sra -O %s/\n'%(fdmp, main_path, gsm,i+1, main_path)
            cat_file1 = cat_file1 + '%s/%s_%s_1.fastq '%(main_path, gsm,i+1)
            cat_file2 = cat_file2 + '%s/%s_%s_2.fastq '%(main_path, gsm,i+1)
        cmd = cmd + 'cat %s> %s/fastq/%s.fastq_R1 \n'%(cat_file1, main_path,gsm)
        cmd = cmd + 'cat %s> %s/fastq/%s.fastq_R2 \n'%(cat_file2, main_path,gsm)
        cmd = cmd + 'rm %s/%s_* \n'%(main_path, gsm)
        # cmd = cmd + 'rmdir %s/ \n' % (main_path)
    return cmd


# get link sequece type and return the path of sbatch file
def gsmInfo(gsm, path):
    path = "%s/%s/"%(path, gsm)
    os.system('echo %s'%gsm)
        # gsm_url is the link of the input GSM data
    try:
        gsm_url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s'%gsm
        gsm_handler = urllib2.urlopen(gsm_url)
        gsm_html = gsm_handler.read()
        # get the ftp location of SRX file and the SRX id
        # ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX951/SRX951932
        srx_regexp = re.compile('https://www.ncbi.nlm.nih.gov/sra\S*"')
        srx_infor = srx_regexp.search(gsm_html)
        if srx_infor:
            srx = srx_infor.group().rstrip('"').lstrip('https://www.ncbi.nlm.nih.gov/sra?term=')
        else:
            os.system('echo "no srx of %s"'%gsm)
        # get the SRR id('>SRR1588518</a></td><td') and find the type of layout
        srx_url = 'http://www.ncbi.nlm.nih.gov/sra?term=%s'%srx
        srx_handler = urllib2.urlopen(srx_url)
        srx_html = srx_handler.read()
        # find the layout type (<div>Layout: <span>SINGLE</span>)
        lay_infor = re.compile('<div>Layout: <span>.{6}</span>')
        lay_type = lay_infor.search(srx_html)
        lay_type = lay_type.group()
        lay_type = lay_type[-13:-7]
        srr_regexp = re.compile('>SRR[0-9]*</a></td><td')
        srr = srr_regexp.findall(srx_html)
        if lay_type == 'SINGLE':
            fastq_file = '%s/fastq/%s.fastq'%(path, gsm)
        elif lay_type == 'PAIRED':
            fastq_file = '%s/fastq/%s.fastq_R1,%s/fastq/%s.fastq_R2'%(path, gsm, path, gsm)
        return srr, lay_type, fastq_file
    except:
        print 'failure gsm : %s' % gsm
        return None


# Generate the command to run salmon
def salmonRun(gsm, species, lay_type, output):
    ### fastq should be a list of fastq files
    cmd = ''
    output = "%s/%s"%(output, gsm)
    fastq = "%s/fastq/%s.fastq" % (output, gsm)
    fastq1 = "%s/fastq/%s.fastq_R1" % (output, gsm)
    fastq2 = "%s/fastq/%s.fastq_R2" % (output, gsm)
    if species == 'hg38':
        cmd = '%s quant -i %s -l -A'%(smn, hg38_index)
    elif species == 'mm10':
        cmd = '%s quant -i %s -l -A' % (smn, mm10_index)
    else:
        print "We do not support other genome index! Select hg38 or mm10 instead!"
    if lay_type == "SINGLE" and os.path.exists(fastq):
        cmd = cmd + ' -r %s -o %s \n'%(fastq, output)
    elif os.path.exists(fastq1) and os.path.exists(fastq2):
        cmd = cmd + ' -1 %s -2 %s -o %s \n'%(fastq1, fastq2, output)
    return cmd


def main():
    try:
        parser = argparse.ArgumentParser(description="""get fastq file from gsm and run salmon maybe""")
        parser.add_argument( '-s', '--species', dest='species', type=str, required=False, choices=['hg38', 'mm10'], help='select either hg38 or mm10')
        parser.add_argument( '-g', '--gsm', dest='gsm', type=str, required=True, help='gsm ID of the file')
        parser.add_argument( '-o', '--output', dest='outputdir', type=str, required=True, help='directory of output file')
        parser.add_argument( '-sal', '--salmon', dest='salmon', action='store_true', help='whether to run salmon for the fastq')
        args = parser.parse_args()
        if args.salmon:
            print "You will run salmon for %s"%args.gsm
            if not args.species:
                sys.stderr.write("Plese include species information if you want to run salmon!\n")
                sys.exit(0)
            else:
                srr, lay_type, fastq_file = gsmInfo(args.gsm, args.outputdir)
                if checkFastq(fastq_file):
                    print "Fastq file of %s exists, I will skip download of fastq file"%args.gsm
                    cmd = salmonRun(args.gsm, args.species, lay_type, args.outputdir)
                else:
                    cmd = catFastq(args.outputdir, args.gsm, srr, lay_type)
                    cmd = cmd + salmonRun(args.gsm, args.species, lay_type, args.outputdir)
        else:
            print "You will only get fastq for %s"%args.gsm
            srr, lay_type, fastq_file = gsmInfo(args.gsm, args.outputdir)
            if checkFastq(fastq_file):
                print "Fastq file of %s exists, I will skip download of fastq file"%args.gsm
                cmd = " "
            else:
                cmd = catFastq(args.outputdir, args.gsm, srr, lay_type)
        cmd = cmd + "rm %s/%s.sh \n"%(args.outputdir, args.gsm)
        cmd_file = open("%s/%s.sh"%(args.outputdir, args.gsm), "w")
        cmd_file.write(cmd)
        cmd_file.close()
        os.system("nohup bash %s/%s.sh > %s.log &"%(args.outputdir, args.gsm, args.gsm))
        # os.system("rm %s/%s.sh"%(args.outputdir, args.gsm))
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)

if __name__ == '__main__':
    main()