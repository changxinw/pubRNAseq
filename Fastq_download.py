import os, sys
import urllib2
import re
# import commands
import argparse

#Set up the path to the fastq-dump here
fdmp = 'fastq-dump'
###This code still require the checking read length script

# transfer the sra files to fastqfiles and cat the multiple fastq files
def catfastq(main_path, gsm, srr, lay_type):
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
        cmd = cmd + 'rmdir %s/ \n' % (main_path)
    return cmd

# get link sequece type and return the path of sbatch file
def LinkPlusDownload(gsm, path, refresh = True):
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
            fastq_file = '%s/fastq/%s.fastq_R1,%s/fastq/%s.fastq_R2'%(path,gsm,path,gsm)
        cmd_file = open("%s/%s.sh"%(path, gsm),"w")
        def checkFastq(X):
            X = X.split(',')
            for q in X:
                if not os.path.exists(q):
                    return False
                elif os.path.exists(q) and (os.path.getsize(q) == 0):
                    return False
                else:
                    pass
            return True # return correct fastq flag
        # check whether need to download fastq again, by fastq existence and the size ?= 0
        if refresh: # new samples need to download fastq
            cmd_file.write(catfastq(path, gsm, srr, lay_type)+'\n')
        elif (not refresh) and (not checkFastq(fastq_file)):# rerun samples but fastq is not correct.
            cmd_file.write(catfastq(path, gsm, srr, lay_type)+'\n')
        else:
            pass
        cmd_file.close()
        return '%s/%s.sh'%(path,gsm),lay_type
    except:
        print(sys.exc_info())
        return None, None
        # print 'failure gsm : %s'%gsm


def main():
    try:
        parser = argparse.ArgumentParser(description="""download srr file and translate to fastq""")
        # parser.add_argument( '-i', dest='inputable', type=str, required=True, help='path of the input table')
        parser.add_argument( '-g', dest='gsm', type=str, required=True, help='gsm ID of the file')
        parser.add_argument( '-o', dest='outputdir', type=str, required=True, help='directory of output file')
        args = parser.parse_args()
        LinkPlusDownload(args.gsm, args.outputdir)
        # os.system("bash %s/%s.sh"%(args.outputdir, args.gsm))
        # os.system("rm %s/%s.sh"%(args.outputdir, args.gsm))
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)

if __name__ == '__main__':
    main()
