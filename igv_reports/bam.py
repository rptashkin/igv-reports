import pysam
from igv_reports.chralias import build_aliastable
import subprocess

class BamReader:

    def __init__(self, filetype, filename, args = None):
        self.filtetype = filetype
        self.filename = filename
        self.args = args

        samargs = "samtools view -H " + filename
        header = str(subprocess.check_output(samargs, shell=True))
        seqnames = parse_seqnames(header)
        self.aliastable = build_aliastable(seqnames)

    # add sam flag for unit tests
    def slice(self, region=None, region2=None,  sam=False):
        if sam:
            samargs = "samtools view -h " + self.filename
        else:
            samargs = "samtools view -b -h " +  self.filename
        if self.filtetype == "cram" and self.args is not None:
            samargs = samargs +  (" -T ")
            samargs = samargs + self.args.fasta

        if self.args is not None and hasattr(self.args, "exclude_flags"):
            if self.args.exclude_flags != 0:
                samargs = samargs + " -F "
                samargs = samargs + str(self.args.exclude_flags)
        else:
            samargs = samargs + "-F "
            samargs = samargs + "1536 "

        if region:
            range_string = self.get_chrname(region['chr']) + ":" + str(region['start']) + "-" + str(region['end'])
            samargs = samargs + " " + range_string
        if region2:
            range_string = self.get_chrname(region2['chr']) + ":" + str(region2['start']) + "-" + str(region2['end'])
            samargs = samargs + " " + range_string
        #print(samargs)
        return subprocess.check_output(samargs, shell=True)

    def get_chrname(self, c):
        if c in self.aliastable:
            return self.aliastable[c]
        else:
            return c



def parse_seqnames(header):
    seqnames = []
    lines = header.split('\n')
    for line in lines:
        if line.startswith('@SQ'):
            idx1 = line.find("SN:")
            if idx1 > 0:
                idx2 = line.find("\t", idx1)
                seqnames.append(line[idx1+3:idx2])

    return seqnames
