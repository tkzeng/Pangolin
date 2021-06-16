import argparse
import gffutils

parser = argparse.ArgumentParser()
parser.add_argument("annotation_file", help="GTF file containing gene annotations. For example, from https://www.gencodegenes.org/")
args = parser.parse_args()

gtf = args.annotation_file
if gtf.endswith(".gtf"):
    prefix = gtf[:-4]
elif gtf.endswith(".gtf.gz"):
    prefix = gtf[:-7]
else:
    exit("ERROR, annotation_file should be a GTF file.")

gffutils.create_db(gtf, prefix+".db", force=True,
                   disable_infer_genes=True, disable_infer_transcripts=True)
print("Database created: %s.db" % prefix)
