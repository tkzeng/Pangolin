import argparse
import gffutils

parser = argparse.ArgumentParser()
parser.add_argument("annotation_file", help="GTF file containing gene annotations. For example, from https://www.gencodegenes.org/")
args = parser.parse_args()

gtf = args.annotation_file
if not gtf.endswith(".gtf"):
    exit("ERROR, annotation_file should be a GTF file.")

gffutils.create_db(gtf, gtf[:-4]+".db", force=True, verbose=True,
                   disable_infer_genes=True, disable_infer_transcripts=True)
print("Database created: %s.db" % gtf[:-4])