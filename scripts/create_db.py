#!/usr/bin/env python3
import argparse
import gffutils

parser = argparse.ArgumentParser()
parser.add_argument("annotation_file", help="GTF file containing gene annotations. For example, from https://www.gencodegenes.org/")
parser.add_argument("--filter", default="Ensembl_canonical", help="Only keep GTF features with the specified tags. Format: tag1,tag2,... or None to keep all features. Default: Ensembl_canonical")
args = parser.parse_args()

gtf = args.annotation_file
if gtf.endswith(".gtf"):
    prefix = gtf[:-4]
elif gtf.endswith(".gtf.gz"):
    prefix = gtf[:-7]
else:
    exit("ERROR, annotation_file should be a GTF file.")

def filter(feat):
    if feat.featuretype not in ["gene","transcript","exon"]:
        return False
    elif args.filter != "None" and feat.featuretype in ["transcript","exon"]:
        present = False
        for tag in args.filter.split(','):
            if "tag" in feat.attributes and tag in feat["tag"]:
                present = True
        if not present:
            return False
    return feat

db = gffutils.create_db(gtf, prefix+".db", force=True,
                        disable_infer_genes=True, disable_infer_transcripts=True,
                        transform=filter)

print("Database created: %s.db" % prefix)
