
"""
this script puts together the config file (run.json) for the RRBS or C&R pipelines

make_config.py \
  --prefix prefix \
  --template Templates/template_RRBS.json \
  --meta meta.tab \
  --exp_type RRBS 

make_config.py \
  --prefix prefix \
  --template Templates/template_CutAndRun.json \
  --meta pairs.tab \
  --exp_type CutAndRun

"""
import glob
import argparse
import os
import sys
import re
from csv import reader
import json


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        required=True,
                        help="Experiment name")
    ## Mandatory options
    parser.add_argument("--template",
                        required=True,
                        help="template json file")
    parser.add_argument("--meta",
                        required=True,
                        help="Tabfile holding samples information")
    parser.add_argument("--exp_type",
                        required=True,
                        help="Experiment type")
## other info:
    parser.add_argument("--annotation",
                        default="mm10",
                        help="Genome version")
    parser.add_argument("--gsize",
                        default="mm",
                        help="Genome version")
    parser.add_argument("--workpath",
                        default=".",
                        help="Workpath")
    parser.add_argument("--clust_config",
                        default="cluster.json",
                        help="cluster settings")
    parser.add_argument("--username",
                        default="tillodc@nih.gov",
                        help="email")
    args = parser.parse_args()
    return vars(args)



def extract_template(template):
    with open(template) as f:
        template_data = json.load(f)
    return template_data


def parse_samples_RRBS(pairs):
    samp = []
    samp_vec = []
    my_file=open(pairs)
    read_file=reader(my_file,delimiter="\t")
    samp_list=list(read_file)
    for row in samp_list[1:]:
        samp.append(row[0])
        samp_vec.append(int(row[1]))
    return samp,samp_vec

def parse_samples_CutAndRun(pairs):
    chips = []
    ctrls = []
    my_file=open(pairs)
    read_file=reader(my_file,delimiter="\t")
    samp_list=list(read_file)
    for row in samp_list[1:]:
        chips.append(row[0])
        ctrls.append(row[1])
    return chips,ctrls




if __name__ == "__main__":
    args = make_parser()
    project={}
    project["annotation"]=args["annotation"]
    project["cluster"]=args["clust_config"]
    project["experiment_name"]=args["prefix"]
    project["gsize"]=args["gsize"]
    project["data_type"]=args["exp_type"]
    project["username"]=args["username"]
    project["workpath"]=os.path.abspath(args["workpath"])
    if args["exp_type"]=="RRBS":
        samp,samp_vec = parse_samples_RRBS(args["meta"])
        project["samples"]=samp
        project["sample_vector"]=samp_vec
    elif args["exp_type"]=="CutAndRun":
        chips, ctrls = parse_samples_CutAndRun(args["meta"])
        project["chips"]=chips
        project["controls"]=ctrls
    run=extract_template(args["template"])
    run["project"]=project
    fout=open("run.json","w")
    fout.write(json.dumps(run,indent = 4,sort_keys=True))
    fout.close()

