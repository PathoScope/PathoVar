import os
import sys
import argparse

import snp_utils

argparser = argparse.ArgumentParser(prog="pathovar.snp_caller.snp_utils")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("--test", action = "store_true", required = False, help="Enter IPython Interactive Session after execution completes [Development Only]")

sub_commands = argparser.add_subparsers(help="sub-command help")

filter_sam_file = sub_commands.add_parser("filter-sam", help="Filter a .sam file")
filter_sam_file_actions = filter_sam_file.add_mutually_exclusive_group(required = True)
filter_sam_file_actions.add_argument('--remove', action='store', type=str, help="A regular expression matching the reference to be removed")
filter_sam_file_actions.add_argument('--rename', action='store', nargs=2, type=str, help="A regular expression matching the reference to be renamed, followed by the new name")
filter_sam_file_actions.add_argument('-o', '--output', action='store', type=str, help="Name of file output file to save to")
filter_sam_file.add_argument("sam_file", action='store', help=".sam file to operate on")

filer_vcf_file = sub_commands.add_parser("filter-vcf", help="Filter a .vcf file")


def filter_sam_remove(args):
    parser = snp_utils.SAMParser(args.sam_file)
    parser.parse_file()
    target = args.remove.decode('string-escape')
    parser.remove(target)





def main(args):
    pass

if __name__ == '__main__':
    args = argparser.parse_args()
    if args.verbose: print(args)
    main(args)