#!/usr/bin/env python3
"""Read project settings. This first script only prints a preview."""

# argparse reads command-line options and provides the --help page.
import argparse
# os provides file and folder checks that work on Linux and other systems.
import os


def build_parser():
    """Describe the arguments the user can enter at the command line."""
    # __doc__ is the description in the triple quotes at the top of this file.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', help='Protein family name')
    parser.add_argument('--taxon', help='Taxonomic group name or NCBI taxonomic ID')
    parser.add_argument('--fasta', help='Local unaligned protein FASTA file')
    parser.add_argument('--preset', choices=['aves'], help='Glucose-6-phosphatase in birds')
    parser.add_argument('--out', required=True, help='New output folder under an existing folder')
    # type=int converts input text to an integer and rejects values like 'ten'.
    # default supplies a value when the user leaves the option out.
    parser.add_argument('--max-download', type=int, default=1000, help='Maximum input records (default: 1000)')
    parser.add_argument('--max-align', type=int, default=100, help='Maximum alignment sequences (default: 100)')
    parser.add_argument('--seed', type=int, default=42, help='Sampling seed (default: 42)')
    parser.add_argument('--window', type=int, default=15, help='Conservation window size (default: 15)')
    # store_true makes a switch False by default and True when supplied.
    # For example, --skip-motifs requests a partial run without motif scanning.
    parser.add_argument('--skip-conservation', action='store_true')
    parser.add_argument('--skip-motifs', action='store_true', help='Partial run only')
    parser.add_argument('--skip-properties', action='store_true')
    return parser


def check_args(args):
    """Check and tidy settings; raise ValueError if something needs fixing.

    args is the collection of values returned by argparse. There is no custom
    settings class. A later GUI can supply the same values to this function.
    """
    # Expand the named example into a family and taxon. Reject conflicting
    # inputs rather than silently choosing which instruction to follow.
    if args.preset:
        if args.family is not None or args.taxon is not None or args.fasta is not None:
            raise ValueError('Use --preset alone, without --family, --taxon or --fasta.')
        args.family = 'glucose-6-phosphatase'
        args.taxon = 'Aves'
        # The preset has now been expanded into ordinary search settings.
        args.preset = None

    # None means the argument was not supplied. strip() removes surrounding
    # spaces, so a name such as " Aves " becomes "Aves".
    if args.family is not None:
        args.family = args.family.strip()
    if args.taxon is not None:
        args.taxon = args.taxon.strip()

    # The user can supply local sequences OR ask for a database search.
    # Mixing the two routes would make the intended input unclear.
    if args.fasta is not None:
        if args.family is not None or args.taxon is not None:
            raise ValueError('Choose a FASTA file OR a family and taxon search.')
        # expanduser replaces ~ with the home folder; abspath gives a full path.
        args.fasta = os.path.abspath(os.path.expanduser(args.fasta))
        if not os.path.isfile(args.fasta):
            raise ValueError('The FASTA file does not exist or is not a file.')
        # Actual FASTA contents are checked by quality.py.
    elif not args.family or not args.taxon:
        raise ValueError('Supply both --family and --taxon, or use --fasta.')

    # These are checks on settings, not the actual dataset. Later stages must
    # also check how many sequences were found and how long the alignment is.
    if args.max_download < 1 or args.window < 1:
        raise ValueError('Maximum input records and window size must be positive.')
    if args.max_align < 2:
        raise ValueError('Maximum alignment size must be at least 2.')
    if args.seed < 0:
        raise ValueError('The sampling seed cannot be negative.')
    if args.skip_conservation and args.skip_motifs and args.skip_properties:
        raise ValueError('Keep at least one analysis enabled.')

    # Protect previous results. lexists also detects a broken symbolic link,
    # which still occupies a name even though its target no longer exists.
    args.out = os.path.abspath(os.path.expanduser(args.out))
    if os.path.lexists(args.out):
        raise ValueError('Output path already exists. Choose a new folder name.')
    if not os.path.isdir(os.path.dirname(args.out)):
        raise ValueError('The parent folder for the output must already exist.')
    # Return the checked argument collection; values are accessed as args.taxon,
    # args.out, etc. This function has not created any folders or files.
    return args


def main():
    parser = build_parser()
    # parse_args reads what the user typed after the script name. Hyphens in
    # option names become underscores: --max-align becomes args.max_align.
    args = parser.parse_args()
    try:
        args = check_args(args)
    # ValueError covers invalid settings; OSError covers file-system problems.
    # parser.error shows a readable message and stops with a failure exit code.
    except (ValueError, OSError) as error:
        parser.error(str(error))

    print('Settings preview only: no downloads or analyses have run.\n')
    # vars() gives a dictionary of argument names and their values.
    for name, value in vars(args).items():
        print(name, '=', value)


# Run main only when this script is launched directly. Importing its functions
# from a future pipeline or GUI should not start the command-line interface.
if __name__ == '__main__':
    main()
