#!/usr/bin/env python3
"""Download candidate proteins with EDirect, then run the existing QC functions.

Example: python3 retrieval.py --preset aves --out aves_run
Requires esearch, efetch and xtract on PATH. No BioPython or shell pipelines are used.
"""

import csv
import json
import os
import shutil
import subprocess

from config import build_parser, check_args
from quality import read_fasta, check_sequences, make_report


def run_edirect(command, input_text=None):
    """Run one EDirect command and return its text output.

    command is a list, e.g. ['esearch', '-db', 'taxonomy', '-query', 'Aves'].
    input_text passes a previous search result to efetch, replacing a shell pipe.
    """
    try:
        # No shell=True: user input is passed as an argument, not executed as code.
        # capture_output keeps stdout and stderr separate. A timeout prevents an
        # unavailable service from leaving this script waiting indefinitely.
        result = subprocess.run(command, input=input_text, text=True,
                                capture_output=True, timeout=300, check=True)
    except FileNotFoundError:
        raise ValueError(command[0] + ' was not found. Load EDirect on the MSc server.')
    except subprocess.TimeoutExpired:
        raise ValueError(command[0] + ' timed out after 5 minutes. Try again later.')
    except subprocess.CalledProcessError as error:
        raise ValueError(command[0] + ' failed: ' + (error.stderr or 'No error details returned.').strip())
    # Treat extraction diagnostics as failures rather than accepting partial rows.
    if command[0] == 'xtract' and result.stderr.strip():
        raise ValueError('xtract reported a problem: ' + result.stderr.strip())
    if result.stderr.strip():
        print('EDirect message:', result.stderr.strip())
    if not result.stdout.strip():
        raise ValueError(command[0] + ' returned no output. Check the query or connection.')
    return result.stdout


def check_xml_errors(xml_text, record_name):
    """Use xtract to check the expected response container for service errors.

    A normal response produces three NONE placeholders. A missing container,
    reported error, or unexpected output stops the run. This is not a separate
    full XML validator: extraction diagnostics and later record checks also matter.
    """
    command = ['xtract', '-pattern', record_name, '-def', 'NONE',
               '-element', 'ERROR', 'Error', 'ErrorList']
    result = run_edirect(command, xml_text).strip()
    if result != 'NONE\tNONE\tNONE':
        raise ValueError('Unexpected NCBI response or service error: ' + result)


def search_count(search_xml):
    """Ask xtract for the match count, then check that it is one whole number."""
    check_xml_errors(search_xml, 'ENTREZ_DIRECT')
    command = ['xtract', '-pattern', 'ENTREZ_DIRECT', '-element', 'Count']
    count = run_edirect(command, search_xml).strip()
    if not count.isdigit():
        raise ValueError('The search response did not contain a valid record count.')
    return int(count)


def resolve_taxon(taxon):
    """Resolve exactly one taxonomic ID, so a broad/ambiguous name is not guessed."""
    # Numeric IDs are still searched to confirm that the taxon exists.
    if taxon.isdigit():
        query = taxon + '[UID]'
    else:
        if '"' in taxon or '[' in taxon or ']' in taxon:
            raise ValueError('Enter a taxon name or numeric ID, not a taxonomy query.')
        query = '"' + taxon + '"[All Names]'
    result = run_edirect(['esearch', '-db', 'taxonomy', '-query', query])
    if search_count(result) != 1:
        raise ValueError('Taxon must match exactly one record. Try its NCBI taxonomic ID.')
    identifier = run_edirect(['efetch', '-format', 'uid'], result).strip()
    if not identifier.isdigit():
        raise ValueError('NCBI did not return one valid taxonomic ID.')
    return identifier


def extract_proteins(xml_text):
    """Ask xtract for one tab-separated row per protein record.

    -pattern selects the repeating record; -element lists its required fields.
    -def supplies a missing-value marker so absent fields do not shift columns.
    xtract also checks the response for service errors before extracting proteins.
    """
    check_xml_errors(xml_text, 'GBSet')
    command = [
        'xtract', '-pattern', 'GBSeq', '-def', 'MISSING',
        '-element', 'GBSeq_accession-version', 'GBSeq_organism',
        'GBSeq_definition', 'GBSeq_length', 'GBSeq_sequence',
    ]
    return run_edirect(command, xml_text)


def parse_proteins(table_text, expected_count):
    """Check xtract's five columns and return ordinary protein dictionaries.

    Each row keeps an accession, its metadata and its sequence together. No
    separate downloads are joined by row position. Missing data stops the run.
    """
    proteins = []
    seen = set()
    for line_number, line in enumerate(table_text.splitlines(), start=1):
        # Split specifically on tabs, not whitespace: descriptions contain spaces.
        fields = line.split('\t')
        if len(fields) != 5:
            raise ValueError('xtract row ' + str(line_number) + ' must contain five columns.')
        accession, organism, description, length, sequence = fields
        accession = accession.strip()
        organism = organism.strip()
        description = description.strip()
        length = length.strip()
        sequence = ''.join(sequence.split()).upper()
        for field in [accession, organism, description, length, sequence]:
            if not field or field == 'MISSING':
                raise ValueError('xtract row ' + str(line_number) + ' has a missing required field.')
        if accession in seen:
            raise ValueError('Duplicate accession in the download: ' + accession)
        if not length.isdigit() or int(length) != len(sequence):
            raise ValueError('Downloaded length does not match the sequence: ' + accession)
        seen.add(accession)
        proteins.append({'accession': accession, 'organism': organism,
                         'description': description, 'length': len(sequence),
                         'sequence': sequence})
    # Never silently continue with a partial download.
    if len(proteins) != expected_count:
        raise ValueError('Expected ' + str(expected_count) + ' proteins but received ' +
                         str(len(proteins)) + '. The download may be incomplete.')
    return proteins


def save_proteins(proteins, output_dir):
    """Save FASTA and a tab-separated metadata table; return the FASTA path."""
    fasta_path = os.path.join(output_dir, 'proteins.fasta')
    with open(fasta_path, 'x', encoding='utf-8') as handle:
        for protein in proteins:
            # A trailing organism label is understood by our quality.py reader.
            handle.write('>' + protein['accession'] + ' ' + protein['description'] +
                         ' [' + protein['organism'] + ']\n')
            sequence = protein['sequence']
            # Wrap at 60 letters for human readability; this changes no residues.
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + '\n')
    with open(os.path.join(output_dir, 'metadata.tsv'), 'x', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['accession', 'organism', 'description', 'length'])
        for protein in proteins:
            writer.writerow([protein['accession'], protein['organism'],
                             protein['description'], protein['length']])
    return fasta_path


def retrieve(args):
    """Search, download and check candidates. Return the quality-checked records.

    args comes from config.check_args(). This function handles search input only;
    a later orchestrator will also offer the existing local-FASTA route.
    """
    if args.fasta is not None:
        raise ValueError('For local input, run quality.py --fasta YOUR_FILE instead.')
    for program in ['esearch', 'efetch', 'xtract']:
        if shutil.which(program) is None:
            raise ValueError(program + ' was not found. Load EDirect on the MSc server.')
    # Create a new folder exclusively. If it already exists, stop without changing it.
    os.mkdir(args.out)
    try:
        print('Resolving taxonomic group:', args.taxon)
        taxon_id = resolve_taxon(args.taxon)
        # The family field is an NCBI protein search expression. Parentheses keep
        # OR alternatives inside the family clause and preserve the taxon filter.
        query = '(' + args.family + ') AND txid' + taxon_id + '[Organism:exp]'
        print('Searching NCBI Protein:', query)
        search = run_edirect(['esearch', '-db', 'protein', '-query', query])
        with open(os.path.join(args.out, 'search.xml'), 'x', encoding='utf-8') as handle:
            handle.write(search)
        count = search_count(search)
        # Retain settings and query even if the count is zero or too large.
        with open(os.path.join(args.out, 'retrieval.json'), 'x', encoding='utf-8') as handle:
            json.dump({'query': query, 'taxon_id': taxon_id, 'count': count,
                       'settings': vars(args)}, handle, indent=2)
        print('Matching records:', count)
        if count == 0:
            raise ValueError('No matches. Try another family name or check the taxon.')
        if count > args.max_download:
            raise ValueError('Too many matches. Narrow the family/taxon or explicitly increase --max-download.')

        print('Downloading protein records and metadata...')
        xml_text = run_edirect(['efetch', '-format', 'gp', '-mode', 'xml'], search)
        with open(os.path.join(args.out, 'proteins.xml'), 'x', encoding='utf-8') as handle:
            handle.write(xml_text)
        print('Extracting protein fields with xtract...')
        table_text = extract_proteins(xml_text)
        # Preserve the raw extraction for troubleshooting; it has no header row.
        with open(os.path.join(args.out, 'extracted_proteins.tsv'), 'x', encoding='utf-8') as handle:
            handle.write(table_text)
        proteins = parse_proteins(table_text, count)
        fasta_path = save_proteins(proteins, args.out)
        records = read_fasta(fasta_path, args.max_download)
        warnings = check_sequences(records)
        report = make_report(records, warnings)
        with open(os.path.join(args.out, 'quality_report.txt'), 'x', encoding='utf-8') as handle:
            handle.write(report)
        print(report)
        print('Download and QC complete. Review candidates before analysis; no analyses have run.')
        return records
    except (ValueError, OSError, UnicodeError) as error:
        # Preserve partial files for diagnosis instead of silently deleting them.
        # Even a failure to write this note must not hide the original error.
        try:
            with open(os.path.join(args.out, 'ERROR.txt'), 'x', encoding='utf-8') as handle:
                handle.write(str(error) + '\nThis run is incomplete. Use a new folder when retrying.\n')
        except OSError:
            pass
        raise


def main():
    # Reuse the same argument names/defaults as config.py, not a second copy.
    parser = build_parser()
    parser.description = __doc__
    args = parser.parse_args()
    try:
        args = check_args(args)
        retrieve(args)
    except (ValueError, OSError, UnicodeError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
