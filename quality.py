#!/usr/bin/env python3
"""Read unaligned protein FASTA and produce a basic quality report.

Uses functions, lists, dictionaries and loops; no BioPython required.
This stage reports problems but does not remove or change biological records.
"""

import argparse
# statistics supplies median(), the middle value of the sequence lengths.
import statistics


def read_fasta(filename, max_records=1000):
    """Return a list of dictionaries containing id, header and sequence.

    Sequences can span several lines. Blank lines are ignored; letters are
    converted to uppercase. Duplicate identifiers and empty records are errors.
    """
    if max_records < 1:
        raise ValueError('The maximum record count must be positive.')
    # records holds every protein as a dictionary. seen_ids stores unique IDs
    # so we can spot duplicates. current refers to the record being assembled.
    records = []
    seen_ids = set()
    current = None

    # with closes the file automatically, including when an error occurs.
    # enumerate supplies line numbers starting at 1 for useful error messages.
    with open(filename, encoding='utf-8') as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            # A > line begins a new record, so first check that the preceding
            # record received a sequence before starting the next one.
            if line.startswith('>'):
                if current is not None and not current['sequence']:
                    raise ValueError('Empty sequence: ' + current['id'])
                # [1:] removes the leading > while retaining the description.
                header = line[1:].strip()
                if not header:
                    raise ValueError('Empty FASTA header on line ' + str(line_number))
                # FASTA identifiers are the first whitespace-separated word.
                # Keep the full header too, for descriptions and organism labels.
                identifier = header.split()[0]
                if identifier in seen_ids:
                    raise ValueError('Duplicate sequence ID: ' + identifier)
                if len(records) >= max_records:
                    raise ValueError('Too many records. Increase --max-records or use a smaller dataset.')
                seen_ids.add(identifier)
                current = {'id': identifier, 'header': header, 'sequence': ''}
                # The list and current refer to the SAME dictionary. Adding
                # sequence letters below therefore updates this list entry too.
                records.append(current)
            else:
                if current is None:
                    raise ValueError('Sequence found before the first >header.')
                # Join removes spaces/tabs used to format some FASTA files.
                # += appends each sequence line to the letters already read.
                current['sequence'] += ''.join(line.split()).upper()

    if not records:
        raise ValueError('No FASTA records were found.')
    # The final record has no following header to trigger the earlier check.
    # [-1] selects the last list entry, so check it separately at end of file.
    if not records[-1]['sequence']:
        raise ValueError('Empty sequence: ' + records[-1]['id'])
    return records


def check_sequences(records):
    """Check characters and return human-readable warnings.

    B, Z, J and X represent ambiguous/unknown residues. U and O are accepted
    special amino acids but may need extra care in downstream programs.
    A nucleotide-like alphabet is a warning, not proof: DNA letters are also
    amino-acid letters. Never translate or reject solely on this heuristic.
    """
    # Call this after read_fasta: it expects non-empty records and sequences.
    # Sets store distinct letters, letting us check the alphabet in one step.
    allowed = set('ACDEFGHIKLMNPQRSTVWYBXZJUO')
    warnings = []
    seen_sequences = {}
    lengths = []
    for record in records:
        sequence = record['sequence']
        # Set subtraction leaves only characters absent from the allowed set.
        invalid = set(sequence) - allowed
        if invalid:
            raise ValueError(record['id'] + ': invalid characters ' +
                             ', '.join(sorted(invalid)) +
                             '. Supply unaligned proteins without gaps or stop symbols.')
        lengths.append(len(sequence))

    # Finish collecting all lengths before comparing any individual record.
    # A second loop can now use the same dataset median for every protein.
    median_length = statistics.median(lengths)
    for record in records:
        sequence = record['sequence']
        identifier = record['id']
        # issubset asks whether EVERY distinct letter belongs to this alphabet.
        # This catches some DNA/RNA-like inputs, but cannot prove sequence type.
        if set(sequence).issubset(set('ACGTUN')):
            warnings.append(identifier + ': could be nucleotide data; confirm the source is protein.')
        # Reset the counter for each protein, then count ambiguous residues.
        ambiguous = 0
        for residue in sequence:
            if residue in 'BXZJ':
                ambiguous += 1
        if ambiguous:
            warnings.append(identifier + ': ' + str(ambiguous) + ' ambiguous/unknown residues.')
        if 'U' in sequence or 'O' in sequence:
            warnings.append(identifier + ': contains U/O; check downstream tool support.')
        if 'fragment' in record['header'].lower():
            warnings.append(identifier + ': header labels this as a fragment.')
        # 0.7 and 1.3 define a simple review threshold, not a deletion rule.
        if len(sequence) < 0.7 * median_length or len(sequence) > 1.3 * median_length:
            warnings.append(identifier + ': length differs by more than 30% from the median (review heuristic).')
        # This dictionary maps a complete sequence to its first record ID.
        # Identical proteins can occur in different organisms, so retain them.
        if sequence in seen_sequences:
            warnings.append(identifier + ': identical sequence to ' + seen_sequences[sequence] + '; retained.')
        else:
            seen_sequences[sequence] = identifier
    # A warning requests review; unlike raise ValueError, it does not stop the run.
    return warnings


def make_report(records, warnings):
    """Build one text report for display, saving and later GUI use."""
    lengths = []
    # A set counts each organism label once, even if several records share it.
    organisms = set()
    unknown = 0
    for record in records:
        lengths.append(len(record['sequence']))
        # Only recognise a trailing [organism] label. Do not guess from IDs.
        header = record['header']
        organism = ''
        if header.endswith(']') and '[' in header:
            # Split at the last [, keep its right-hand part, remove the final ].
            # Example: 'id protein [Bird A]' gives the organism label 'Bird A'.
            organism = header.rsplit('[', 1)[1][:-1].strip()
        if organism:
            organisms.add(organism)
        else:
            unknown += 1

    # Build the report as a list of lines first; join them into one string below.
    # str() converts numbers to text so they can be combined with labels.
    lines = ['PROTEIN INPUT QUALITY REPORT',
             'Records: ' + str(len(records)),
             'Distinct organism labels: ' + str(len(organisms)),
             'Records without an organism label: ' + str(unknown),
             'Length minimum / median / maximum: ' +
             str(min(lengths)) + ' / ' + str(statistics.median(lengths)) + ' / ' + str(max(lengths)),
             'Organism labels come from headers; they are not verified species counts.', '']
    if len(records) == 1:
        lines.append('One sequence: skip conservation; motif/property analysis can still run.')
    if len(organisms) < 2:
        lines.append('Multiple species are not established. Review metadata before continuing.')
    lines.append('Warnings: ' + str(len(warnings)))
    for warning in warnings:
        lines.append('- ' + warning)
    lines.append('\nNo records were removed. Review this report before analysis.')
    # \n is a newline. The same returned text can be printed or saved.
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True, help='Unaligned protein FASTA file')
    parser.add_argument('--max-records', type=int, default=1000)
    parser.add_argument('--report', help='Optional new text report file; existing files are protected')
    args = parser.parse_args()
    try:
        # The output of each function becomes the input to the next step:
        # FASTA file -> record dictionaries -> warnings -> readable report.
        records = read_fasta(args.fasta, args.max_records)
        warnings = check_sequences(records)
        report = make_report(records, warnings)
        if args.report:
            # Mode x creates a new file and refuses to overwrite an old report.
            with open(args.report, 'x', encoding='utf-8') as handle:
                handle.write(report)
        print(report)
    # Catch malformed input, file-access errors and text-decoding errors.
    # Report the problem to the user rather than showing a Python traceback.
    except (ValueError, OSError, UnicodeError) as error:
        parser.error(str(error))


# Importing these functions into another script must not start a separate run.
if __name__ == '__main__':
    main()
