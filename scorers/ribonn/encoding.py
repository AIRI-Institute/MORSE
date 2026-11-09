"""Strict full-transcript routing and upstream RiboNN nucleotide encoding."""
import numpy as np

MAX_UTR5 = 1381
MAX_CDS_UTR3 = 11937
ALLOWED = ('five_prime_utr', 'start_codon', 'cds', 'stop_codon', 'three_prime_utr')


def transcript_parts(sequence, context):
    sequence = sequence.upper().replace('U', 'T')
    if not sequence or set(sequence) - set('ACGT'):
        raise ValueError('RiboNN requires a nonempty A/C/G/T/U transcript')
    ranges = context.get('part_ranges', [])
    cursor, previous, found = 0, -1, {}
    for r in ranges:
        name = r['part']
        if name not in ALLOWED or name in found or ALLOWED.index(name) <= previous:
            raise ValueError('Expected complete transcript-body parts in biological order, without poly(A)')
        start, end = int(r['start']), int(r['end'])
        if start != cursor or end <= start or end > len(sequence):
            raise ValueError('Part ranges must cover the routed sequence exactly')
        if r.get('part_start', 0) != 0:
            raise ValueError('Sliced regions are not supported')
        found[name] = (start, end)
        cursor, previous = end, ALLOWED.index(name)
    if cursor != len(sequence) or not all(k in found for k in ('start_codon', 'cds', 'stop_codon')):
        raise ValueError('Explicit start, CDS and terminal stop context is required')
    start = found['start_codon'][0]
    stop = found['stop_codon'][1]
    if found['start_codon'][1] - start != 3 or stop - found['stop_codon'][0] != 3:
        raise ValueError('Start and stop must each contain three nucleotides')
    if sequence[start:start + 3] != 'ATG' or sequence[stop - 3:stop] not in ('TAA', 'TAG', 'TGA'):
        raise ValueError('RiboNN scorer requires AUG and a canonical terminal stop')
    if (stop - start) % 3:
        raise ValueError('Complete CDS length must be divisible by three')
    if start > MAX_UTR5 or len(sequence) - start > MAX_CDS_UTR3:
        raise ValueError('Transcript exceeds upstream limits; truncation is forbidden')
    return sequence, start, stop - start


def encode(sequences, contexts):
    if len(sequences) != len(contexts):
        raise ValueError('One context is required per sequence')
    result = np.zeros((len(sequences), 5, MAX_UTR5 + MAX_CDS_UTR3), dtype=np.float32)
    mapping = {c: i for i, c in enumerate('ATCG')}
    for i, (seq, ctx) in enumerate(zip(sequences, contexts)):
        seq, utr5, cds_size = transcript_parts(seq, ctx)
        offset = MAX_UTR5 - utr5
        result[i, [mapping[c] for c in seq], np.arange(len(seq)) + offset] = 1
        # Upstream includes the start position of the terminal stop codon.
        result[i, 4, MAX_UTR5:MAX_UTR5 + cds_size:3] = 1
    return result

