#!/usr/bin/python

# Author: John Hawkins (jsh)

from collections import defaultdict
import itertools
import logging
import math
import multiprocessing
import numpy
import optparse
import os
import string
import sys

import HTSeq

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


# Yes, GLOBAL VARIABLES are horrible, but so is passing flags through seven
# levels of hell.  I'm compromising on these.
MINIMUM_BINOMIAL_COVERAGE = 8
NON_TERMINAL_LENGTH_CUTOFF = 25


def main():
  logging.info('Parsing command line.')
  usage = '%prog [options] input_shard [...] counts_file output_file_base'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option(
      '--max_read_len', type='int',
      action='store', dest='max_read_len',
      default=45,
      help='Max len of input reads.  Too big: Wasted RAM.  Too small: Crash.')
  parser.add_option(
      '--normalization_count', type='int',
      action='store', dest='normalization_count',
      default=10 * 1000 * 1000,
      help='Scale per-location counts by normalization_count/total_aligned.')
  parser.add_option(
      '--parallelism', type='int',
      action='store', dest='parallelism',
      default=6,
      help='More workers use more CPU/RAM, but do more work.')
  parser.add_option(
      '--desired_confidence', type='int',
      action='store', dest='desired_confidence',
      default=0.95,
      help='Report bottom of interval with specified confidence.')
  (options, args) = parser.parse_args()
  input_shards = args[:-2]
  count_file = args[-2]
  output_file_base = args[-1]
  if not input_shards:
    binary = os.path.basename(sys.argv[0])
    logging.fatal(
        '{0} called with no input files specified.'.format(binary))
    return 2
  if not os.path.exists(count_file):
    logging.fatal('No such file: {0}'.format(count_file))
    return 2
  try:
    total_count = int(open(count_file).read().strip())
  except ValueError, e:
    logging.fatal('Count file {0} could not be read as an integer.')
    raise
  for input_shard in input_shards:
    if not os.path.exists(input_shard):
      logging.fatal('No such file: {0}'.format(input_shard))
      return 2
    if not '.shard.' in input_shard:
      complaint = 'Shard file names must contain ".shard.".  {0} does not.'
      logging.fatal(complaint.format(input_shard))
      return 2
  process_input_shards(input_shards,
                       total_count,
                       options.parallelism,
                       options.normalization_count,
                       options.max_read_len,
                       output_file_base)


def chrom_dir_pos_dict(innertype):
  """Returns a defaultdict of defaultdicts of defaultdicts of innertypes.

  So you can do:
  foo = chrom_dir_pos_dict(int)
  foo['totally']['new']['key'] += 3

  Though of course the intent is more like:
  foo['chrXII']['+'][203448] += 1
  """
  return defaultdict(lambda : defaultdict(lambda : defaultdict(innertype)))


class PositionalStats(object):
  def __init__(self):
    self.terminal_coverage = 0
    self.terminal_mismatch = 0
    self.non_terminal_coverage = 0
    self.non_terminal_mismatch = 0

  def __add__(self, other):
    new = PositionalStats()
    new.terminal_coverage = self.terminal_coverage + other.terminal_coverage
    new.terminal_mismatch = self.terminal_mismatch + other.terminal_mismatch
    new.non_terminal_coverage = (self.non_terminal_coverage +
                                 other.non_terminal_coverage)
    new.non_terminal_mismatch = (self.non_terminal_mismatch +
                                 other.non_terminal_mismatch)
    return new

  def __iadd__(self, other):
    self.terminal_coverage += other.terminal_coverage
    self.terminal_mismatch += other.terminal_mismatch
    self.non_terminal_coverage += other.non_terminal_coverage
    self.non_terminal_mismatch += other.non_terminal_mismatch
    return self


def process_input_shards(input_shards,
                         total_count,
                         parallelism,
                         normalization_count,
                         max_read_len,
                         output_file_base):
  """Generate output shards for each input chromosome shard."""
  results = []
  pool = multiprocessing.Pool(processes=parallelism)
  for shard in input_shards:
    results.append(pool.apply_async(process_input_shard,
                                    [shard,
                                     total_count,
                                     normalization_count,
                                     max_read_len,
                                     output_file_base]))
  pool.close()
  pool.join()
  for r in results:
    r.get()


def process_input_shard(input_shard_name,
                        total_count,
                        normalization_count,
                        max_read_len,
                        output_file_base):
  """Generate output shards for one chromosome shard."""
  input_shard_base, shard_name = input_shard_name.split('.shard.')
  n_seen = 0
  counts_array = chrom_dir_pos_dict(PositionalStats)
  logging.info('Processing file: {0}'.format(input_shard_name))
  sam_alignments = HTSeq.SAM_Reader(input_shard_name)
  groups = itertools.groupby(sam_alignments, lambda x: x.read.name)
  for name, grouper in groups:
    n_seen += 1
    if n_seen % 100000 == 0:
      logging.info(
          'Now analyzing {shard_name} group #{n_seen} {name}.'.format(**vars()))
    alignments = list(grouper)
    features = features_for_alignment_group(alignments, max_read_len)
    if not features:
      continue
    for pos in features:
      c, d, p = pos
      if counts_array[c][d].has_key(p):
        counts_array[c][d][p] += features[pos]
      else:
        counts_array[c][d][p] = features[pos]
  write_track_files(counts_array,
                    float(normalization_count)/total_count,
                    shard_name,
                    output_file_base)


def features_for_alignment_group(alignments, max_read_len):
  """Return features of unique alignment, or None.

  None is returned if the position cannot be uniquely identified.
  """
  alignments = [x for x in alignments if x.aligned]
  alignments = [x for x in alignments if not '_' in x.iv.chrom]
  if len(alignments) == 1:
    return features_for_alignment(alignments[0], max_read_len)
  else:
    return None


def features_for_alignment(alignment, max_read_len):
  """Returns {genomic position: feature counts} mapping.

  Returns:
    features: {position->PositionalStats}
  """
  # TODO(jsh): This code needs to account for the cigar gap.
  # TODO(jsh): HTSeq returns the WHOLE alignment interval, so a cigar string
  # TODO(jsh):   like 10M2000N10M gives an interval of width 2020, which is
  # TODO(jsh):   extremely slow and generates incorrect statistics.  This is
  # TODO(jsh):   primarily relevant for tophat alignment.
  features = dict()
  (iv, read, md_string) = correctly_complement_alignment(alignment)
  positions = [(p.chrom, p.strand, p.pos)
                  for p in debugged_genomic_interval_xrange_d(iv)]
  sequence = read.seq
  scores = read.qual
  # We care about offset from terminal location, hence reverse-enumeration.
  for (i, p) in enumerate(reversed(positions)):
    features[p] = PositionalStats()
    if i == 0:
      features[p].terminal_coverage += 1
    elif i < NON_TERMINAL_LENGTH_CUTOFF:
      features[p].non_terminal_coverage += 1
    else:
      continue  # We ignore trailing sequence elements.
  for (idx, seq, ref, prob) in scored_mismatches(md_string,
                                                 sequence,
                                                 scores):
    # Ignore "mismatches" which are really just read-failures.
    if seq == 'N':
      continue
    # (We ignore prob because it's almost always *very* nearly 1.0.)
    # First compute offset for idx.
    offset_index = (len(sequence) - 1) - idx
    if offset_index == 0:
      features[positions[idx]].terminal_mismatch += 1
    elif offset_index < NON_TERMINAL_LENGTH_CUTOFF:
      features[positions[idx]].non_terminal_mismatch += 1
    else:
      continue  # We ignore trailing sequence elements.
  return features


def debugged_genomic_interval_xrange_d(gi, step=1):
  base_range = xrange(gi.start, gi.end, step)
  if gi.strand == '-':
    base_range = reversed(base_range)
  for pos in base_range:
    yield HTSeq.GenomicPosition(gi.chrom, pos, gi.strand)


def tokens_from_md(md_string):
  tokens = list()
  while md_string:
    if md_string[0].isdigit():
      token_end = find_if_not(str.isdigit, md_string)
      token, md_string = md_string[:token_end], md_string[token_end:]
    elif md_string[0] == '^':
      md_string = md_string[1:]
      if not md_string[0] in 'GCTAN':
        raise ValueError, '"^" followed by nonsense: ^{0}'.format(md_string)
      token_end = find_if_not(lambda x: x in 'GCTAN', md_string)
      deleted, md_string = md_string[:token_end], md_string[token_end:]
      token = '^' + deleted
    elif md_string[0] in 'GCTAN':
      token, md_string = md_string[:1], md_string[1:]
    else:
      raise ValueError, 'Unrecognized leading token: {0}'.format(md_string)
    tokens.append(token)
  return tokens


def revcomp_md(md_string):
  """Return reverse complement of md_string field."""
  def fliphats(token):
    if token.startswith('^'):
      return '^' + ''.join(reversed(token[1:]))
    else:
      return token
  md_string = md_string.translate(string.maketrans('GCAT', 'CGTA'))
  tokens = tokens_from_md(md_string)
  return ''.join(reversed(map(fliphats, tokens_from_md(md_string))))


def correctly_complement_alignment(alignment):
  """Fix HTSeq's failure to revcomp MD fields, then return revcomp(alignment)

  We do this in one step instead of two, because it can save us one reverse
  complement operation on the MD strings.

  Alignments have two oddities.  First, the MD string is always relative to
  the positive strand.  (This is actually true of the sequence, too, but that
  is correctly complemented where necessary by HTSeq.)  Second, our sequences
  are RNA, but what we're really trying to align is the DNA from which they
  were derived.  So we take the reverse complement of the MD string IFF the
  alignment was to the minus side, and then we reverse complement the whole
  alignment.

  Returns:
    (iv, read, md_string) -- HTSeq.SAM_Alignment isn't very mutable, sadly.
  """
  if alignment.iv.strand == '+':
    # revcomp(x) one time for RNA -> DNA mapping.
    md_string = revcomp_md(alignment.optional_field('MD'))
  elif alignment.iv.strand == '-':
    # revcomp(x) for mapping, and again to fix the inconsistency in SAM format.
    # revcomp(revcomp(x)) == x, so just use the original.
    md_string = alignment.optional_field('MD')
  else:
    raise ValueError, 'Alignment had strand "{0}".'.format(alignment.iv.strand)
  complement = alignment.read.seq.translate(string.maketrans('GCAT', 'CGTA'))
  read = alignment.read.get_reverse_complement()
  read.name = alignment.read.name
  flip_map = {'+': '-', '-': '+'}
  iv = alignment.iv
  iv.strand = flip_map[alignment.iv.strand]
  return (iv, read, md_string)


def wilson_95_from_tally(k, n):
  """Simple wrapper to compute Wilson interval based on raw counts."""
  return compute_wilson_interval_95(float(k)/n, n)


def compute_wilson_interval_95(p_est, n):
  """Returns wilson interval for p_true value at confidence z = 0.95.

  Args:
    p_est: Fraction of seen trials which were positive.
    n: Number of seen trials.
  Returns:
    (p_low, p_high): With confidence z, p_true is between p_low and p_high.
  """
  p_est = float(p_est)
  n = float(n)
  z = 1.96  # The Z-value for the 95% confidence interval.
  base = p_est + (1/(2*n)) * z*z
  radius = z*math.sqrt(p_est*(1-p_est)/n + z*z/(4*n*n))
  denom = 1 + (1/n)*z*z
  p_low = (base - radius)/denom
  p_high = (base + radius)/denom
  return (p_low, p_high)


def find_if_not(tagger, things):
  """Return the index of the first thing NOT tagged, or len(things)."""
  return next((i for i in range(len(things)) if not tagger(things[i])),
              len(things))


def scored_mismatches(md_string, sequence, scores):
  """Generator yielding alignment mismatches in sequence.

  Mismatches are returned as (pos, seq, ref, prob) where:

  idx: Sequence-relative zero-based index of mismatch.
  seq: Sequence letter at pos.
  ref: Reference letter aligned to pos.
  prob: Probability that mismatch is actually a mismatch.
  """
  idx = 0
  for (seq, ref) in aligned_fragments(md_string, sequence):
    if len(seq) != 1 or len(ref) != 1 or seq == ref:
      idx += len(seq)
    else:
      log = scores[idx]/-10
      prob = 1.0 - math.pow(10, log)
      yield (idx, seq, ref, prob)
      idx += len(seq)  # Yeah, it'll always be one, but hey.


def aligned_fragments(md_string, sequence):
  """Generator yielding matched fragments of sequence and reference.

  Parses md_string, which has the format specified in:

      http://samtools.sourceforge.net/SAM1.pdf  (page 6, footnote 2)

  Yields: (seq_fragment, ref_fragment) where
      seq_fragment: The next chunk of the provided sequence.
      ref_fragment: The corresponding chunk of the implied reference.
  """
  for token in tokens_from_md(md_string):
    if token[0].isdigit():
      match_len = int(token)
      fragment, sequence = sequence[:match_len], sequence[match_len:]
      yield (fragment, fragment)
    elif token[0] == '^':
      yield ('', token[1:])
    elif token in 'GCTAN':
      seq_base, sequence = sequence[0], sequence[1:]
      if seq_base == token:
        # If we ever see a mismatch between, say, 'A' and 'A', something broke.
        complaint = '{0} implies a non-mismatched mismatch for {1}.'
        raise ValueError, complaint.format(md_string, sequence)
      yield (seq_base, token)
    else:
      raise ValueError, 'Unrecognized leading token: {0}'.format(md_string)
  if sequence:
    raise ValueError, 'Sequence continues past md_string: {0} / {1}'.format(
        sequence, md_string)


def terminal_coverage(features):
  return features.terminal_coverage

def terminal_mismatch(features):
  return features.terminal_mismatch

def terminal_mismatch_wilson_lower(features):
  """Compute lower Wilson bound on terminal coverage distribution.

  Returns None for small sample size, so we can ignore them.
  """
  k = terminal_mismatch(features)
  n = terminal_coverage(features)
  if n < MINIMUM_BINOMIAL_COVERAGE:
    return None
  return wilson_95_from_tally(k, n)[0]

def non_terminal_coverage(features):
  return features.non_terminal_coverage

def non_terminal_mismatch(features):
  return features.non_terminal_mismatch

def non_terminal_mismatch_wilson_lower(features):
  """Compute lower Wilson bound on non-terminal coverage distribution.

  Returns None for small sample size, so we can ignore them.
  """
  k = non_terminal_mismatch(features)
  n = non_terminal_coverage(features)
  if n < MINIMUM_BINOMIAL_COVERAGE:
    return None
  return wilson_95_from_tally(k, n)[0]


def write_track_files(cdp_dict,
                      normalization_factor,
                      shard_name,
                      output_file_base):
  """Dump the coverage data to various wiggle_0 files.

  Coverage data is scaled by the normalization factor.
  """
  direction_ext = dict()
  direction_ext['+'] = '.plus'
  direction_ext['-'] = '.minus'
  direction_files = dict()
  # NOTE(jsh): IF YOU EXPAND mapper_to_norm_factor YOU MUST...
  mapper_to_norm_factor= dict()
  mapper_to_norm_factor[terminal_coverage] = normalization_factor
  mapper_to_norm_factor[non_terminal_mismatch_wilson_lower] = 100
  mapper_to_norm_factor[terminal_mismatch_wilson_lower] = 100
  # NOTE(jsh): ...ALSO EXPAND mapper_to_ext!!  (And vice-versa.)
  mapper_to_ext = dict()
  mapper_to_ext[terminal_coverage] = ''
  mapper_to_ext[non_terminal_mismatch_wilson_lower] = '.mismatch.non_terminal'
  mapper_to_ext[terminal_mismatch_wilson_lower] = '.mismatch.terminal'
  mapper_to_file_by_direction = dict()
  mapper_to_file_by_direction['+'] = dict()
  mapper_to_file_by_direction['-'] = dict()
  for direction in mapper_to_file_by_direction:
    for mapper, ext in mapper_to_ext.iteritems():
      fname = output_file_base + ext + direction_ext[direction] + '.wig'
      fname += '.shard.' + shard_name
      f = open(fname, 'w')
      mapper_to_file_by_direction[direction][mapper] = f
      logging.info(
          'Prepping {f.name}.'.format(**vars()))
      f.write('track type=wiggle_0\n')
  for chrom in sorted(cdp_dict.keys()):
    if chrom in ['chrMito', '2-micron']:
      continue
    if '_' in chrom:
      # These should have already been cleaned out by now.
      raise RuntimeError(
          'Encountered "chromosome" {chrom} in final array.'.format(**vars()))
    logging.info(
        'Writing wiggle entries for {chrom}.'.format(**vars()))
    for direction in mapper_to_file_by_direction:
      for f in mapper_to_file_by_direction[direction].values():
        f.write('variableStep chrom={chrom}\n'.format(**vars()))
      if not cdp_dict[chrom].has_key(direction):
        continue
      values_by_pos = cdp_dict[chrom][direction]
      sorted_keys = sorted(values_by_pos.keys())
      for pos in sorted_keys:
        value = values_by_pos[pos]
        # Ignore statistics for noisy trailing read locations.
        for mapper, f in mapper_to_file_by_direction[direction].iteritems():
          # We add one to start/end here because wiggle file is one-based.
          wiggle_pos = pos + 1
          v = mapper(value)
          if v:
            norm_value = v * mapper_to_norm_factor[mapper]
            f.write('{wiggle_pos} {norm_value:f}\n'.format(**vars()))


##############################################
if __name__ == "__main__":
    sys.exit(main())
