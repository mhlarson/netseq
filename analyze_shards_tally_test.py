#!/usr/bin/env python

# Author: John Hawkins (jsh)

import analyze_shards as lib

import logging
import unittest
import os
import sys
import tempfile

import HTSeq

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class TestShardProcessing(unittest.TestCase):
  def setUp(self):
    pass

  def testFindIfNot(self):
    tagger = lambda x: x % 2 == 1
    things = [1, 3, 2, 5]
    self.assertEqual(lib.find_if_not(tagger, things), 2)
    things = [1, 3, 5, 9, 7]
    self.assertEqual(lib.find_if_not(tagger, things), 5)

  def testRevComplementMdString(self):
    md_string = '2C4^GT0T0'
    self.assertEqual('0A0^AC4G2', lib.revcomp_md(md_string))

  def testAlignedFragmentsTrivial(self):
    md_string = '2C4^GG0A0'
    sequence = 'TTTTTTTT'
    generator = lib.aligned_fragments(md_string, sequence)
    self.assertEqual(generator.next(), ('TT', 'TT'))
    self.assertEqual(generator.next(), ('T', 'C'))
    self.assertEqual(generator.next(), ('TTTT', 'TTTT'))
    self.assertEqual(generator.next(), ('', 'GG'))
    self.assertEqual(generator.next(), ('', ''))
    self.assertEqual(generator.next(), ('T', 'A'))
    self.assertEqual(generator.next(), ('', ''))
    self.assertRaises(StopIteration, generator.next)

  def testComputeWilsonInterval95(self):
    (low, high) = lib.compute_wilson_interval_95(0.2, 100)
    self.assertAlmostEqual(low, 0.133367, 4)
    self.assertAlmostEqual(high, 0.288829, 4)
    (low, high) = lib.compute_wilson_interval_95(0.25, 4)
    self.assertAlmostEqual(low, 0.045587, 4)
    self.assertAlmostEqual(high, 0.699358, 4)

  def testScoredMismatches(self):
    md_string = '2C4^GG0A0'
    sequence = 'TTTTTTTT'
    scores = [80, 70, 60, 50, 40, 30, 20, 10]
    generator = lib.scored_mismatches(md_string, sequence, scores)
    self.assertEqual(generator.next(), (2, 'T', 'C', 0.999999))
    self.assertEqual(generator.next(), (7, 'T', 'A', 0.9))
    self.assertRaises(StopIteration, generator.next)

  def testCorrectlyComplementAlignmentPlus(self):
    line_chunks = ('foo', '0', 'chrI', '10101', '255', '10M', '*', '0', '0',
                   'GTTTTTTTTT', '0123456789',
                   'XA:i:2', 'MD:Z:2C4^GT0A2', 'NM:i:2')
    line = '\t'.join(line_chunks)
    alignment = HTSeq.SAM_Alignment.from_SAM_line(line)
    (iv, read, md_string) = lib.correctly_complement_alignment(alignment)
    self.assertEqual(iv.strand, '-')
    self.assertEqual(read.seq, 'AAAAAAAAAC')
    self.assertEqual(read.qual[0], 24)
    self.assertEqual(read.qual[-1], 15)
    self.assertEqual(md_string, '2T0^AC4G2')

  def testCorrectlyComplementAlignmentMinus(self):
    line_chunks = ('foo', '16', 'chrI', '10101', '255', '10M', '*', '0', '0',
                   'GTTTTTTTTT', '0123456789',
                   'XA:i:2', 'MD:Z:2C4^GT0A2', 'NM:i:2')
    line = '\t'.join(line_chunks)
    alignment = HTSeq.SAM_Alignment.from_SAM_line(line)
    (iv, read, md_string) = lib.correctly_complement_alignment(alignment)
    self.assertEqual(iv.strand, '+')
    self.assertEqual(read.seq, 'GTTTTTTTTT')
    self.assertEqual(read.qual[0], 15)
    self.assertEqual(read.qual[-1], 24)
    self.assertEqual(md_string, '2C4^GT0A2')

  def testFeaturesForAlignmentPlus(self):
    line_chunks = ('foo', '0', 'chrI', '10101', '255', '10M', '*', '0', '0',
                   'GTTTTTTTTT', '0123456789',
                   'XA:i:2', 'MD:Z:2C4^GT0A2', 'NM:i:2')
    line = '\t'.join(line_chunks)
    alignment = HTSeq.SAM_Alignment.from_SAM_line(line)
    features = lib.features_for_alignment(alignment)
    # for pos in sorted(features.keys()):
    #   stats = features[pos]
    #   logging.info('{pos}: {stats}'.format(**vars()))
    pos = ('chrI', '-', 10100)
    self.assertEqual(0, features[pos].non_terminal_coverage)
    self.assertEqual(1, features[pos].terminal_coverage)
    self.assertEqual(0, features[pos].terminal_mismatch)
    pos = ('chrI', '-', 10101)
    self.assertEqual(1, features[pos].non_terminal_coverage)
    self.assertEqual(0, features[pos].non_terminal_mismatch)
    self.assertEqual(0, features[pos].terminal_coverage)
    pos = ('chrI', '-', 10102)
    self.assertEqual(1, features[pos].non_terminal_coverage)
    self.assertEqual(1, features[pos].non_terminal_mismatch)
    self.assertEqual(0, features[pos].terminal_coverage)

  def testFeaturesForAlignmentMinus(self):
    line_chunks = ('foo', '16', 'chrI', '10101', '255', '10M', '*', '0', '0',
                   'GTTTTTTTTT', '0123456789',
                   'XA:i:2', 'MD:Z:2C4^GT0A2', 'NM:i:2')
    line = '\t'.join(line_chunks)
    alignment = HTSeq.SAM_Alignment.from_SAM_line(line)
    features = lib.features_for_alignment(alignment)
    # for pos in sorted(features.keys()):
    #   stats = features[pos]
    #   logging.info('{pos}: {stats}'.format(**vars()))
    pos = ('chrI', '+', 10109)
    self.assertEqual(0, features[pos].non_terminal_coverage)
    self.assertEqual(1, features[pos].terminal_coverage)
    self.assertEqual(0, features[pos].terminal_mismatch)
    pos = ('chrI', '+', 10108)
    self.assertEqual(1, features[pos].non_terminal_coverage)
    self.assertEqual(0, features[pos].non_terminal_mismatch)
    self.assertEqual(0, features[pos].terminal_coverage)
    pos = ('chrI', '+', 10107)
    self.assertEqual(1, features[pos].non_terminal_coverage)
    self.assertEqual(1, features[pos].non_terminal_mismatch)
    self.assertEqual(0, features[pos].terminal_coverage)

  def testTallyWriter(self):
    test_dir = tempfile.mkdtemp(dir='/tmp')
    logging.info('Working in {0}.'.format(test_dir))
    file_base = os.path.join(test_dir, 'test')
    loc = ('chrXV', '-', 149453)
    ps = lib.PositionalStats()
    ps.non_terminal_coverage = 9
    ps.non_terminal_mismatch = 6
    ps.non_terminal_mismatch_by_base['A'] = 6
    ps.non_terminal_mismatch_by_base['C'] = 0
    ps.non_terminal_mismatch_by_base['G'] = 0
    ps.non_terminal_mismatch_by_base['T'] = 0
    ps.terminal_coverage = 0
    ps.terminal_mismatch = 0
    ps.terminal_mismatch_by_base['A'] = 0
    ps.terminal_mismatch_by_base['C'] = 0
    ps.terminal_mismatch_by_base['G'] = 0
    ps.terminal_mismatch_by_base['T'] = 0
    writer = lib.TallyFileWriter('chrXV', file_base)
    writer.start_chrom('chrXV')
    writer.write(loc, ps)
    writer.close()
    with open(file_base + '.minus.tallies.shard.chrXV', 'r') as tally_file:
      header = tally_file.readline()
      tally = tally_file.read().strip().split()
      self.assertEqual(15, len(tally))
      self.assertEqual('chrXV', tally[0])
      self.assertEqual('-', tally[1])
      self.assertEqual(149454, int(tally[2]))
      self.assertEqual(9, int(tally[3]))
      self.assertEqual(6, int(tally[4]))
      self.assertEqual(6, int(tally[5]))
      self.assertEqual(0, int(tally[6]))
      self.assertEqual(0, int(tally[7]))
      self.assertEqual(0, int(tally[8]))
      self.assertEqual(0, int(tally[9]))
      self.assertEqual(0, int(tally[10]))
      self.assertEqual(0, int(tally[11]))
      self.assertEqual(0, int(tally[12]))
      self.assertEqual(0, int(tally[13]))
      self.assertEqual(0, int(tally[14]))


  def testTrackWriter(self):
    test_dir = tempfile.mkdtemp(dir='/tmp')
    logging.info('Working in {0}.'.format(test_dir))
    file_base = os.path.join(test_dir, 'test')
    loc = ('chrXV', '-', 149453)
    ps = lib.PositionalStats()
    ps.non_terminal_coverage = 10
    ps.non_terminal_mismatch = 6
    ps.non_terminal_mismatch_by_base['A'] = 6
    ps.non_terminal_mismatch_by_base['C'] = 0
    ps.non_terminal_mismatch_by_base['G'] = 0
    ps.non_terminal_mismatch_by_base['T'] = 0
    ps.terminal_coverage = 0
    ps.terminal_mismatch = 0
    ps.terminal_mismatch_by_base['A'] = 0
    ps.terminal_mismatch_by_base['C'] = 0
    ps.terminal_mismatch_by_base['G'] = 0
    ps.terminal_mismatch_by_base['T'] = 0
    writer = lib.TrackFileWriter(190, 'chrXV', file_base)
    writer.start_chrom('chrXV')
    writer.write(loc, ps)
    writer.close()
    file_name = file_base + '.mismatch.non_terminal.minus.wig.shard.chrXV'
    with open(file_name, 'r') as track_file:
      header = track_file.readline() + track_file.readline()
      logging.info('Header: ' + header)
      line = track_file.read().strip().split()
      self.assertEqual(2, len(line))
      self.assertEqual(149454, int(line[0]))


if __name__ == '__main__':
  unittest.main()