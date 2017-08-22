# -*- coding: UTF-8 -*-
from unittest import TestCase
from nose import with_setup
from niptplus.map_check import PQFile
from os import path

dirlocate = path.dirname(path.dirname(path.abspath(__file__)))
dirlocate = '/'.join([dirlocate, 'data', 'input', '151015_NS500132_0225_AHGT37BGXX'])
pqfile = PQFile('/'.join([dirlocate, '15HR03765.R1.clean.fastq.gz.PQ.txt']), 11000000)
#print pqfile
"""
每个检测中不能直接用pqfile，因为那是类外pqfile的一个全复制
"""
pqfile.parse_file()
pqfile.count_content()
pqfile.count_lessdata()
class TestPQFile(TestCase):

    def test_parse_file(self):
        #pqfile.parse_file()
        #print pqfile
        self.assertEqual(pqfile.unimap_rds, 16321440)
        #print vars(pqfile)

    def test_count_content(self):
        #print pqfile
        #pqfile.parse_file()
        #pqfile.count_content()
        self.assertEqual(pqfile.content, '-')

    def test_count_lessdata(self):
        #pqfile.parse_file()
        #pqfile.count_content()
        #pqfile.count_lessdata()
        self.assertEqual(pqfile.lessdata, 'NO')


