#!python
# -*- coding: utf-8 -*-

##############################################################
#   Adapted Script to compute HSSP and filter out false-positive
#   sequences from BLAST tabular output format 6
#
##############################################################
import sys
import math
#from logger import Logger
import os
import argparse
import unittest

i=0
minHVAL = -10000000

def hssp(minAliLength, pi, length):
	returnValue = minHVAL
	if length >= minAliLength:
		if length < 11:
			returnValue = pi - 100
		if length > 450:
			returnValue = pi - 19.5
		else:
			exp = -0.32 * (1 + math.exp(- length / 1000.0 ))
			returnValue = pi - (480.0 * (length ** exp))
	return returnValue


class BlastallOutput:
	def __init__(self,Paths):
		self.Paths = Paths
		try:
			file = open(self.Paths, 'r')
			prevFirstId = ""
			prevSecondId = ""
			hVALList = []
			minAliLength=11
			hvalRow = []
			
			lines = file.readlines()
			for x, line in enumerate(lines):
				lineSplit = line.rstrip("\n").split("\t")
				firstId = int(lineSplit[0])
				secondId = int(lineSplit[1])
				seqId = float(lineSplit[2])
				aliLen = int(lineSplit[3])
				mms = int(lineSplit[4])
				numIdenticals = int(((seqId / 100) * aliLen) + 0.5 )
				numGaps = aliLen - numIdenticals - mms
				aliLenHssp = aliLen - numGaps
				percentIdenticals = (float(numIdenticals) / aliLenHssp) * 100
				hval = hssp(minAliLength, percentIdenticals, aliLenHssp)
				#only print lines with HSSP equal or above to 5
				if hval >= 5:
					print(line)
		except:
			raise
		finally:
			i=1


class TestHSSPMethods(unittest.TestCase):
	def MyOwntest1(self):
		self.assertEqual(hssp(0,25.6,125), -0.585732313395)
		self.assertEqual(hssp(0,26.190476190476192,126), 0.0945359650884)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Compute HSSP and filter out false-positive sequences from BLAST tabular output format 6')
	parser.add_argument('-i', '--id-file', action="store", dest="id_file", default="", type=str)
	
	if len(sys.argv) < 2:
		parser.print_usage()
		sys.exit(1)
	
	args = parser.parse_args()
	u = BlastallOutput(args.id_file)

