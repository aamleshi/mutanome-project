import re
import os.path
import urllib
import pprint
from Bio import SwissProt

class PDB:
	def __init__(self, name, chain, start, end, length, num_variants, variants, var_pdb_pos, var_up_pos, protein, gene, freq, project):
		self.name = name
		self.chain = chain
		self.start = start
		self.end = end
		self.length = length
		self.num_variants = num_variants
		self.variants = variants
		self.var_pdb_pos = var_pdb_pos
                self.var_up_pos = var_up_pos
                self.protein = protein
                self.gene = gene
                self.freq = freq
                self.project = project
                self.distances = []
                self.wap = -1
	def __repr__(self):
		return repr((self.name, self.chain, self.start, self.end, self.length, self.num_variants, self.variants, self.var_pdb_pos, self.var_up_pos, self.protein, self.gene, self.freq, self.project))

class ProtVar:
	def __init__(self, protid, chr, pos, consequence, freq, gene, aa_pos, aa_len, aa_var, project):
                self.protid = protid
		self.chr = chr
		self.pos = pos
		self.consequence = consequence
                self.freq = freq
		self.gene = gene
		self.aa_pos = aa_pos
		self.aa_len = aa_len
		self.aa_var = aa_var
                self.project = project
	def __repr__(self):
		return repr((self.protid, self.chr, self.pos, self.consequence, self.freq, self.gene, self.aa_pos, self.aa_len, self.aa_var, self.project))

