import unittest
from treetoreads import TreeToReads

#class IntegrationTestCase(unittest.TestCase):
#    """Tests that it run`."""

 #   def can_run(self):
 #       """Does it run?"""
 #       self.assertTrue(TreeToReads(configfi='tests/input/test1.cfg',run=1))

 #   def correct_snp_num(self):
 #   	"""Is the expected number of snps getting created?"""
 #   	ttr=TreeToReads(configfi='tests/input/test1.cfg',run=0)
 #   	ttr.selectMutsites() 
#        assert(len(ttr.rands))


if __name__ == '__main__':
    ttr=TreeToReads(configfi='tests/input/test1.cfg',run=0)
    ttr.runART()
    ttr=TreeToReads(configfi='seqsim.cfg',run=0)
    ttr.runART()