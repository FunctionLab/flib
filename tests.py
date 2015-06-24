'''
Created on Feb 4, 2013

@author: atadych
'''
import unittest
from xdslparser import CptNodesHolder
from tempfile import TemporaryFile


class TestCptNodes(unittest.TestCase):

    
    def setUp(self):
        """
        Creates 2 example files with cpt nodes entry for parsing 
        """
        self.tmpFile1=TemporaryFile()

        self.tmpFile1.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        self.tmpFile1.write('<smile version="1.0" id="Unnamed" numsamples="1000">\n')
        self.tmpFile1.write('<nodes>\n')
        self.tmpFile1.write('<cpt id="FR">\n')
        self.tmpFile1.write('<state id="FR00" />\n')
        self.tmpFile1.write('<state id="FR01" />\n')
        self.tmpFile1.write('<probabilities>0.9869151711463928 0.01308482233434916</probabilities>\n')
        self.tmpFile1.write('</cpt>\n')
        self.tmpFile1.write('</nodes>\n')
        self.tmpFile1.write('</smile>\n')
        self.tmpFile1.seek(0)
                
        self.states1 = ['FR00','FR01']   
        self.pos1 = [0.9869151711463928]
        self.neg1 = [0.01308482233434916]
        

        xmlFile2 = '<smile version="1.0" id="Unnamed" numsamples="1000">\n\
                    <nodes>\n\
                    <cpt id="biogrid">\n\
                        <state id="biogrid00" />\n\
                        <state id="biogrid01" />\n\
                        <state id="biogrid02" />\n\
                        <state id="biogrid03" />\n\
                        <parents>FR</parents>\n\
                        <probabilities>0.9999541044235229 4.262849324732088e-05 2.680929355847184e-06 5.241560643298726e-07 0.9927160143852234 0.00523342052474618 0.001486745080910623 0.0005638484144583344</probabilities>\n\
                    </cpt>\n\
                </nodes>\
                </smile>'
        self.tmpFile2=TemporaryFile()        
        self.tmpFile2.write(xmlFile2)
        self.tmpFile2.seek(0)
        self.states2 = ['biogrid00','biogrid01','biogrid02','biogrid03']   
        self.pos2 = [0.9999541044235229,4.262849324732088e-05,2.680929355847184e-06,5.241560643298726e-07]
        self.neg2 = [0.9927160143852234,0.00523342052474618,0.001486745080910623,0.0005638484144583344]
        self.parents =['FR']        

    def tearDown(self):
        self.tmpFile1.close()
        self.tmpFile2.close()

    def testParseFile1(self):
        nodeHolder = CptNodesHolder(filename=self.tmpFile1)
        self.assertIsNotNone(nodeHolder)
        node = nodeHolder.get_node('FR')
        self.assertIsNotNone(node)
        self.assertEqual(node.states, self.states1,'State list is not correct')
        self.assertEqual(node.parents,[], 'Parents list is incorrect')
        self.assertEqual(nodeHolder.get_probabilities('FR'),[self.pos1,self.neg1])

    def testParseFile2(self):
        nodeHolder = CptNodesHolder(filename=self.tmpFile2)
        self.assertIsNotNone(nodeHolder)
        node = nodeHolder.get_node('biogrid')
        self.assertIsNotNone(node)
        self.assertEqual(node.states, self.states2,'State list is not correct')
        self.assertEqual(node.parents,['FR'], 'Parents list is incorrect')
        self.assertEqual(nodeHolder.get_probabilities('biogrid'),[self.pos2,self.neg2])

    
    def testParseFile1CheckNonExisting(self):
        nodeHolder = CptNodesHolder(filename=self.tmpFile2)
        self.assertIsNotNone(nodeHolder)
        node = nodeHolder.get_node('INCORRECT')
        self.assertIsNone(node)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
