# importing sys
import sys
import os

# adding pyWGCNA to the system path
sys.path.insert(0, '/Users/nargesrezaie/Documents/MortazaviLab/PyWGCNA')
import PyWGCNA


ABI3 = PyWGCNA.WGCNA(name='ABI3', geneExpPath='5xFAD_ABI3/expressionList_sorted', sep='\t', save=True)
ABI3 = ABI3.runWGCNA()
ABI3.save_WGCNA()

ABI3.analyseWGCNA()





