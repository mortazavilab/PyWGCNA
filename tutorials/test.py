
import PyWGCNA
geneExp = '5xFAD_paper/expressionList_sorted'
pyWGCNA_5xFAD = PyWGCNA.WGCNA(name='5xFAD', geneExpPath=geneExp, sep='\t', save=True)
pyWGCNA_5xFAD.preprocess()

pyWGCNA_5xFAD.findModules()


pyWGCNA_5xFAD.addSampleInfo(path='5xFAD_paper/metaData', sep='\t')
# add color for metadata
pyWGCNA_5xFAD.setMetadataColor('Sex', {'Female': 'green',
                                       'Male': 'yellow'})
pyWGCNA_5xFAD.setMetadataColor('Genotype', {'5xFADWT': 'darkviolet',
                                            '5xFADHEMI': 'deeppink'})
pyWGCNA_5xFAD.setMetadataColor('Age', {'4mon': 'thistle',
                                       '8mon': 'plum',
                                       '12mon': 'violet',
                                       '18mon': 'purple'})
pyWGCNA_5xFAD.setMetadataColor('Tissue', {'Hippocampus': 'red',
                                          'Cortex': 'blue'})
pyWGCNA_5xFAD = pyWGCNA_5xFAD.analyseWGCNA()





