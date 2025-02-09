import os
from MMinte.MMinte.site.widget4 import totalEXRxns,createEXmodel,createReverseEXmodel, addEXMets2SpeciesEX, replaceRxns,replaceMets,createCommunityModel,allPairComModels,createAllPairs
from MMinte.MMinte.site.widget5 import calculateGR
from MMinte.MMinte.site.widget6 import evaluateInteractions
import gurobipy


my_dir=os.getcwd()
single_model_folder = my_dir+'/singlemodel/'
pair_model_folder = my_dir+'/pairedmodel/'
outputGRs=my_dir+'/outputGRs.txt'
com=pair_model_folder
data_folder = my_dir+'/data/'
analysis_folder = my_dir+'/analyse/'
#(single_model_folder,outputFolder = 'pairsList.txt')
#pair_model_filenames = allPairComModels(listOfPairs='pairsList.txt',comFolder=pair_model_folder,modelFolder=single_model_folder)
#calculateGR(diet="data/test_diet.tsv", comFolder=com, outputGRs="outputGRs.txt")
evaluateInteractions(inGRs=analysis_folder+'outputGRs.txt',outInter=analysis_folder+'/random_name.tsv')
