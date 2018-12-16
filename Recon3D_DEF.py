#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:43:44 2018

@author: marie
"""

import cobra
import re
import random

'''
FICHIER CONTENANT LES DEFINITIONS POUR ANALYSE RECON3D
'''

############### Introduction ###############
'''
Calculer fba_init et récupérer le minimal medium 1 associé.
Mettre les réactions de ce medium à 0.
Recalculer la fba et récupérer le nouveau médium.
Répéter la boucle 10 fois.
Objectif, prouver que recon est très robuste est possède donc des voies
redondantes permettant sa survie.
'''
def Minimedium(model):

    fba_init = model.optimize()
    mini_medium = cobra.medium.minimal_medium(model, fba_init.f)
    mini_medium_1 = dict(mini_medium)
    return(fba_init, mini_medium_1)


############### Tester modèle ###############

#### gene deletion analysis ####
def Drd(model, mini_medium_1):
    import cobra.test
    drd = cobra.flux_analysis.double_reaction_deletion(model,\
                                    mini_medium_1.keys(), return_frame=True)
    return(drd)  

##### tester robustesse
def Robustess(model, fba_init, mini_medium_1):
    
    minimediumL = [(None,None) for i in range(9)]
    minimediumL[0] = (mini_medium_1, fba_init.f)
    fbaT=[fba_init.f]
    nbrRKO = [0, len(mini_medium_1.keys())]
    
    for i in range(len(minimediumL)-1):
        for key in minimediumL[i][0].keys():
            model.reactions.get_by_id(key).lower_bound = 0
            model.reactions.get_by_id(key).upper_bound = 0
            #print(modelepi.reactions.get_by_id(key).bounds)       
        fba = model.optimize()
        fbaT.append(fba.f)
        new_mini_medium = cobra.medium.minimal_medium(model, fba.f)
        minimediumL[i+1]= (dict(new_mini_medium), fba.f)
        nbrRKO.append(nbrRKO[i+1]+len(dict(new_mini_medium).keys()))
        
    return(minimediumL, fbaT, nbrRKO)
    
    

###### PARTIE II : FRAGILISER RECON POUR LE RENDRE PLUS REACTIF #######
'''
idée : 
récupérer toutes les réactions extérieures de recon
récupérer le mdedium minimal
mettre toutes les reactions extérieurs à 0 sauf le minimal medium
'''

def FragilisationExt(model, reactextrecon, reac_com) :
    fba_init = model.optimize()
    
    for r in reactextrecon.keys() :
        model.reactions.get_by_id(r).lower_bound = 0
        model.reactions.get_by_id(r).upper_bound = 0
        #print(r, modelepi.reactions.get_by_id(r).bounds)
        for r, value in reac_com.items() :
            model.reactions.get_by_id(r).lower_bound = -int(value)
            model.reactions.get_by_id(r).upper_bound = 1000
    fba_fragile = model.optimize()
    return(fba_init, fba_fragile)
    
    
'''
    reactions = [["",""]]
    for r in model.reactions:
        reactions[-1][0] = str(r.id)
        reactions[-1][1] = str(model.reactions.get_by_id(r.id).bounds)
        reactions += ["",""]
 '''
def FragilisationInt(model) :
    random.seed()

    reaction = []
    reaction_aleatoire = []
    for r in model.reactions:
        if re.search("^[^EX]", r.id):
            reaction.append(r.id)

    for i in range(int(len(reaction)/100)):
        reaction_aleatoire.append(random.choice(reaction))

    for r in reaction_aleatoire :
        model.reactions.get_by_id(r).lower_bound = 0
        model.reactions.get_by_id(r).upper_bound = 0
    
    fba = model.optimize()
    
    return(reaction_aleatoire, fba)


'''
permet de verifier si les reactions de mon minimal medium sont bien
presentes dans mes reactions totales.

influenceD = {}
for r in ext_1.keys() :
    influenceD[r] = 1
for r in mini_medium_1.keys():
    if not ext_1.get(r):
        influenceD[r]=2
    else:
        influenceD[r]=influenceD[r]+2
print(influenceD)
print(len(influenceD.keys()))

influenceL = [0 for i in range(3)]
for i in range(len(influenceL)):
    for r in influenceD:
        if influenceD[r] == 1+i:
            influenceL[i] += 1
print("reactions que fusion = ",influenceL[0])
print("reactions que cell epi = ",influenceL[1])
print("reactions commun = ", influenceL[2])
print("reactions tot ext fusion = ", len(mini_medium_1.keys()))
'''

  
    
    
    
############ PARTIE III : INTERACTION FUSION ET RECON ####

def Nomenclature(dicofusioninit) :
    
    dicofusion = {}
    for key, value in dicofusioninit.items():
        dicofusion[key.replace("_e","[e]")] = value
    return(dicofusion)


def Commun(dicofusion, reactionextrecon, model):
            
    influenceD = {}
    for r in dicofusion.keys() :
        influenceD[r] = 1
    for r in reactionextrecon:
        if not dicofusion.get(r):
            influenceD[r]=2
        else:
            influenceD[r]=influenceD[r]+2

    influenceL = [0 for i in range(3)]
    for i in range(len(influenceL)):
        for r in influenceD:
            if influenceD[r] == 1+i:
                influenceL[i] += 1
                
    influence_communD = {}
    for key, value in influenceD.items():
        if value == 3 :
            influence_communD[key] = value
                
    return(influenceD, influenceL, influence_communD)
    
    
####################################
def traduire_nom_reacs(listereacs,model):
    listenom=[]
    for r in listereacs:
        reac=model.reactions.get_by_id(r)
        for m in reac.metabolites.keys():
            nom=m.name[:20]
        listenom.append(nom)
    return(listenom)

def traduire_nom_reac(reac,model):
    if type(reac) == str:
        reac=model.reactions.get_by_id(reac)
    for m in reac.metabolites.keys():
        nom=m.name    
    return(nom)
    
   


