#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 14:38:38 2018

@author: marie
"""

import Recon3D_DEF as rd
import cobra
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import re

model = cobra.io.load_matlab_model('Recon3D_301.mat')

############### Introduction ###############


fba_init, minimedium1 = rd.Minimedium(model)
print("La fba initiale de Recon3D est : ","\n", fba_init.f, "\n \n",\
      "Le minimal medium associé est : ","\n", minimedium1)



############### Tester modèle ###############


##### double reaction deletion sur le minimal medium montre la robustesse
'''
Effectue une double reaction deletion sur le minimal medium
'''
drd = rd.Drd(model, minimedium1)



##### tester robustesse
'''
Faire KO sur le minimal medium. Recupérer fba. Récupérer le nouveau minimal
medium. Recommencer.
'''
minimediumL, fbaT, nbrRKO = rd.Robustess(model, fba_init, minimedium1)
'''
minimediumL = [[dico_minimedium1, fba1], [dico_minimedium2, fba2], ...]
nbr_reactions = [len(dico_minimedium1), len(dico_minimedium2), ...]
'''
for i in range(len(minimediumL)):
    print("Réactions du minimal medium :")
    for r in minimediumL[i][0].keys():
        print(model.reactions.get_by_id(r).name)
    print("\n","Nous allons les bloquer.")
    print("La fba avec ce minimal medium est :")
    print(minimediumL[i][1],"\n")

# nombre total de réactions KO en fonction des étapes de la loop 
nbrRKO.pop(-1)
plt.bar([i+1 for i in range(9)], height = nbrRKO)
plt.ylabel('nombre de réactions d\'échange KO') 
plt.xlabel('step')
plt.title("Nombre total de réactions d'échange KO \n" \
          "en fonction des étapes de la loop")
plt.savefig("recon3D_KO.png")
plt.show()

# fba de recon3D en fonction du nombres de réactions mises KO
plt.plot([i+1 for i in range(9)], fbaT)
plt.ylabel('fba de recon3D') 
plt.xlabel('step : KO progressif du minimal medium')
plt.title('fba de recon3D en fonction du nombre de réactions  KO')
plt.savefig('recon3D_fba.png')
plt.show()

for r in model.reactions:
    model.reactions.get_by_id(r.id).lower_bound = -1000
    model.reactions.get_by_id(r.id).upper_bound = 1000

######### Réactions en commun entre communauté bactérienne et recon3D #########


#Medium après simulation fusion (R ext de fusion)
R_B = {'EX_12ppd__S_e': 1000.0, 'EX_glc__D_e': 1000.0, 'EX_h2o_e': 1000.0, 'EX_23cump_e': 1000.0, 'EX_26dap__M_e': 919.2, 'EX_h_e': 6000.0, 'EX_leu__L_e': 1000.0, 'EX_3amp_e': 1000.0, 'EX_3cmp_e': 1000.0, 'EX_3gmp_e': 1000.0, 'EX_3mb_e': 1000.0, 'EX_3ump_e': 1000.0, 'EX_4hbz_e': 1000.0, 'EX_ala__L_e': 1000.0, 'EX_cl_e': 995.768, 'EX_arab__L_e': 1000.0, 'EX_4abut_e': 1000.0, 'EX_acald_e': 1000.0, 'EX_ac_e': 5990.249, 'EX_chol_e': 1000.0, 'EX_pi_e': 2900.389, 'EX_adn_e': 1000.0, 'EX_ins_e': 1000.0, 'EX_nh4_e': 1000.0, 'EX_mal__L_e': 1000.0, 'EX_alaala_e': 1000.0, 'EX_asn__L_e': 1000.0, 'EX_cys__L_e': 1000.0, 'EX_ala__D_e': 1000.0, 'EX_gln__L_e': 929.639, 'EX_gly_e': 5978.0, 'EX_arg__L_e': 1, 'EX_alltn_e': 1, 'EX_amp_e': 971.908, 'EX_cgly_e': 1000.0, 'EX_malt_e': 1000.0, 'EX_malttr_e': 1000.0, 'EX_starch1200_e': 1000.0, 'EX_amylose300_e': 1000.0, 'EX_arab__D_e': 1000.0, 'EX_asp__L_e': 1000.0, 'EX_k_e': 944.289, 'EX_cellb_e': 1, 'EX_meoh_e': 1000.0, 'EX_but_e': 1000.0, 'EX_bz_e': 1000.0, 'EX_ca2_e': 995.768, 'EX_cell4_e': 1000.0, 'EX_chols_e': 1000.0, 'EX_cit_e': 846.381, 'EX_mg2_e': 995.613, 'EX_mn2_e': 995.969, 'EX_cobalt2_e': 995.996, 'EX_zn2_e': 995.985, 'EX_succ_e': 1000.0, 'EX_for_e': 1000.0, 'EX_cmp_e': 1000.0, 'EX_co2_e': 6000.0, 'EX_cu2_e': 995.968, 'EX_o2_e': 928.731, 'EX_glu__L_e': 1000.0, 'EX_cytd_e': 1000.0, 'EX_dad_2_e': 1000.0, 'EX_din_e': 1000.0, 'EX_dca_e': 1000.0, 'EX_dcyt_e': 1000.0, 'EX_dgsn_e': 1000.0, 'EX_dha_e': 1000.0, 'EX_drib_e': 1000.0, 'EX_ser__D_e': 1000.0, 'EX_duri_e': 1000.0, 'EX_cell500_e': 1000.0, 'EX_enlipa_e': 1000.0, 'EX_galur_e': 1000.0, 'EX_etha_e': 1000.0, 'EX_eths_e': 1000.0, 'EX_etoh_e': 1, 'EX_h2_e': 1000.0, 'EX_fe2_e': 995.511, 'EX_enter_e': 1000.0, 'EX_feenter_e': 1000.0, 'EX_fru_e': 1000.0, 'EX_fuc__L_e': 1000.0, 'EX_fum_e': 1000.0, 'EX_g3pe_e': 1000.0, 'EX_gal_e': 1000.0, 'EX_galman4_e': 1000.0, 'EX_galman600_e': 1000.0, 'EX_galman6_e': 1000.0, 'EX_galt_e': 1000.0, 'EX_glcman4_e': 1000.0, 'EX_glcman600_e': 1000.0, 'EX_glcman6_e': 1000.0, 'EX_glcur_e': 1000.0, 'EX_glucan1500_e': 1000.0, 'EX_glucan6_e': 1000.0, 'EX_glucan4_e': 1000.0, 'EX_glyclt_e': 1000.035, 'EX_glycogen1500_e': 1000.0, 'EX_glyc_e': 1829.132, 'EX_gm2lipa_e': 1000.0, 'EX_gsn_e': 1000.0, 'EX_gua_e': 1, 'EX_h2s_e': 1000.0, 'EX_4hphac_e': 1000.0, 'EX_hxa_e': 1000.0, 'EX_ile__L_e': 1000.0, 'EX_inost_e': 1000.0, 'EX_kdo2lipid4_e': 1000.0, 'EX_Larab_e': 1000.0, 'EX_lcts_e': 1000.0, 'EX_lac__L_e': 1000.0, 'EX_lmn2_e': 1000.0, 'EX_lmn30_e': 1000.0, 'EX_lyx__L_e': 1000.0, 'EX_malthp_e': 650.477, 'EX_manttr_e': 1000.0, 'EX_udcpdp_e': 1000.0, 'EX_no2_e': 1, 'EX_no3_e': 1, 'EX_uri_e': 1000.0, 'EX_octa_e': 1000.0, 'EX_pac_e': 1000.0, 'EX_ppa_e': 1000.0, 'EX_progly_e': 1000.0, 'EX_pullulan1200_e': 1000.0, 'EX_pyr_e': 1000.0, 'EX_so4_e': 1000.0, 'EX_val__L_e': 1000.0, 'EX_raffin_e': 1000.0, 'EX_spmd_e': 1646.2939999999999, 'EX_rmn_e': 1000.0, 'EX_salchs4fe_e': 1000.0, 'EX_salchs4_e': 1000.0, 'EX_sbt__D_e': 1000.0, 'EX_s_e': 1133.404, 'EX_tartr__D_e': 1000.0, 'EX_thymd_e': 1000.0, 'EX_tol_e': 1000.0, 'EX_udcpp_e': 1000.0, 'EX_xmp_e': 1000.0, 'EX_xtsn_e': 1000.0, 'EX_xylan4_e': 1000.0, 'EX_xyl3_e': 1000.0, 'EX_xyl__D_e': 1000.0, 'EX_xylb_e': 1000.0, 'EX_2pglyc_e': 1000.0, 'EX_LalaDgluMdap_e': 1000.0, 'EX_LalaDgluMdapDala_e': 1000.0, 'EX_dextrin_e': 1000.0, 'EX_araban__L_e': 1000.0, 'EX_abt__L_e': 1000.0, 'EX_acnam_e': 1000.0, 'EX_ribflv_e': 995.991, 'EX_ade_e': 1000.0, 'EX_ser__L_e': 1000.0, 'EX_thr__L_e': 981.841, 'EX_anhgm_e': 1000.0, 'EX_phe__L_e': 987.232, 'EX_pro__L_e': 1000.0, 'EX_ala_B_e': 1000.0, 'EX_cell6_e': 1000.0, 'EX_coa_e': 1000.0, 'EX_fol_e': 1000.0, 'EX_glcn_e': 1000.0, 'EX_glyc3p_e': 1, 'EX_glyb_e': 1000.0, 'EX_his__L_e': 986.073, 'EX_hxan_e': 1000.0, 'EX_melib_e': 1000.0, 'EX_met__L_e': 978.87, 'EX_mmet_e': 1000.0, 'EX_uaagmda_e': 998.0, 'EX_ump_e': 1000.0, 'EX_udcpo5_e': 1000.0, 'EX_ptrc_e': 1000.0, 'EX_urate_e': 1819.0, 'EX_trp__L_e': 1000.0, 'EX_tyr__L_e': 1000.0, 'EX_rib__D_e': 1000.0, 'EX_skm_e': 1000.0, 'EX_sucr_e': 1000.0, 'EX_ura_e': 1000.0, 'EX_xan_e': 1000.0, 'EX_3oxoadp_e': 1000.0, 'EX_4abz_e': 1000.0, 'EX_acser_e': 1000.0, 'EX_airs_e': 1000.0, 'EX_akg_e': 1000.0, 'EX_fe3_e': 995.841, 'EX_lys__L_e': 960.021, 'EX_orn_e': 1000.0, 'EX_chor_e': 1000.0, 'EX_csn_e': 1000.0, 'EX_lac__D_e': 1000.0, 'EX_g1p_e': 1, 'EX_gthrd_e': 852.721, 'EX_glu__D_e': 1000.0, 'EX_glyald_e': 1, 'EX_icit_e': 1000.0, 'EX_id3acald_e': 1000.0, 'EX_man1p_e': 1, 'EX_murein5p5p_e': 1000.0, 'EX_murein5px4p_e': 1000.0, 'EX_murein5px3p_e': 1000.0, 'EX_murein5p5p5p_e': 1000.0, 'EX_murein5px4px4p_e': 1000.0, 'EX_murein4px4p_e': 1000.0, 'EX_murein4px4px4p_e': 1000.0, 'EX_murein5p4p_e': 1000.0, 'EX_murein4p4p_e': 1000.0, 'EX_murein4p3p_e': 1000.0, 'EX_murein5p3p_e': 1000.0, 'EX_rnam_e': 995.954, 'EX_pdima_e': 1000.0, 'EX_pnto__R_e': 995.988, 'EX_rbt_e': 1000.0, 'EX_urea_e': 1000.0, 'EX_srb__L_e': 1000.0, 'EX_thm_e': 995.995, 'EX_thym_e': 995.998, 'EX_tre_e': 1000.0}

#Dico des réactions echanges fusion
'''
Recupérer les reactions extérieurs de la fusion : influence le milieu.
La nomenclature des reactions changent pour s'adapter à recon3D.
"_e" -> " [e]
'''
dicofusion = rd.Nomenclature(R_B)

#Récupérer les réactions extérieurs de Recon3D
ext_r_recon = []
for r in model.reactions:
    if re.search("^EX", r.id) :
        ext_r_recon.append(r.id)
print(len(ext_r_recon))

#Quelles sont les réactions en communs?
'''
influenceD = dico des réactions extérieurs de communauté bactérienne et recon3D.
Si value = 1 alors réactions uniques à la communauté
Si value = 2 alors réactions uniques à recon3D
Si value = 3 alors réactions partagées

influenceL :
influenceL[0] = nombre de réactions uniques à la communauté
influenceL[1] = nombre de réactions uniques à recon3D
influenceL[2] = nombre de réactions partagées

influence_communD = dico contenant que les réactions partagées
'''
influenceD, influenceL, influence_communD = rd.Commun(dicofusion, ext_r_recon, model)

# Diagramme de Venn permettant la visualisation des réactions
v=venn2(subsets = influenceL, set_labels=("Communauté bactérienne","Recon3D"))

#Quel est le nom des metabolites en commun?
nom_metabolites_communs = rd.traduire_nom_reacs(influence_communD.keys(), model)



############### Epurer le modèle ###############
'''
Bloquer toutes les réactions extérieurs de Recon3D à 0.
Sauf les réactions extérieurs en commun entre recon3D et communauté bactérienne
'''

for r in ext_r_recon :
    model.reactions.get_by_id(r).lower_bound = 0
    model.reactions.get_by_id(r).upper_bound = 0
for r in influence_communD.keys() :
    model.reactions.get_by_id(r).lower_bound = -500
    model.reactions.get_by_id(r).upper_bound = 1000
fba_fragile = model.optimize()
print(fba_fragile.f)


############### Plus de détail sur le medium ###############
'''
Connaitre les métabolites qui augmentent, diminuent ou restent constant en
présence de la communauté bactérienne
'''

#Récupérer les métabolites partagés qui augmentent dans le medium en présence
#dela communauté bactériennes, ceux qui diminuent ou ceux qui restent inchangés
R_augm_com = {}
R_dim_com = {}
R_same_com = {}
for r in influence_communD.keys():
    for f, medium in dicofusion.items():
        if f == r:
            if medium > 1100:
                R_augm_com[r] = medium
            if medium < 900:
                R_dim_com[r] = medium
            if medium > 900 and medium < 1100:
                R_same_com[r] = medium
nom_metabolites_augm = rd.traduire_nom_reacs(R_augm_com.keys(), model)
nom_metabolites_dim = rd.traduire_nom_reacs(R_dim_com.keys(), model)
nom_metabolites_same = rd.traduire_nom_reacs(R_same_com.keys(), model)
print(nom_metabolites_augm)
print(nom_metabolites_dim)
print(nom_metabolites_same)



############### Influence du medium sur Recon3D ###############

#Metabolites qui augmentent après fusion
bounds_augm = [0 for i in range(1000)]
T_augm = [0 for i in range(1000)]

for i in range(1000):
    for r in R_augm_com.keys():
        model.reactions.get_by_id(r).lower_bound = -(i+1)
    bounds_augm[i] = (i+1)
    fba_augm = model.optimize()
    T_augm[i] = fba_augm.f
    
plt.plot(bounds_augm,T_augm)
plt.ylabel('Biomass production of Recon3D') 
plt.xlabel('metabolites uptake')
plt.title('Biomasse de recon3D en fonction des métabolites augmentés \n'\
          'dans le milieu par la communauté bactérienne')
plt.savefig("augm.png")
plt.show()
 
#Metabolites qui diminuent après fusion
bounds_dim = [0 for i in range(1000)]
T_dim = [0 for i in range(1000)]

for i in range(1000):
    for r in R_dim_com.keys():
        model.reactions.get_by_id(r).lower_bound = -(i+1)
    bounds_dim[i] = (-i+1)
    fba_dim = model.optimize()
    T_dim[i] = fba_dim.f
    
plt.plot(bounds_dim,T_dim)
plt.ylabel('Biomass production of Recon3D') 
plt.xlabel('metabolites uptake')
plt.title('Biomasse de recon3D en fonction des métabolites diminués \n'\
          'dans le milieu par la communauté bactérienne')
plt.savefig("dim.png")
plt.show()