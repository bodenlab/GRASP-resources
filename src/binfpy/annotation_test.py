import annotations
import phylo
tree = phylo.parseNewick("(Paenibacillus_thiaminolyticus:4.0,(((bacterium_endosymbiont_of_Mortierella_elongata_FMR23_6:4.0,(Pandoraea_faecigallinarum:4.0,Pandoraea_vervacti:4.0,Pandoraea_oxalativorans:4.0):4.0,(Burkholderia_sp_b14:4.0,Burkholderia_sp_b13:4.0,(Burkholderia_pseudomallei_406e:4.0,Burkholderia_pseudomallei_1710a:4.0):4.0):4.0):4.0,(Chromobacterium_amazonense:4.0,(Microvirgula_sp_AG722:4.0,Microvirgula_aerodenitrificans:4.0):4.0):4.0):4.0,(Candidatus_Endobugula:4.0,Moritella_sp_PE36:4.0,(Enterovibrio_nigricans:4.0,Photobacterium_iliopiscarium:4.0,Vibrio_campbellii:4.0):4.0,(((Pantoea_sp_AMG_501:4.0,Pantoea_wallisii:4.0,Pantoea_rodasii:4.0):4.0,(Erwinia_sp_ErVv1:4.0,Erwinia_toletana:4.0,Erwinia_mallotivora:4.0):4.0):4.0,(Candidatus_Fukatsuia:4.0,Rahnella_aquatilis:4.0,(Yersinia_pekkanenii:4.0,Yersinia_entomophaga:4.0,Yersinia_mollaretii:4.0,(Yersinia_wautersii:4.0,Yersinia_similis:4.0,Yersinia_pseudotuberculosis:4.0,Yersinia_pestis:4.0):4.0,Yersinia_enterocolitica:4.0):4.0):4.0,(Cosenzaea_myxofaciens:4.0,(Photorhabdus_laumondii:4.0,Photorhabdus_bodei:4.0,Photorhabdus_sp_HUG-39:4.0,Photorhabdus_sp_CRCIA-P01:4.0,Photorhabdus_namnaonensis:4.0,Photorhabdus_khanii:4.0,Photorhabdus_heterorhabditis:4.0,Photorhabdus_temperata:4.0,Photorhabdus_asymbiotica:4.0,Photorhabdus_australis:4.0,Photorhabdus_thracensis:4.0,Photorhabdus_luminescens:4.0):4.0,(Xenorhabdus_ishibashii:4.0,Xenorhabdus_khoisanae:4.0,Xenorhabdus_mauleonii:4.0,Xenorhabdus_miraniensis:4.0,Xenorhabdus_vietnamensis:4.0,Xenorhabdus_stockiae:4.0,Xenorhabdus_szentirmaii:4.0,Xenorhabdus_budapestensis:4.0,Xenorhabdus_bovienii:4.0,Xenorhabdus_nematophila:4.0):4.0,(Proteus_sp_TJ1640:4.0,Proteus_sp_TJ1636:4.0,Proteus_sp_FJ2001126-3:4.0,Proteus_columbae:4.0,Proteus_alimentorum:4.0,Proteus_genomosp_6_str._ATCC_51471:4.0,Proteus_genomosp_4_str._ATCC_51469:4.0,Proteus_cibarius:4.0,Proteus_hauseri:4.0,Proteus_penneri:4.0,Proteus_vulgaris:4.0):4.0,(Morganella_sp_HMSC11D09:4.0,Morganella_sp_EGD-HP17:4.0,Morganella_morganii:4.0):4.0):4.0,(Escherichia_sp_ESNIH1:4.0,Mangrovibacter_phragmitis:4.0,(Enterobacter_sp_DC4:4.0,Enterobacter_sp_BIDMC_26:4.0):4.0,Kosakonia_sacchari:4.0,Pseudescherichia_vulneris:4.0):4.0):4.0,(Pseudomonas_kribbensis:4.0,Pseudomonas_lactis:4.0,Pseudomonas_paralactis:4.0,Pseudomonas_helleri:4.0,Pseudomonas_weihenstephanensis:4.0,Pseudomonas_coleopterorum:4.0,Pseudomonas_endophytica:4.0,Pseudomonas_granadensis:4.0,Pseudomonas_prosekii:4.0,Pseudomonas_brassicacearum:4.0,Pseudomonas_deceptionensis:4.0,Pseudomonas_baetica:4.0,Pseudomonas_simiae:4.0,Pseudomonas_moraviensis:4.0,Pseudomonas_batumici:4.0,Pseudomonas_antarctica:4.0,Pseudomonas_rhizosphaerae:4.0,Pseudomonas_lini:4.0,Pseudomonas_kilonensis:4.0,Pseudomonas_psychrophila:4.0,Pseudomonas_abietaniphila:4.0,Pseudomonas_thivervalensis:4.0,Pseudomonas_jessenii:4.0,Pseudomonas_plecoglossicida:4.0,Pseudomonas_agarici:4.0,(Pseudomonas_cichorii:4.0,Pseudomonas_syringae:4.0):4.0,Pseudomonas_sp:4.0,(Pseudomonas_lundensis:4.0,Pseudomonas_fragi:4.0):4.0,(Pseudomonas_poae:4.0,Pseudomonas_mediterranea:4.0,Pseudomonas_extremorientalis:4.0,Pseudomonas_orientalis:4.0,Pseudomonas_libanensis:4.0,Pseudomonas_synxantha:4.0,Pseudomonas_corrugata:4.0,Pseudomonas_fluorescens:4.0):4.0):4.0):4.0):4.0);")
# tree = phylo.readNewick("/Users/gabefoley/Dropbox/PhD/Projects/Phylo Island/Species trees/species_tree.nwk")

# tree = phylo.readNewick("/Users/gabefoley/Dropbox/PhD/Smaller Projects/GRASP tree/non_unique.nwk")

# print (tree)
# unique_tree = get_unique_tree(tree)
#
# print (unique_tree)



# working_dir = "/Users/gabefoley/Dropbox/PhD/Smaller Projects/Nexus_colouring/Read_annotations/"
#
# tree = phylo.read_nexus(working_dir + "annotation_simple.nexus")
#
# print (tree)
# print (tree.nexus_annotations.annotations)
#
# tree.swap_annotations("PDB")
#
# print (tree)
# print (tree.nexus_annotations.annotations)
#
# tree.write_to_nexus(working_dir + "output.nexus")

# nexus_annotations = annotations.NexusAnnotations2()

# for node in tree.getNodes():
#     if node.isLeaf():
#         nexus_annotations.add_annotation(node, "node", "annotation")
#
# print (nexus_annotations.annotations)


# REGEX CODE
# import re
#
# pattern = r"<person>(.*?)</person>"
# re.findall(pattern, str, flags=0) #you may need to add flags= re.DOTALL if your str is multiline