import sys
import os

#This is the name of the bed file containing all the contigs
file = sys.argv[1]
#These are the focal species we want to mask.  It should be all species that share the trait.
to_mask = sys.argv[2]

#This is the name of the reference species that PhyloP will be computed for.  It should generally be the species with the most complete genome.
phylop_species = sys.argv[3]
#This is the folder to output everything to
folder_name = sys.argv[4]


#Make folders
o = open(file)
os.mkdir(folder_name)
os.mkdir(folder_name + "/run1")
os.mkdir(folder_name + "/All")
os.mkdir(folder_name + "/GTF")
os.mkdir(folder_name + "/Errored")

#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
#Replace REFERENCE_REPLACE with phylop_species
#Replace BED_NAME with name of the bed file containing the region to extract
command_phylop1 = '/scratch/users/astarr97/Common_Software/cactus-bin-v2.6.13/bin/hal2maf --refTargets BED_NAME --refGenome REFERENCE_REPLACE --onlyOrthologs --targetGenomes "Solenodon_paradoxus,Erinaceus_europaeus,Crocidura_indochinensis,Sorex_araneus,Uropsilus_gracilis,Condylura_cristata,Scalopus_aquaticus,Megaderma_lyra,Craseonycteris_thonglongyai,Hipposideros_armiger,Hipposideros_galeritus,Rhinolophus_sinicus,Macroglossus_sobrinus,Eidolon_helvum,Pteropus_vampyrus,Pteropus_alecto,Rousettus_aegyptiacus,Noctilio_leporinus,Pteronotus_parnellii,Mormoops_blainvillei,Carollia_perspicillata,Artibeus_jamaicensis,Anoura_caudifer,Tonatia_saurophila,Micronycteris_hirsuta,Desmodus_rotundus,Pipistrellus_pipistrellus,Eptesicus_fuscus,Lasiurus_borealis,Murina_feae,Myotis_myotis,Myotis_brandtii,Myotis_lucifugus,Myotis_davidii,Miniopterus_natalensis,Miniopterus_schreibersii,Tadarida_brasiliensis,Vicugna_pacos,Camelus_dromedarius,Camelus_ferus,Camelus_bactrianus,Hippopotamus_amphibius,Eubalaena_japonica,Balaenoptera_bonaerensis,Balaenoptera_acutorostrata,Eschrichtius_robustus,Kogia_breviceps,Platanista_gangetica,Mesoplodon_bidens,Ziphius_cavirostris,Inia_geoffrensis,Lipotes_vexillifer,Neophocaena_asiaeorientalis,Phocoena_phocoena,Delphinapterus_leucas,Monodon_monoceros,Tursiops_truncatus,Orcinus_orca,Tragulus_javanicus,Moschus_moschiferus,Bubalus_bubalis,Bos_taurus,Bos_indicus,Bos_mutus,Bison_bison,Beatragus_hunteri,Ammotragus_lervia,Hemitragus_hylocrius,Capra_hircus,Capra_aegagrus,Ovis_aries,Ovis_canadensis,Pantholops_hodgsonii,Saiga_tatarica,Okapia_johnstoni,Giraffa_tippelskirchi,Antilocapra_americana,Odocoileus_virginianus,Rangifer_tarandus,Elaphurus_davidianus,Catagonus_wagneri,Sus_scrofa,Cryptoprocta_ferox,Mungos_mungo,Suricata_suricatta,Helogale_parvula,Hyaena_hyaena,Paradoxurus_hermaphroditus,Panthera_tigris,Panthera_pardus,Panthera_onca,Felis_catus,Felis_catus_fca126,Felis_nigripes,Puma_concolor,Acinonyx_jubatus,Canis_lupus_VD,CanFam4,Canis_lupus_familiaris,Canis_lupus_dingo,Canis_lupus_orion,Lycaon_pictus,Nyctereutes_procyonoides,Vulpes_lagopus,Otocyon_megalotis,Pteronura_brasiliensis,Enhydra_lutris,Mustela_putorius,Mellivora_capensis,Ailurus_fulgens,Spilogale_gracilis,Zalophus_californianus,Odobenus_rosmarus,Leptonychotes_weddellii,Mirounga_angustirostris,Neomonachus_schauinslandi,Ailuropoda_melanoleuca,Ursus_maritimus,Manis_javanica,Manis_pentadactyla,Equus_caballus,Equus_przewalskii,Equus_asinus,Tapirus_indicus,Tapirus_terrestris,Dicerorhinus_sumatrensis,Diceros_bicornis,Ceratotherium_simum,Ceratotherium_simum_cottoni,Papio_anubis,Papio_papio,Papio_hamadryas,Papio_ursinus,Papio_kindae,Papio_cynocephalus,Lophocebus_aterrimus,Theropithecus_gelada,Cercocebus_atys,Cercocebus_lunulatus,Cercocebus_torquatus,Cercocebus_chrysogaster,Mandrillus_leucophaeus,Mandrillus_sphinx,Macaca_leonina,Macaca_silenus,Macaca_siberu,Macaca_nemestrina,Macaca_nigra,Macaca_maura,Macaca_tonkeana,Macaca_fuscata,Macaca_cyclopis,Macaca_mulatta,Macaca_fascicularis,Macaca_assamensis,Macaca_thibetana,Macaca_radiata,Macaca_arctoides,Cercopithecus_albogularis,Cercopithecus_nictitans,Cercopithecus_ascanius,Cercopithecus_cephus,Cercopithecus_petaurista,Cercopithecus_lowei,Cercopithecus_mona,Cercopithecus_pogonias,Cercopithecus_neglectus,Cercopithecus_diana,Cercopithecus_roloway,Cercopithecus_hamlyni,Miopithecus_ogouensis,Allochrocebus_preussi,Allochrocebus_lhoesti,Allochrocebus_solatus,Chlorocebus_aethiops,Chlorocebus_pygerythrus,Chlorocebus_sabaeus,Erythrocebus_patas,Allenopithecus_nigroviridis,Trachypithecus_phayrei,Trachypithecus_crepusculus,Trachypithecus_obscurus,Trachypithecus_melamera,Trachypithecus_auratus,Trachypithecus_cristatus,Trachypithecus_germaini,Trachypithecus_francoisi,Trachypithecus_leucocephalus,Trachypithecus_laotum,Trachypithecus_hatinhensis,Trachypithecus_geei,Trachypithecus_pileatus,Semnopithecus_entellus,Semnopithecus_hypoleucos,Semnopithecus_priam,Semnopithecus_schistaceus,Semnopithecus_johnii,Semnopithecus_vetulus,Pygathrix_cinerea,Pygathrix_nigripes_b,Pygathrix_nigripes_a,Rhinopithecus_strykeri,Rhinopithecus_bieti,Rhinopithecus_roxellana,Nasalis_larvatus,Presbytis_mitrata,Presbytis_comata,Colobus_polykomos,Colobus_guereza,Colobus_angolensis,Piliocolobus_tephrosceles,Piliocolobus_gordonorum,Piliocolobus_kirkii,Piliocolobus_badius,Pan_paniscus,Pan_troglodytes,Homo_sapiens,Gorilla_gorilla,Gorilla_beringei,Pongo_abelii,Pongo_pygmaeus,Hylobates_pileatus_b,Hylobates_pileatus_a,Hylobates_abbotti,Hylobates_muelleri,Hylobates_agilis,Hylobates_klossii,Hoolock_leuconedys,Nomascus_siki_b,Nomascus_siki_a,Nomascus_gabriellae,Nomascus_annamensis,Nomascus_concolor,Symphalangus_syndactylus,Mico_schneideri,Mico_argentatus,Mico_humeralifer,Callibella_humilis,Cebuella_pygmaea,Cebuella_niveiventris,Callithrix_kuhlii,Callithrix_jacchus,Callithrix_geoffroyi,Callimico_goeldii,Leontopithecus_rosalia,Leontopithecus_chrysomelas,Saguinus_bicolor,Saguinus_midas,Saguinus_oedipus,Saguinus_geoffroyi,Saguinus_mystax,Saguinus_labiatus,Saguinus_inustus,Saguinus_imperator,Leontocebus_illigeri,Leontocebus_fuscicollis,Leontocebus_nigricollis,Aotus_vociferans,Aotus_griseimembra,Aotus_nancymaae,Aotus_azarae,Aotus_trivirgatus,Cebus_olivaceus,Cebus_albifrons,Cebus_unicolor,Sapajus_macrocephalus,Sapajus_apella,Saimiri_oerstedii,Saimiri_ustus,Saimiri_cassiquiarensis,Saimiri_macrodon,Saimiri_sciureus,Saimiri_boliviensis,Alouatta_seniculus,Alouatta_juara,Alouatta_puruensis,Alouatta_macconnelli,Alouatta_nigerrima,Alouatta_caraya,Alouatta_belzebul,Alouatta_discolor,Alouatta_palliata,Ateles_marginatus,Ateles_chamek,Ateles_belzebuth,Ateles_paniscus,Ateles_geoffroyi_b,Ateles_geoffroyi_a,Lagothrix_lagothricha,Cacajao_hosomi,Cacajao_ayresi,Cacajao_melanocephalus,Cacajao_calvus,Chiropotes_sagulatus,Chiropotes_israelita,Chiropotes_albinasus,Pithecia_albicans,Pithecia_vanzolinii,Pithecia_pissinattii,Pithecia_mittermeieri,Pithecia_hirsuta,Pithecia_pithecia,Pithecia_chrysocephala,Plecturocebus_grovesi,Plecturocebus_moloch,Plecturocebus_bernhardi,Plecturocebus_dubius,Plecturocebus_caligatus,Plecturocebus_cupreus,Plecturocebus_brunneus,Plecturocebus_miltoni,Plecturocebus_cinerascens,Plecturocebus_hoffmannsi,Cheracebus_torquatus,Cheracebus_regulus,Cheracebus_lucifer,Cheracebus_lugens,Cephalopachus_bancanus,Carlito_syrichta,Tarsius_lariang,Tarsius_wallacei,Eulemur_rufus,Eulemur_collaris,Eulemur_fulvus,Eulemur_albifrons,Eulemur_sanfordi,Eulemur_flavifrons,Eulemur_macaco,Eulemur_coronatus,Eulemur_mongoz,Eulemur_rubriventer,Hapalemur_griseus,Hapalemur_meridionalis,Hapalemur_gilberti,Hapalemur_alaotrensis,Hapalemur_occidentalis,Prolemur_simus,Lemur_catta,Varecia_rubra,Varecia_variegata,Propithecus_coquerelli,Propithecus_verreauxi,Propithecus_tattersalli,Propithecus_coronatus,Propithecus_diadema,Propithecus_perrieri,Propithecus_edwardsi,Avahi_laniger,Avahi_peyrierasi,Indri_indri,Microcebus_murinus,Mirza_zaza,Cheirogaleus_major,Cheirogaleus_medius,Lepilemur_septentrionalis,Lepilemur_ankaranensis,Lepilemur_dorsalis,Lepilemur_ruficaudatus,Daubentonia_madagascariensis,Nycticebus_coucang,Nycticebus_bengalensis,Nycticebus_pygmaeus,Loris_tardigradus,Loris_lydekkerianus,Perodicticus_ibeanus,Perodicticus_potto,Arctocebus_calabarensis,Otolemur_garnettii,Otolemur_crassicaudatus,Galago_senegalensis,Galago_moholi,Galagoides_demidoff,Galeopterus_variegatus,Tupaia_tana,Tupaia_chinensis,Oryctolagus_cuniculus,Lepus_americanus,Ochotona_princeps,Ctenodactylus_gundi,Petromus_typicus,Thryonomys_swinderianus,Heterocephalus_glaber,Fukomys_damarensis,Dolichotis_patagonum,Hydrochoerus_hydrochaeris,Cavia_tschudii,Cavia_porcellus,Cavia_aperea,Dasyprocta_punctata,Cuniculus_paca,Octodon_degus,Ctenomys_sociabilis,Myocastor_coypus,Capromys_pilorides,Chinchilla_lanigera,Dinomys_branickii,Hystrix_cristata,Rattus_norvegicus,Mus_musculus,Mus_spretus,Mus_caroli,Mus_pahari,Acomys_cahirinus,Meriones_unguiculatus,Psammomys_obesus,Mesocricetus_auratus,Cricetulus_griseus,Microtus_ochrogaster,Ondatra_zibethicus,Ellobius_talpinus,Ellobius_lutescens,Sigmodon_hispidus,Peromyscus_maniculatus,Onychomys_torridus,Cricetomys_gambianus,Nannospalax_galili,Jaculus_jaculus,Allactaga_bullata,Zapus_hudsonius,Perognathus_longimembris,Dipodomys_stephensi,Dipodomys_ordii,Castor_canadensis,Xerus_inauris,Spermophilus_dauricus,Ictidomys_tridecemlineatus,Marmota_marmota,Aplodontia_rufa,Muscardinus_avellanarius,Glis_glis,Graphiurus_murinus,Dasypus_novemcinctus,Tolypeutes_matacus,Chaetophractus_vellerosus,Tamandua_tetradactyla,Myrmecophaga_tridactyla,Choloepus_didactylus,Choloepus_hoffmanni,Trichechus_manatus,Procavia_capensis,Heterohyrax_brucei,Loxodonta_africana,Microgale_talazaci,Echinops_telfairi,Chrysochloris_asiatica,Elephantulus_edwardii,Orycteropus_afer" /scratch/users/astarr97/PhyloP/hg38.447way.hal MAF_NAME_REPLACE.maf\n'

#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
command_phylop2 = "/scratch/users/astarr97/Common_Software/mafTools/bin/mafDuplicateFilter -m MAF_NAME_REPLACE.maf > MAF_NAME_REPLACE.dedup.maf\n"

#Replace CONTIG_NAME with the name of the contig
#Replace MASK_SPECIES_REPLACE with the set of species to be masked
#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
command_phylop3 = "/scratch/users/astarr97/Common_Software/phast/bin/maf_parse --features ../GTF/CONTIG_NAME.gtf --mask-features MASK_SPECIES_REPLACE MAF_NAME_REPLACE.dedup.maf > MAF_NAME_REPLACE.dedup.Masked.maf\n"

#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
command_phylop4 = "/home/groups/hbfraser/Common_Software/phast/bin/maf_parse --split 1000 --out-root MAF_PARSE_MAF_NAME_REPLACE_ MAF_NAME_REPLACE.dedup.Masked.maf\n"
command_phylop5 = "for file in MAF_PARSE_MAF_NAME_REPLACE_*.maf;\ndo\n\techo $file\n\t"
command_phylop6 = "/scratch/users/astarr97/Common_Software/phast/bin/phyloP --no-prune --chrom CHROM --msa-format MAF --method LRT --mode CONACC -d 6 --wig-scores -g /scratch/users/astarr97/astarr_scripts/AccelConvDist/fullTreeAnc239.100kb.mod $file > ${file::-4}.PhyloP.wig\n"

#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
#Replace CONTIG_NAME with the name of the contig
#Unused
command_phylop7 = "\tpython /scratch/users/astarr97/astarr_scripts/AccelConvDist/fix_wig.py ${file::-4}.PhyloP.wig CONTIG_NAME\n"

#Replace MAF_NAME_REPLACE with the name of the MAF file to create without .maf in it
command_phylop8 = "\t/scratch/users/astarr97/astarr_scripts/AccelConvDist/bin/convert2bed --do-not-sort -i wig < ${file::-4}.PhyloP.wig > ${file::-4}.PhyloP.bed\ndone\n\n"
command_phylop9 = "python /scratch/users/astarr97/astarr_scripts/AccelConvDist/get_species_support_447_new.py MAF_NAME_REPLACE.dedup.Masked.maf MAF_NAME_REPLACE.dedup.Masked.SpecSup.txt PHYLOP_SPECIES\n"
command_phylop10 = "python /scratch/users/astarr97/astarr_scripts/AccelConvDist/check_error_write.py\nchmod 777 move_error.sh\n./move_error.sh\n\n"
command_phylop11 = "cat MAF_PARSE_MAF_NAME_REPLACE*.PhyloP.bed > MAF_NAME_REPLACE.dedup.Masked.PhyloP.bed\nrm MAF_PARSE*\n\n"

def write_beg(out):
    out.write("#!/bin/bash\n#SBATCH --time=168:00:00\n#SBATCH -p hns,hbfraser\n#SBATCH --mem=16GB\n\n")
    out.write("ml load python/3.9\n\n")

#Function to write out one 10 megabase block
def write_block(out, fout, contig_name, bed_name):

    out.write(command_phylop1.replace("BED_NAME", bed_name).replace("REFERENCE_REPLACE", phylop_species).replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop2.replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop3.replace("MAF_NAME_REPLACE", fout).replace("MASK_SPECIES_REPLACE", to_mask).replace("CONTIG_NAME", contig_name))
    out.write("echo MAF_NAME_REPLACE.dedup.Masked.maf\n".replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop4.replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop5.replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop6.replace("MAF_NAME_REPLACE", fout).replace("--chrom CHROM", "--chrom " + contig_name))
    #out.write(command_phylop7.replace("MAF_NAME_REPLACE", fout).replace("CONTIG_NAME", contig_name))
    out.write(command_phylop8.replace("MAF_NAME_REPLACE", fout))
    out.write(command_phylop9.replace("MAF_NAME_REPLACE", fout).replace("PHYLOP_SPECIES", phylop_species))
    out.write(command_phylop10)
    out.write(command_phylop11.replace("MAF_NAME_REPLACE", fout))
    out.write("\n")

def write_end(out, run_num):
    out.write("cat *.SpecSup.txt > " + "Run." + str(run_num) + ".dedup.Masked.SpecSup.txt\n")
    out.write("cat *.PhyloP.bed > " + "Run." + str(run_num) + ".dedup.Masked.PhyloP.bed\n")
    out.write("mv " + "Run." + str(run_num) + ".dedup.Masked.SpecSup.txt ../All\n")
    out.write("mv " + "Run." + str(run_num) + ".dedup.Masked.PhyloP.bed ../All\n")
    

#Open the first script
out = open(folder_name + "/run1/" + "run1.sh", 'w')
write_beg(out)
base_sum = 0
c = 1
    
run_contigs = 0
for line in o:
    run_contigs += 1
    l = line.replace("\n", "").split("\t")
    gtf = open(folder_name + "/GTF/" + l[0] + ".gtf", 'w')
    gtf.write("\t".join([l[0], l[0], ".", "1", str(l[2])]) + "\n")
    gtf.close()
    
    cur_len = int(l[2])
    if cur_len < 1000000:
        bed = open(folder_name + "/run" + str(c) + "/" + l[0] + ".0.bed", 'w')
        bed.write("\t".join([l[0], l[1], l[2]]) + "\n")
        bed.close()
        
        fout = l[0] + "_" + str(0) + "-" + str(int(l[2]))
        write_block(out, fout, l[0], l[0] + "." + "0.bed")
        base_sum += cur_len
    else:
        s_end = 1000000
        while s_end < cur_len:
            if base_sum > 10000000:
                write_end(out, c)
                out.close()
                c += 1
                os.mkdir(folder_name + "/" + "run" + str(c))
                out = open(folder_name + "/" + "run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
                write_beg(out)
                base_sum = 0
                run_contigs = 0
            bed = open(folder_name + "/run" + str(c) + "/" + l[0] + "." + str((s_end - 1000000)//1000000) + ".bed", 'w')
            bed.write("\t".join([l[0], str(s_end - 1000000), str(s_end)]) + "\n")
            bed.close()
            fout = l[0] + "_" + str(s_end - 1000000) + "-" + str(s_end)
            write_block(out, fout, l[0], l[0] + "." + str((s_end - 1000000)//1000000) + ".bed")
            base_sum += 1000000
            s_end += 1000000
            
        bed = open(folder_name + "/run" + str(c) + "/" + l[0] + "." + str((s_end - 1000000)//1000000) + ".bed", 'w')
        bed.write("\t".join([l[0], str(s_end - 1000000), str(l[2])]) + "\n")
        bed.close()
        fout = l[0] + "_" + str(s_end - 1000000) + "-" + str(l[2])
        write_block(out, fout, l[0], l[0] + "." + str((s_end - 1000000)//1000000) + ".bed")
        base_sum += int(l[2]) - (s_end - 1000000)
        

    if base_sum > 10000000 or run_contigs >= 300:
        write_end(out, c)
        out.close()
        c += 1
        os.mkdir(folder_name + "/" + "run" + str(c))
        out = open(folder_name + "/" + "run" + str(c) + "/" + "run" + str(c) + ".sh", 'w')
        write_beg(out)
        base_sum = 0
        run_contigs = 0
write_end(out, c)
out.close()
