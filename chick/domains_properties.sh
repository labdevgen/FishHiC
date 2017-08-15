HMAP="mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.IC"
nlog="--nolog"
#nlog=""
#fulllog="--fulllog"
fulllog=""

ChEFdomains=( 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation" 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation.shuffled" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation.shuffled"
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation" 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation.shuffled" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation.shuffled" 
	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_All.bed
#	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_ConSym.bed
#	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_Divirgent.bed
	)
#ChEFFPKM="/mnt/storage/home/vsfishman/tmp/Distr/Chick_expression/Fibs/genes.fpkm_tracking"
#CHEFGENES="/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/E1_input/ChEF-all-HindIII-100k.hm.eig"
#CHEFE1="/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/E1_input/ChEF-all-HindIII-100k.hm.eig"

echo "ChEF"
for i in "${ChEFdomains[@]}"
do 
	echo ""
#	~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i ${nlog} --FPKM ${ChEFFPKM} ${fulllog}
	~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i ${nlog} ${fulllog}
#	~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i --genes=${CHEFGENES} ${nlog} ${fulllog}
#	~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i --E1=${CHEFE1} ${nlog} ${fulllog}
done;

# echo "ChEFbstrap"
# for i in "${ChEFdomains[@]}"
# do 
	# if [[ $i == *"shuffled"* ]]; then
		# for j in {1..21}
		# do
			# ~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i$j ${nlog} --FPKM ${ChEFFPKM} ${fulllog}
		# done;
	# fi;
# done;

BloodFPKM="/mnt/storage/home/vsfishman/tmp/Distr/Chick_expression/Erythro/tophat_out/genes.fpkm_tracking"
BloodDomains=( 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation" 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation.shuffled" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation.shuffled" 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation" 
#	"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation.shuffled" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation" 
#	"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation.shuffled"
	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_All.bed
#	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_ConSym.bed
#	${HOME}/HiC/tutorial_Fishman/chick/CTCF_input_Miroslav/CTCF_Divirgent.bed
	)

echo "Blood"
HMAP="mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.IC"
#BLOODGENES="/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/E1_input/Blood-all-HindIII-100k.hm.eig"
#BLOODE1="/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/E1_input/Blood-all-HindIII-100k.hm.eig"
for i in "${BloodDomains[@]}"
do 
	# echo ""
	# ~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i ${nlog} --FPKM ${BloodFPKM} ${fulllog}
	~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i ${nlog} ${fulllog}
	# ~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i --genes ${BLOODGENES} ${nlog} ${fulllog}
	# ~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i --E1 ${BLOODE1} ${nlog} ${fulllog}
done;

# echo "Bloodbstrap"
# for i in "${BloodDomains[@]}"
# do 
	# if [[ $i == *"shuffled"* ]]; then
		# for j in {1..21}
		# do
			# ~/PFILES/Python-2.7/bin/python2.7 domains_properties.py --hmap ${HMAP} --domains $i$j ${nlog} --FPKM ${BloodFPKM} ${fulllog}
		# done;
	# fi;
# done;