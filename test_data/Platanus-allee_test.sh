greenhill \
	-c Platanus-allee_result/out_nonBubbleOther.fa \
	-b Platanus-allee_result/out_primaryBubble.fa Platanus-allee_result/out_secondaryBubble.fa \
	-IP1 reads/PE_*.fq.gz \
	-OP2 reads/MP5k_*.fq.gz \
	-OP3 reads/MP9k_*.fq.gz \
	-p reads/longread.fq.gz \
	-HIC reads/HIC_1.fq.gz reads/HIC_2.fq.gz \
	-o Platanus-allee \
	>Platanus-allee.log 2>&1
