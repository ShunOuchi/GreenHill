greenhill \
	-cph FALCON-Unzip_result/cns_p_ctg.fa FALCON-Unzip_result/cns_h_ctg.fa \
	-p reads/longread.fq,gz \
	-HIC reads/HIC_1.fq.gz reads/HIC_2.fq.gz \
	-o FALCON-Unzip \
	>FALCON-Unzip.log 2>&1
