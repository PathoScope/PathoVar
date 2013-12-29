select gi, sub_gi, taxon, gene, stbp, edbp  from (select gi, taxon from giAnnoT where taxon = ($tax_id)) 
	as myGiAnnoT join giDelimT on myGiAnnoT.gi = giDelimT.sub_gi or myGiAnnoT.gi = giDelimT.gi
	where stbp =< ($snp_loc) and edbp => ($snp_loc) ;