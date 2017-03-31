use List::Util qw/ min max sum shuffle/;
use Statistics::TTest;
#use Statistics::Test::WilcoxonRankSum;


$,="\t";
my $ng=39;
my $nrep=3;
#@ALLgenes=("Fas","Olr1","Pax6","Chek2","FN1EDB","FN1EDA","NUMB","BCL2L1","Rac1","MSTR1","MADD","CASP2","APAF1","MINK1","MAP3K7","MAP4K3","NOTCH3","MAP4K2","PKM2","VEGFA","MCL1","CFLAR","CASP9","SYK","DIABLO","BMF","CCNE1","CCND1","H2AFY","STAT3","GADD45A","BIM","BIRC5","LMNA","SMN1","SMN2","PHF19");

@ALLgenes=("Fas","Olr1","Pax6","Chek2","FN1EDB","FN1EDA","NUMB","BCL2L1","Rac1","MSTR1","MADD","CASP2","APAF1","MINK1","MAP3K7","MAP4K3","NOTCH3","MAP4K2","PKM2","VEGFA","MCL1","CFLAR","CASP9","SYK","DIABLO","BMF","CCNE1","CCND1","H2AFY","STAT3","GADD45A","BIM_L_EL","BIM_EL" ,"BIRC5_2B","BIRC5_D3", "LMNA","SMN1","SMN2","PHF19");
@MEDWEIGHTS=(10,6,12,20,75,55,22,16,38,8,19,15,10,5,13,11,15,13,13,6,22,25,7,1,24,2,20,12,37,38,2,1,1,3,3,7,2,2,290);

@genes=@ALLgenes[0..$ng-1];

#open (IN, "LABCHIPS/labchip37_wdrugs_new.txt")||die;
#open (IN, "LABCHIPS/labchip_19_01_12.txt")||die;
#open (IN, "LABCHIPS/labchip37_wdrugs_iron_corrected.txt")||die;
#open (IN, "LABCHIPS/labchip37_02_04_13.txt")||die;
#open (IN, "LABCHIPS/labchip37_09_04_13.txt")||die;	#Fall back here for Splicing Factors
#open (IN, "LABCHIPS/labchip37_25_04_13.txt")||die;
#open (IN, "LABCHIPS/labchip37_23_05_13.txt")||die;
#open (IN, "LABCHIPS/labchip37_26_06_13.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip37_26_06_13_bare.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip37_26_06_13_SplDrugs.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip37_26_06_13_Iron.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip37_26_06_13_KATERINA.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip37_22_07_13_ALL_Iron_Hypoxia.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip_37_22_07_13.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
#open (IN, "LABCHIPS/labchip_37_01_08_14.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  
open (IN, "LABCHIPS/labchip_37_29_02_16.txt")||die;	#Added Fas Screening Factors, IK/SMU1 in HEK, Stress Factors, MK6 treatments  



while (<IN>){
$line=$_;
chomp $line;
$line =~ s/\r//g;

if ($line=~/^\d+\t+[^\t]{1,100}\t[^\t]{1,100}\t[^\t]{2,100}\t[\d\.\-]+\t/i){
@mat=split /\t/,$line;
next if $mat[3] eq "No cells";
#next if ($mat[0]>=289 );		#415 for all the new conditions, 474 for Isab,
#next if ($mat[0]>=289 && $mat[0] < 474 );		#415 for all the new conditions, 474 for Isab,
next if ($mat[0]>=289 && ($mat[0] < 470 || $mat[0] > 473) );		#415 for all the new conditions, 474 for Isab,

$SF{$mat[3]}=defined;
$Nid{$mat[3]}=$mat[0];
$Nid2SF{$mat[0]}=$mat[3];
$LUHRMANN{$mat[3]}=$mat[1];
$ROB_REED{$mat[3]}=$mat[2];


$rep{$mat[3]}=$rep{$mat[3]}+$nrep;
	for $n (0..$ng-1){
	$gene=$genes[$n];
		for $r (0..$nrep-1){
		
		next if ($mat[4+2*$nrep*$n+2*$r]<$MEDWEIGHTS[$n]*0.1 && $mat[4+2*$nrep*$n]+$mat[4+2*$nrep*$n+2]+$mat[4+2*$nrep*$n+4] > $MEDWEIGHTS[$n]*0.4);	#Before only first condition. Changed 28-06-13
		next if $mat[4+2*$nrep*$n+2*$r]<0.01; 
		$mat[4+2*$nrep*$n+2*$r]=0.01 if $mat[4+2*$nrep*$n+2*$r]<0.01;		
		
		push @{$ALL_TOTALS{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r];
		push @{$ALL_INDICES{$genes[$n]}},$mat[5+2*$nrep*$n+2*$r];
		
		push @{$TOTALS{$mat[3]}{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r];
		push @{$SINDEX{$mat[3]}{$genes[$n]}},$mat[5+2*$nrep*$n+2*$r];
				
		if ($n>=0) {push @{$WEIGHTS{$mat[3]}{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r]; push @{$GWEIGHTS{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r]}
		elsif ($mat[4+2*$nrep*$n+2*$r]<=10 )	{push @{$WEIGHTS{$mat[3]}{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r]; push @{$GWEIGHTS{$genes[$n]}},$mat[4+2*$nrep*$n+2*$r]}
		elsif ($mat[4+2*$nrep*$n+2*$r]>10 && $mat[4+2*$nrep*$n+2*$r]<=20)	{push @{$WEIGHTS{$mat[3]}{$genes[$n]}},10; push @{$GWEIGHTS{$genes[$n]}},10}
		elsif ($mat[4+2*$nrep*$n+2*$r]>20 && $mat[4+2*$nrep*$n+2*$r]<=30)	{push @{$WEIGHTS{$mat[3]}{$genes[$n]}},30-$mat[4+2*$nrep*$n+2*$r]; push @{$GWEIGHTS{$genes[$n]}},30-$mat[4+2*$nrep*$n+2*$r]}
		else {push @{$WEIGHTS{$mat[3]}{$genes[$n]}},0.01; push @{$GWEIGHTS{$genes[$n]}},0.01}				
		}
	}
}

}
close IN;

  
foreach $sf (keys %SF){
	foreach $n (0..$ng-1){
		$gene=$genes[$n];
		next if (sum( @{$WEIGHTS{$sf}{$gene}} )<0.5 );	#was 0.01. Changed 02-08-13
		@CONTROLS=(@{$SINDEX{"Mock SiRNA"}{$gene}}, @{$SINDEX{"Untransfected"}{$gene}});
		$MD{$sf}{$gene}=median(\@{$SINDEX{$sf}{$gene}})+0.1;
		$CMED{$gene}=median(\@CONTROLS)+0.1;
		$CMAD{$gene}=MADN(\@CONTROLS)+0.1;
		$MN{$sf}{$gene}=mean(\@{$SINDEX{$sf}{$gene}})+0.1;
		$WM{$sf}{$gene}=wmean(\@{$SINDEX{$sf}{$gene}},\@{$WEIGHTS{$sf}{$gene}})+0.1;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if scalar(@{$WEIGHTS{$sf}{$gene}})<3;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if sum( @{$WEIGHTS{$sf}{$gene}} )< $MEDWEIGHTS[$n]*0.125;;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if min( @{$WEIGHTS{$sf}{$gene}} )< $MEDWEIGHTS[$n]*0.5;;			
		$STDV{$sf}{$gene}=stddev(@{$SINDEX{$sf}{$gene}})+0.1;
		$WSDV{$sf}{$gene}=wstd_dev(\@{$SINDEX{$sf}{$gene}},\@{$WEIGHTS{$sf}{$gene}})+0.1;
		$MADN{$sf}{$gene}=MADN(\@{$SINDEX{$sf}{$gene}})+0.1;

		#The following sufficient statistics are only for informational purposed and used nowhere downstream:
		if ($sf eq "Untransfected"){
		%{$UNT_suff{$gene}}=(
    	'count' =>scalar(@{$WEIGHTS{$sf}{$gene}}),
    	'mean' =>$MD{$sf}{$gene},
    	'variance' =>$MADN{$sf}{$gene}**2);
		}
		elsif ($sf=~/Sudemycin/i){
		%{$SUD_suff{$gene}}=(
    	'count' =>scalar(@{$WEIGHTS{$sf}{$gene}}),
    	'mean' =>$MD{$sf}{$gene},
    	'variance' =>$MADN{$sf}{$gene}**2);
		}
		elsif ($sf=~/Splic/i){
		%{$SPL_suff{$gene}}=(
    	'count' =>scalar(@{$WEIGHTS{$sf}{$gene}}),
    	'mean' =>$MD{$sf}{$gene},
    	'variance' =>$MADN{$sf}{$gene}**2);
		}
		elsif ($sf=~/iron/i){
		%{$IRON_suff{$gene}}=(
    	'count' =>scalar(@{$WEIGHTS{$sf}{$gene}}),
    	'mean' =>$MD{$sf}{$gene},
    	'variance' =>$MADN{$sf}{$gene}**2);
		}
	}
}


for $g (0..$ng-1){
$gene=$genes[$g];

$GMD[$g]=median(\@{$ALL_INDICES{$genes[$g]}})+0.1;	
$GMN[$g]=mean(\@{$ALL_INDICES{$genes[$g]}})+0.1;
$GWM[$g]=wmean(\@{$ALL_INDICES{$genes[$g]}},\@{$GWEIGHTS{$genes[$g]}})+0.1;

$GMADN[$g]=MADN(\@{$ALL_INDICES{$genes[$g]}})+0.1;
$GSTDV[$g]=stddev(@{$ALL_INDICES{$genes[$g]}})+0.1;
$GWSDV[$g]=wstd_dev(\@{$ALL_INDICES{$genes[$g]}},\@{$GWEIGHTS{$genes[$g]}})+0.1;

$GMADN[$g]=0.5 * $GSTDV[$g]  if $GMADN[$g]<0.5 && $GSTDV[$g]>=1;

%{$Msuff{$genes[$g]}}=(
        'count' =>scalar(@{$GWEIGHTS{$genes[$g]}}),
        'mean' =>$GMD[$g],
        'variance' =>$GMADN[$g]**2         );
}



#Start printing of Z_P table:
$,="\t";
print "SF";
for (0..$#genes){
print "\t$genes[$_]\t$genes[$_]";	
}
print "\n";

		#0		1		2		3		4			5		6		7	8		9		10			11		12		13			14		15		16		17		18	19
@COND=("GEN","diron1","siron1","diron2","siron2","stopo","dtopo","Sud","Spl","SF3B1_LV","ATT_P27","HYPOX","Actin","LV_Drugs","OsmStr","DYRK1A","HEK","FasKD","MK6","Stress",
		#20			#21				#22			#23			#24		#25			#26			#27			#28		#29			#30		#31		#32	
		"SplF_PCNT","IK_SMU1_T24","IK_SMU1_T48","Mito","LV_Drugs2","LV_Drugs3","CYRILLE","LDrugs","Luisa_16","CycTunBort","RECQ5","Konig","H2O2_4","H2O2_24","HS_43","AS","AS_43");

@{$CN{$COND[0]}}=("Mock SiRNA","Untransfected");	#GENERAL
@{$CN{$COND[1]}}=("CN_iron_DMSO_rep1","DFOA_rep1","Hemin_rep1","CPX_rep1");	#IRON drugs
@{$CN{$COND[2]}}=("CN_iron_Mock_rep1","CN_iron_Mock_rep1","ACO1_rep1","FTL_rep1","ZNF326");	#IRON siRNAs												
@{$CN{$COND[3]}}=("CN_iron_DMSO_rep2","DFOA_rep2","Hemin_rep2","CPX_rep2");	#IRON drugs 2nd batch
@{$CN{$COND[4]}}=("CN_iron_Mock_rep2","CN_iron_Mock_rep2","ACO1_rep2","FTL_rep2");	#IRON siRNAs 2nd batch
@{$CN{$COND[5]}}=("Mock_siRNA_TOP","Mock_siRNA_TOP","TOP1","TOPBP1","TOP2A");	#siRNA Topoisomerases
@{$CN{$COND[6]}}=("DNA_damage_drugs_0h","DNA_damage_drugs_0h","CPT_3h","Etopo_3h","Doxor_3h","TOP2ci_3h");	#Drugs for Topoisomerases
@{$CN{$COND[7]}}=("CN Sudemycin","CN spliceostatin","Scrambled");	#Sudemycin
@{$CN{$COND[8]}}=("CN spliceostatin");	#Spliceostatin
@{$CN{$COND[9]}}=("Scrambled");	#SF3B1 siRNA 2nd batch (Luisa)
@{$CN{$COND[10]}}=("Mock_siRNA_P27","Mock_siRNA_P27","P27","roscovitine","Hypoxia");	#Attila P27
@{$CN{$COND[11]}}=("Mock_siRNA_P27","Mock_siRNA_P27","P27","roscovitine","Hypoxia","Hypoxia");	#Hypoxia
@{$CN{$COND[12]}}=("ActD_0h","Sta_0h");	#Actinomycin
@{$CN{$COND[13]}}=("Control_DMSO2","Control_DMSO1","Cyclosporin","Meayamycin_L","TG003");	#Luisa's 2nd batch of Drugs
@{$CN{$COND[14]}}=("OS_CN");	#Caterina's Osmotic Stress
@{$CN{$COND[15]}}=("CN_DYRK1A");	#Caterina's DYRK1A
@{$CN{$COND[16]}}=("CL_HEK");	#HEK Context for IK, SMU1
@{$CN{$COND[17]}}=("Mock_siRNA_Fas");	#Fas Regulators from Screening KDs 
@{$CN{$COND[18]}}=("GFP_OE");	#Caterina's MKK6
@{$CN{$COND[19]}}=("CN_DMSO");	#Other Stress Factors
@{$CN{$COND[20]}}=("Control_siRNA","Control_72h","Control_siRNA","Control_72h","Control_48h");	#New Splicing Factors and Paper Controls / #IK_SMU_T72h / #Keiko's ZRSR2
@{$CN{$COND[21]}}=("Control_24h","Control_24h","Control_24h","Control_24h","Control_48h","Control_48h","IK_24h","SMU1_24h");	#IK_SMU_T24h
@{$CN{$COND[22]}}=("Control_48h","Control_48h","Control_24h","Control_72h");	#IK_SMU_T48h
@{$CN{$COND[23]}}=("Control_Mito","Control_Mito","Drug_Mito","Control_48h","Control_24h");	#Eli's Mitochondrial Treatments
@{$CN{$COND[24]}}=("DMSO_8h","DMSO_8h","CyclosporinA_8h","DMSO_3h","Control_Mito","Control_48h");	#Luisa's Km136 and new CyclosporinA
@{$CN{$COND[25]}}=("DMSO_3h","DMSO_3h","Meayamycin_3h","DMSO_8h","Control_Mito","Control_48h");	#Luisa's Sudem, Meayam new timepoints
@{$CN{$COND[26]}}=("CN_cyrille","untransfected_cyrille","CN_juan");	#Cyrille

@{$CN{$COND[27]}}=("LD_028","LD_296","LD_297","LD_344353","LD_343878","LD_DMSO","LD_DMSO");	#Luhrman drugs
@{$CN{$COND[28]}}=("PPIH_mix","DMSO_8h_2","DMSO_8h_2","D_TG003","D_TG009","DMSO_8h_3","DMSO_8h_3","Control_si_2a","Control_si_2b","Control_si_3a","Control_si_3b","GFP_1","MMTP1");	#Luisa_16

for $i (434..457){
push @{$CN{$COND[29]}}, $Nid2SF{$i};	#CycTunBort
push @{$CN{$COND[29]}}, $Nid2SF{$i} if $Nid2SF{$i}=~/UNT/;	#CycTunBort
}

for $i (458..469){
push @{$CN{$COND[30]}}, $Nid2SF{$i};	#RECQ5,Kinetin,Km249
}
push @{$CN{$COND[30]}},("siCTRL_1","untransfected_16","siCTRL_2","CTRL_3","DMSO_16");


@{$CN{$COND[31]}}=("Contr_16","Contr_16","Contr_16","PTB_KD","nPTB_KD");	#JulianKonig.

for $i (474..481){
push @{$CN{$COND[32]}}, $Nid2SF{$i};	#H2O2_4h
}
push  @{$CN{$COND[32]}},("ctr_4h","ctr_4h","ctr_4h","ctr_4h","ctr_24h","ctr_37");

for $i (474..481){
push @{$CN{$COND[33]}}, $Nid2SF{$i};	#H2O2_24h
}
push @{$CN{$COND[33]}},("ctr_4h","ctr_24h","ctr_24h","ctr_24h","ctr_24h","ctr_37");

for $i (474..481){
push @{$CN{$COND[34]}}, $Nid2SF{$i};	#HShock 43
}
push @{$CN{$COND[34]}},("ctr_37","ctr_37","ctr_37","21AS_37","21AS_37","ctr_4h","ctr_24h");

for $i (474..481){
push @{$CN{$COND[35]}}, $Nid2SF{$i};	#21AS
}
push @{$CN{$COND[35]}},("ctr_37","ctr_37","ctr_37","ctr_37","ctr_4h","ctr_24h");

for $i (474..481){
push @{$CN{$COND[36]}}, $Nid2SF{$i};	#21AS_43
}
push @{$CN{$COND[36]}},("ctr_37","ctr_37","21AS_37","21AS_37","21AS_37","ctr_4h","ctr_24h");

#$,="\t#";
#print "\n\n";
#print @{$CN{$COND[36]}};
#print "\n\n";
#die;

foreach $sf (keys %SF){
	foreach $n (0..$ng-1){
	$gene=$genes[$n];
		for $cnd (0..$#COND){
			@MED_CNT=();
			@CONTROLS=@{$CN{$COND[$cnd]}};
			for (0..$#CONTROLS) {push @MED_CNT,$MD{$CONTROLS[$_]}{$gene}};
			$med_CNT=median(\@MED_CNT);
			$dist=( ($med_CNT-$MD{$sf}{$gene})/$GMADN[$n] );
			push @{$CN_DIST{$cnd}{$sf}},$dist;
			push @{$CN_DIST_SQ{$cnd}{$sf}},$dist**2+0.005;
			$DIST{$cnd}{$gene}+=$dist**2;
		}
	@CNSdiron=($MD{"CN_iron_DMSO_rep1"}{$gene},$MD{"DFOA_rep1"}{$gene},$MD{"Hemin_rep1"}{$gene},$MD{"CPX_rep1"}{$gene});	
	$CN_dIRON_d=(median(\@CNSdiron)-$MD{$sf}{$gene})/$GMADN[$n];	#IRON drugs
	push @{$CN_dIRON_DIST{$sf}},($CN_dIRON_d)**2+0.005;#  unless ($GMD[$n]<2 || $GMD[$n]>98 );
	push @{$CN_dIRON_DISTr{$sf}},($CN_dIRON_d);#  unless ($GMD[$n]<2 || $GMD[$n]>98 );
	$IRONd{$gene}+=$CN_dIRON_d**2;
	}
}



for $cnd (0..$#COND){
	foreach $sf (keys %SF){
	$AbsMedDist=abs(median(\@{$CN_DIST{$cnd}{$sf}}));
	$MedDist=median(\@{$CN_DIST{$cnd}{$sf}});
	$MedDistSQ=median(\@{$CN_DIST_SQ{$cnd}{$sf}});
	$weight=(1/sqrt($MedDistSQ));  
		if ($MedDistSQ<0.25  && $AbsMedDist<0.125){	#0.3 0.15
		#print "$COND[$cnd]\t$sf\t$AbsMedDist\t$MedDistSQ\t$weight\n";
		push @{$WEIs_CNs{$COND[$cnd]}},$weight;
			foreach $n (0..$ng-1){
			$gene=$genes[$n];
			#push @{$VALs_CNs{$COND[$cnd]}{$gene}},$MD{$sf}{$gene};
			push @{$VALs_CNs{$COND[$cnd]}{$gene}},($MD{$sf}{$gene}-0.099) x int($weight) if $MD{$sf}{$gene}>0;			
			}

		}
		if ($MedDistSQ<0.1  && $AbsMedDist<0.05){
			foreach $n (0..$ng-1){
			$gene=$genes[$n];
			push @{$VALs_CNs{$COND[$cnd]}{$gene}}, @{$SINDEX{$sf}{$gene}} x int($weight);			
			}

		}
	}
}


#die;



foreach $sf (keys %SF){
next if ($Nid{$sf}>=289 && $Nid{$sf} < 415 );
#next if ($Nid{$sf}>=289);

#print "$Nid{$sf}\t$sf";
print "$sf";

	foreach $n (0..$ng-1){
		$gene=$genes[$n];
		if (sum( @{$WEIGHTS{$sf}{$gene}} )>0.5 ){	#was 0.01. Changed 02-08-13
		$MD{$sf}{$gene}=median(\@{$SINDEX{$sf}{$gene}})+0.1;
		$MN{$sf}{$gene}=mean(\@{$SINDEX{$sf}{$gene}})+0.1;
		$WM{$sf}{$gene}=wmean(\@{$SINDEX{$sf}{$gene}},\@{$WEIGHTS{$sf}{$gene}})+0.1;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if scalar(@{$WEIGHTS{$sf}{$gene}})<3;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if sum( @{$WEIGHTS{$sf}{$gene}} )< $MEDWEIGHTS[$n]*0.5;
		$MD{$sf}{$gene}=$WM{$sf}{$gene} if min( @{$WEIGHTS{$sf}{$gene}} ) < $MEDWEIGHTS[$n]*0.125 ;
		
		$STDV{$sf}{$gene}=stddev(@{$SINDEX{$sf}{$gene}})+0.1;
		$WSDV{$sf}{$gene}=wstd_dev(\@{$SINDEX{$sf}{$gene}},\@{$WEIGHTS{$sf}{$gene}})+0.1;
		$MADN{$sf}{$gene}=MADN(\@{$SINDEX{$sf}{$gene}})+0.1;

		$base=$GMD[$n];
		$base=wmean(\@{$VALs_CNs{"GEN"}{$gene}},\@{$WEIs_CNs{"GEN"}});
		$base=median(\@{$VALs_CNs{"GEN"}{$gene}});

		$var=$GMADN[$n];
		#$var=16*$var if $gene eq "SMN1";	#Rescaling for SMN1 to match previous (paper) data
		#print "\nCHECKIT: $gene\t$base\t$var\t$GMD[$n]\t$GMADN[$n]\tM:$MD{'Mock SiRNA'}{$gene}\tU:$MD{'Untransfected'}{$gene}\n";



		$,="\t";
		if ($Nid{$sf}>=375){
		$cnd=26 if ($Nid{$sf}>=375 && $Nid{$sf}<=382);
		$cnd=20 if (($Nid{$sf}>=383 && $Nid{$sf}<=386) || ($Nid{$sf}>=395 && $Nid{$sf}<=406));
		$cnd=21 if ($Nid{$sf}>=387 && $Nid{$sf}<=390);
		$cnd=22 if ($Nid{$sf}>=391 && $Nid{$sf}<=394);
		$cnd=23 if ($Nid{$sf}>=407 && $Nid{$sf}<=408);
		$cnd=24 if ($Nid{$sf}>=409 && $Nid{$sf}<=411);
		$cnd=25 if ($Nid{$sf}>=412 && $Nid{$sf}<=414);
		#$cnd=26 if ($Nid{$sf}>=412 && $Nid{$sf}<=414);
		$cnd=27 if ($Nid{$sf}>=415 && $Nid{$sf}<=420);
		$cnd=28 if ($Nid{$sf}>=421 && $Nid{$sf}<=433);
		$cnd=29 if ($Nid{$sf}>=434 && $Nid{$sf}<=457);
		$cnd=30 if ($Nid{$sf}>=458 && $Nid{$sf}<=469);
		$cnd=31 if ($Nid{$sf}>=470 && $Nid{$sf}<=473);
		$cnd=32 if ($Nid{$sf}>=474 && $Nid{$sf}<=475);
		$cnd=33 if ($Nid{$sf}>=476 && $Nid{$sf}<=477);
		$cnd=34 if ($Nid{$sf}>=478 && $Nid{$sf}<=479);
		$cnd=35 if ($Nid{$sf}==480);
		$cnd=36 if ($Nid{$sf}==481);
		$cond=$COND[$cnd];

		@vals=@{$VALs_CNs{$cond}{$gene}};
		@weis=@{$WEIs_CNs{$cond}};
		$est=median(\@vals)+0.1;
		#$est=wmean(\@vals,\@weis);
		#print "\n>>$Nid{$sf}\t$cond\t$sf\t$gene\t$MD{$sf}{$gene}\t$#vals\t$est\t$base\tVVVVV";
		#print @vals;
		#print "VVVVVV<<<\n";
		$corr=$MD{$sf}{$gene} * $base/ $est;
		$corr=0.1 if $corr<0.1;
		$corr=100.1 if $corr>100.1;
		$RZval{$sf}{$gene}=($MD{$sf}{$gene} - $est)  / ($var*($est/$base));  #
		$RZval{$sf}{$gene}="NAn" if abs($est-$base)/$var > 5;
		%{$Csuff{$gene}}=%{$Msuff{$gene}};	
		}




		else {
		$Zval{$sf}{$gene}=($MN{$sf}{$gene}-$GMN[$n])/$GSTDV[$n];
		#$RZval{$sf}{$gene}=($MD{$sf}{$gene}-$GMD[$n])/$var;		#Robust Zvalue
		$RZval{$sf}{$gene}=($MD{$sf}{$gene}-$base)/$var;		#Robust Zvalue		
		#print "\nAAAAA $sf\t$gene\t$MD{$sf}{$gene}\t$base\t$var\n";
		$WZval{$sf}{$gene}=($WM{$sf}{$gene}-$GWM[$n])/$GWSDV[$n];		#Weighted Zvalue
		%{$Csuff{$gene}}=%{$Msuff{$gene}};
		}









#		print "\n\n\n$AAA $sf\t$gene\t$RZval{$sf}{$gene}\t$MD{$sf}{$gene}\t$GMD[$n]\n";
			if ( scalar(@{$WEIGHTS{$sf}{$gene}})>1 ){
					%{$suff{$sf}{$gene}}=(
		        	'count' =>scalar(@{$WEIGHTS{$sf}{$gene}}),
		        	'mean' =>$MD{$sf}{$gene},
		        	'variance' =>$MADN{$sf}{$gene}**2
		         	);
						$ttest = new Statistics::TTest::Sufficient; 
						$ttest -> load_data(\%{$Csuff{$gene}},\%{$suff{$sf}{$gene}});
						$ttest->set_significance(99.99);	#Assume Equal Variances Unless Overwhelming Evidence Suggest Otherwise		
						$probT{$sf}{$gene} = $ttest->{t_prob};
				}
			
			else {
				%{$suff{$sf}{$gene}}=(
	        	'count' =>2,
	        	'mean' =>$WM{$sf}{$gene},
	        	'variance' =>100
	         	);
					$ttest = new Statistics::TTest::Sufficient; 
					$ttest -> load_data(\%{$Csuff{$gene}},\%{$suff{$sf}{$gene}});
					$ttest->set_significance(99.99);	#Assume Equal Variances Unless Overwhelming Evidence Suggest Otherwise		
					$probT{$sf}{$gene} = $ttest->{t_prob};								
				}
		#$wilcox_test = Statistics::Test::WilcoxonRankSum->new();
		#$wilcox_test->load_data(\@{$SINDEX{$sf}{$gene}},\@{$SINDEX{"Untransfected"}{$gene}});
		#$probW{$sf}{$gene} = $wilcox_test->probability();
		}
		
		else{
			($MD{$sf}{$gene},$MN{$sf}{$gene},$WM{$sf}{$gene},
			$STDV{$sf}{$gene},$WSDV{$sf}{$gene},$MADN{$sf}{$gene},
			$Zval{$sf}{$gene},$RZval{$sf}{$gene},$WZval{$sf}{$gene},
			$probT{$sf}{$gene})=("NAn","NAn","NAn","NAn","NAn","NAn","NAn","NAn","NAn","NAn");  }
	
		
	
	#print "\t$RZval{$sf}{$gene}";	

	print "\t$RZval{$sf}{$gene}\t$probT{$sf}{$gene}";	
	}
print "\n";
}

















sub mean {
  my $ar=shift;
  my $result;
  foreach (@$ar) {$result+= $_ }
  return $result / scalar(@$ar);
}

sub median {
my ($ar) = shift;
my $count = scalar @$ar;
# Sort a COPY of the array, leaving the original untouched
my @array = sort { $a <=> $b } @$ar;
if ($count % 2) {
return $array[int($count/2)];
} else {
return ($array[$count/2] + $array[$count/2 - 1]) / 2;
}
}



sub stddev {
	return unless @_;
	return 0 unless @_ > 1;
	return sqrt variance (@_);
}



sub variance {
	return unless @_;
	return 0 unless @_ > 1;
	my $mean= mean (\@_);
	return (sum map { ($_ - $mean)**2 } @_) / $#_;
}




sub wmean {
	my $ar=shift;
	my $wr=shift;
	if ($#ar!=$#wr) {print "\n!ERROR!: a,w vectors are not of the same length in Weighted Mean Calculation\n"; die;}
	if (min(@$wr)<0) {print "\n!ERROR!: You cannot have negative entries in w vector in Weighted Mean Calculation\n"; die;}
	if (sum(@$wr)==0) {print "\n!ERROR!: No non-zero entries in w vector in Weighted Mean Calculation\n"; die;}
	my $result;
	for (0..$#$ar) {$result+=$$ar[$_]*$$wr[$_]}
	return $result/sum(@$wr);
}



sub wstd_dev {
	my $ar=shift;
	my $wr=shift;
	my $nz=0;
	return 0 unless @$wr > 1;
	if ($#ar!=$#wr) {print "\n!ERROR!: a,w vectors are not of the same length in Weighted Mean Calculation\n"; die;}
	if (min(@$wr)<0) {print "\n!ERROR!: You cannot have negative entries in w vector in Weighted Mean Calculation\n"; die;}
	if (sum(@$wr)==0) {print "\n!ERROR!: No non-zero entries in w vector in Weighted Mean Calculation\n"; die;}
	my $wmean = wmean($ar,$wr);
	for (0..$#$ar) {$result+=$$wr[$_]*(($$ar[$_]-$wmean)**2);$nz++ if $$wr[$_]>0}
    $result=$result / ( (($nz-1)/$nz)*sum(@$wr) );
return sqrt($result);
}




sub MAD{
 my($x) = @_;
 my @ad = ();
 my($median) = median($x);
 foreach my $xi (@$x)
 {
  my($adi) = abs($xi - $median);
  push(@ad,$adi);
 }
 my($mad) = median(\@ad) + 0.0;
 return $mad;
}



# The rescaled Median Absolute Deviation
sub MADN{
 my($x) = @_; 
 my $mad = MAD($x)*1.4826; #$mad /= qnorm(0.75);
 return $mad;
}







