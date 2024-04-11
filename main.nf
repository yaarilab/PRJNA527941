$HOSTNAME = ""
params.outdir = 'results'  

evaluate(new File("${params.projectDir}/nextflow_header.config"))
params.metadata.metadata = "${params.projectDir}/tools.json"

if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).into{g_6_mate_g_8;g_6_mate_g0_0;g_6_mate_g0_5;g_6_mate_g0_7;g_6_mate_g1_11;g_6_mate_g1_9;g_6_mate_g1_12;g_6_mate_g3_12;g_6_mate_g3_15;g_6_mate_g3_19;g_6_mate_g4_10;g_6_mate_g4_12;g_6_mate_g4_14}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_7_reads_g_8}
 } else {  
	g_7_reads_g_8 = Channel.empty()
 }



process unizp {

input:
 set val(name),file(reads) from g_7_reads_g_8
 val mate from g_6_mate_g_8

output:
 set val(name),file("*.fastq")  into g_8_reads0_g0_0

script:

readArray = reads.toString().split(' ')	
R1 = readArray[0]
R2 = readArray[1]

"""
case "$R1" in
*.gz | *.tgz ) 
        gunzip -c $R1 > R1.fastq
        ;;
*)
        cp $R1 ./R1.fastq
        echo "$R1 not gzipped"
        ;;
esac

case "$R2" in
*.gz | *.tgz ) 
        gunzip -c $R2 > R2.fastq
        ;;
*)
        cp $R2 ./R2.fastq
        echo "$R2 not gzipped"
        ;;
esac
"""
}


process Filter_Sequence_Quality_filter_seq_quality {

input:
 set val(name),file(reads) from g_8_reads0_g0_0
 val mate from g_6_mate_g0_0

output:
 set val(name), file("*_${method}-pass.fast*")  into g0_0_reads0_g3_12
 set val(name), file("FS_*")  into g0_0_logFile1_g0_5
 set val(name), file("*_${method}-fail.fast*") optional true  into g0_0_reads22
 set val(name),file("out*") optional true  into g0_0_logFile33

script:
method = params.Filter_Sequence_Quality_filter_seq_quality.method
nproc = params.Filter_Sequence_Quality_filter_seq_quality.nproc
q = params.Filter_Sequence_Quality_filter_seq_quality.q
n_length = params.Filter_Sequence_Quality_filter_seq_quality.n_length
n_missing = params.Filter_Sequence_Quality_filter_seq_quality.n_missing
fasta = params.Filter_Sequence_Quality_filter_seq_quality.fasta
//* @style @condition:{method="quality",q}, {method="length",n_length}, {method="missing",n_missing} @multicolumn:{method,nproc}

if(method=="missing"){
	q = ""
	n_length = ""
	n_missing = "-n ${n_missing}"
}else{
	if(method=="length"){
		q = ""
		n_length = "-n ${n_length}"
		n_missing = ""
	}else{
		q = "-q ${q}"
		n_length = ""
		n_missing = ""
	}
}

readArray = reads.toString().split(' ')	

fasta = (fasta=="true") ? "--fasta" : ""

if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R1_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	FilterSeq.py ${method} -s $R2 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_R2_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}else{
	R1 = readArray[0]
	"""
	FilterSeq.py ${method} -s $R1 ${q} ${n_length} ${n_missing} --nproc ${nproc} --log FS_${name}.log --failed ${fasta} 2>&1 | tee -a out_${R1}_FS.log
	"""
}


}


process Filter_Sequence_Quality_parse_log_FS {

input:
 set val(name), file(log_file) from g0_0_logFile1_g0_5
 val mate from g_6_mate_g0_5

output:
 file "*table.tab"  into g0_5_logFile0_g0_7, g0_5_logFile0_g0_16

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray}  -f ID QUALITY
"""

}


process Filter_Sequence_Quality_report_filter_Seq_Quality {

input:
 val mate from g_6_mate_g0_7
 file log_files from g0_5_logFile0_g0_7

output:
 file "*.rmd"  into g0_7_rMarkdown0_g0_16


shell:

if(mate=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]

	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	quality_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, quality_log_2, titles=plot_titles, sizing="figure")
	```
	
	`r figures("quality")`
		
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	name = R1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read")#params$quality_titles
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("quality", 
	        paste("Mean Phred quality scores for",  plot_titles[1],
	              "The dotted line indicates the average quality score under which reads were removed."))
	```
	
	```{r, echo=FALSE}
	quality_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Quality Scores
	
	Quality filtering is an essential step in most sequencing workflows. pRESTO’s
	FilterSeq tool remove reads with low mean Phred quality scores. 
	Phred quality scores are assigned to each nucleotide base call in automated 
	sequencer traces. The quality score (`Q`) of a base call is logarithmically 
	related to the probability that a base call is incorrect (`P`): 
	$Q = -10 log_{10} P$. For example, a base call with `Q=30` is incorrectly 
	assigned 1 in 1000 times. The most commonly used approach is to remove read 
	with average `Q` below 20.
	
	```{r, echo=FALSE}
	plotFilterSeq(quality_log_1, titles=plot_titles[1], sizing="figure")
	```
	
	`r figures("quality")`
	
	EOF
	
	open OUT, ">FSQ_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Filter_Sequence_Quality_presto_render_rmarkdown {

input:
 file rmk from g0_7_rMarkdown0_g0_16
 file log_file from g0_5_logFile0_g0_16

output:
 file "*.html"  into g0_16_outputFileHTML00
 file "*csv" optional true  into g0_16_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_assemble_pairs {

input:
 set val(name),file(reads) from g0_0_reads0_g3_12
 val mate from g_6_mate_g3_12

output:
 set val(name),file("*_assemble-pass.f*")  into g3_12_reads0_g1_11
 set val(name),file("AP_*")  into g3_12_logFile1_g3_15
 set val(name),file("*_assemble-fail.f*") optional true  into g3_12_reads_failed22
 set val(name),file("out*")  into g3_12_logFile33

script:
method = params.Assemble_pairs_assemble_pairs.method
coord = params.Assemble_pairs_assemble_pairs.coord
rc = params.Assemble_pairs_assemble_pairs.rc
head_fields_R1 = params.Assemble_pairs_assemble_pairs.head_fields_R1
head_fields_R2 = params.Assemble_pairs_assemble_pairs.head_fields_R2
failed = params.Assemble_pairs_assemble_pairs.failed
fasta = params.Assemble_pairs_assemble_pairs.fasta
nproc = params.Assemble_pairs_assemble_pairs.nproc
alpha = params.Assemble_pairs_assemble_pairs.alpha
maxerror = params.Assemble_pairs_assemble_pairs.maxerror
minlen = params.Assemble_pairs_assemble_pairs.minlen
maxlen = params.Assemble_pairs_assemble_pairs.maxlen
scanrev = params.Assemble_pairs_assemble_pairs.scanrev
minident = params.Assemble_pairs_assemble_pairs.minident
evalue = params.Assemble_pairs_assemble_pairs.evalue
maxhits = params.Assemble_pairs_assemble_pairs.maxhits
fill = params.Assemble_pairs_assemble_pairs.fill
aligner = params.Assemble_pairs_assemble_pairs.aligner
// align_exec = params.Assemble_pairs_assemble_pairs.// align_exec
// dbexec = params.Assemble_pairs_assemble_pairs.// dbexec
gap = params.Assemble_pairs_assemble_pairs.gap
usearch_version = params.Assemble_pairs_assemble_pairs.usearch_version
assemble_reference = params.Assemble_pairs_assemble_pairs.assemble_reference
head_seqeunce_file = params.Assemble_pairs_assemble_pairs.head_seqeunce_file
//* @style @condition:{method="align",alpha,maxerror,minlen,maxlen,scanrev}, {method="sequential",alpha,maxerror,minlen,maxlen,scanrev,ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="reference",ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec} {method="join",gap} @multicolumn:{method,coord,rc,head_fields_R1,head_fields_R2,failed,nrpoc,usearch_version},{alpha,maxerror,minlen,maxlen,scanrev}, {ref_file,minident,evalue,maxhits,fill,aligner,align_exec,dbexec}, {gap} 

// args
coord = "--coord ${coord}"
rc = "--rc ${rc}"
head_fields_R1 = (head_fields_R1!="") ? "--1f ${head_fields_R1}" : ""
head_fields_R2 = (head_fields_R2!="") ? "--2f ${head_fields_R2}" : ""
failed = (failed=="false") ? "" : "--failed"
fasta = (fasta=="false") ? "" : "--fasta"
nproc = "--nproc ${nproc}"

scanrev = (scanrev=="false") ? "" : "--scanrev"
fill = (fill=="false") ? "" : "--fill"

// align_exec = (align_exec!="") ? "--exec ${align_exec}" : ""
// dbexec = (dbexec!="") ? "--dbexec ${dbexec}" : ""


ref_file = (assemble_reference!='') ? "-r ${assemble_reference}" : ""



args = ""

if(method=="align"){
	args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev}"
}else{
	if(method=="sequential"){
		args = "--alpha ${alpha} --maxerror ${maxerror} --minlen ${minlen} --maxlen ${maxlen} ${scanrev} ${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
	}else{
		if(method=="reference"){
			args = "${ref_file} --minident ${minident} --evalue ${evalue} --maxhits ${maxhits} ${fill} --aligner ${aligner}"
		}else{
			args = "--gap ${gap}"
		}
	}
}


readArray = reads.toString().split(' ')	


if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[1]
	
	if(R1.contains("."+head_seqeunce_file)){
		R1 = readArray[0]
		R2 = readArray[1]
	}else{
		R2 = readArray[0]
		R1 = readArray[1]
	}
	
	"""
	if [ "${method}" != "align" ]; then
		if  [ "${aligner}" == "usearch" ]; then
			wget -q --show-progress --no-check-certificate https://drive5.com/downloads/usearch${usearch_version}_i86linux32.gz
			gunzip usearch${usearch_version}_i86linux32.gz
			chmod +x usearch${usearch_version}_i86linux32
			mv usearch${usearch_version}_i86linux32 /usr/local/bin/usearch2
			align_exec="--exec /usr/local/bin/usearch2"
			dbexec="--dbexec /usr/local/bin/usearch2"
		else
			align_exec="--exec /usr/local/bin/blastn"
			dbexec="--dbexec /usr/local/bin/makeblastdb"
		fi
	else
		align_exec=""
		dbexec=""
	fi

	AssemblePairs.py ${method} -1 ${R1} -2 ${R2} ${coord} ${rc} ${head_fields_R1} ${head_fields_R2} ${args} \$align_exec \$dbexec ${fasta} ${failed} --log AP_${name}.log ${nproc}  2>&1 | tee out_${R1}_AP.log
	"""

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}

}


process Mask_Primer_MaskPrimers {

input:
 val mate from g_6_mate_g1_11
 set val(name),file(reads) from g3_12_reads0_g1_11

output:
 set val(name), file("*_primers-pass.fastq")  into g1_11_reads0_g4_10
 set val(name), file("*_primers-fail.fastq") optional true  into g1_11_reads_failed11
 set val(name), file("MP_*")  into g1_11_logFile2_g1_9
 set val(name),file("out*")  into g1_11_logFile33

script:
method = params.Mask_Primer_MaskPrimers.method
barcode_field = params.Mask_Primer_MaskPrimers.barcode_field
primer_field = params.Mask_Primer_MaskPrimers.primer_field
barcode = params.Mask_Primer_MaskPrimers.barcode
revpr = params.Mask_Primer_MaskPrimers.revpr
mode = params.Mask_Primer_MaskPrimers.mode
failed = params.Mask_Primer_MaskPrimers.failed
nproc = params.Mask_Primer_MaskPrimers.nproc
maxerror = params.Mask_Primer_MaskPrimers.maxerror
umi_length = params.Mask_Primer_MaskPrimers.umi_length
start = params.Mask_Primer_MaskPrimers.start
extract_length = params.Mask_Primer_MaskPrimers.extract_length
maxlen = params.Mask_Primer_MaskPrimers.maxlen
skiprc = params.Mask_Primer_MaskPrimers.skiprc
R1_primers = params.Mask_Primer_MaskPrimers.R1_primers
R2_primers = params.Mask_Primer_MaskPrimers.R2_primers
//* @style @condition:{method="score",umi_length,start,maxerror}{method="extract",umi_length,start},{method="align",maxerror,maxlen,skiprc}, {method="extract",start,extract_length} @array:{method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc} @multicolumn:{method,barcode_field,primer_field,barcode,revpr,mode,failed,nproc,maxerror,umi_length,start,extract_length,maxlen,skiprc}

method = (method.collect().size==2) ? method : [method[0],method[0]]
barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
start = (start.collect().size==2) ? start : [start[0],start[0]]
extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
failed = (failed=="true") ? "--failed" : ""

def args_values = [];
[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
    
    if(m=="align"){
        s = ""
    }else{
        if(bc=="false"){
            s = "--start ${s}"
        }else{
            s = s + ul
            s = "--start ${s}"
        }
    }
    
    el = (m=="extract") ? "--len ${el}" : ""
    mr = (m=="extract") ? "" : "--maxerror ${mr}" 
    ml = (m=="align") ? "--maxlen ${ml}" : "" 
    sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
    
    PRIMER_FIELD = "${pf}"
    
    // all
    bf = (bf=="") ? "" : "--bf ${bf}"
    pf = (pf=="") ? "" : "--pf ${pf}"
    bc = (bc=="false") ? "" : "--barcode"
    rp = (rp=="false") ? "" : "--revpr"
    args_values.add("${m} --mode ${md} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk} ${pf} ${bf}")
    
    
}}

readArray = reads.toString().split(' ')
if(mate=="pair"){
	args_1 = args_values[0]
	args_2 = args_values[1]
	
  


	R1 = readArray[0]
	R2 = readArray[1]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
	
	
	"""
	
	MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} 2>&1 | tee -a out_${R1}_MP.log & \
	MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} 2>&1 | tee -a out_${R1}_MP.log & \
	wait
	"""
}else{
	args_1 = args_values[0]
	
	R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
	
	R1 = readArray[0]

	"""
	echo -e "Assuming inputs for R1\n"
	
	MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} 2>&1 | tee -a out_${R1}_MP.log
	"""
}

}


process Mask_Primer_parse_log_MP {

input:
 val mate from g_6_mate_g1_9
 set val(name), file(log_file) from g1_11_logFile2_g1_9

output:
 file "*table.tab"  into g1_9_logFile0_g1_12, g1_9_logFile0_g1_19

script:
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ID PRIMER BARCODE ERROR
"""

}


process Mask_Primer_try_report_maskprimer {

input:
 file primers from g1_9_logFile0_g1_12
 val mate from g_6_mate_g1_12

output:
 file "*.rmd"  into g1_12_rMarkdown0_g1_19


shell:

if(mate=="pair"){
	readArray = primers.toString().split(' ')	
	primers_1 = readArray[0]
	primers_2 = readArray[1]
	name = primers_1 - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read 1", "Read 2")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1], "(top) and", plot_titles[2], "(bottom),",
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers_1}"))
	primer_log_2 <- loadLogTable(file.path(".", "!{primers_2}"))
	
	primer_log1_error <- any(is.na(primer_log_1[['ERROR']]))
	primer_log2_error<- any(is.na(primer_log_2[['ERROR']]))
	
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	# check the error column exists 
	if(!primer_log1_error && !primer_log2_error)
		plotMaskPrimers(primer_log_1, primer_log_2, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{

	readArray = primers.toString().split(' ')
	primers = readArray[0]
	name = primers - "_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	
	plot_titles<- c("Read")
	print(plot_titles)
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("primers_count", 
	        paste("Count of assigned primers for",  plot_titles[1],
	              "The bar height indicates the total reads assigned to the given primer,
	               stacked for those under the error rate threshold (Pass) and
	               over the threshold (Fail)."))
	figures("primers_hist", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "The error rate is the percentage of mismatches between the primer sequence and the 
	               read for the best matching primer. The dotted line indicates the error threshold used."))
	figures("primers_error", 
	        paste("Distribution of primer match error rates for", plot_titles[1],
	              "broken down by assigned primer. The error rate is the percentage of mismatches between the 
	               primer sequence and the read for the best matching primer. The dotted line indicates the error
	               threshold used."))
	```
	
	```{r, echo=FALSE}
	primer_log_1 <- loadLogTable(file.path(".", "!{primers}"))
	```
	
	# Primer Identification
	
	The MaskPrimers tool supports identification of multiplexed primers and UMIs.
	Identified primer regions may be masked (with Ns) or cut to mitigate downstream
	SHM analysis artifacts due to errors in the primer region. An annotion is added to 
	each sequences that indicates the UMI and best matching primer. In the case of
	the constant region primer, the primer annotation may also be used for isotype 
	assignment.
	
	## Count of primer matches
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles,
	                style="count", sizing="figure")
	```
	
	`r figures("primers_count")`
	
	## Primer match error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="hist", sizing="figure")
	```
	
	`r figures("primers_hist")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotMaskPrimers(primer_log_1, titles=plot_titles, 
	                style="error", sizing="figure")
	```
	
	`r figures("primers_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}
}


process Mask_Primer_presto_render_rmarkdown {

input:
 file rmk from g1_12_rMarkdown0_g1_19
 file log_file from g1_9_logFile0_g1_19

output:
 file "*.html"  into g1_19_outputFileHTML00
 file "*csv" optional true  into g1_19_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process Assemble_pairs_parse_log_AP {

input:
 set val(name),file(log_file) from g3_12_logFile1_g3_15
 val mate from g_6_mate_g3_15

output:
 file "*table.tab"  into g3_15_logFile0_g3_25, g3_15_logFile0_g3_19

script:
field_to_parse = params.Assemble_pairs_parse_log_AP.field_to_parse
readArray = log_file.toString()	

"""
ParseLog.py -l ${readArray}  -f ${field_to_parse}
"""


}


process Assemble_pairs_report_assemble_pairs {

input:
 file log_files from g3_15_logFile0_g3_19
 val matee from g_6_mate_g3_19

output:
 file "*.rmd"  into g3_19_rMarkdown0_g3_25



shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')
	assemble = readArray[0]
	name = assemble-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	```{r, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("assemble_length", "Histogram showing the distribution assembled sequence lengths in 
	                            nucleotides for the Align step (top) and Reference step (bottom).")
	figures("assemble_overlap", "Histogram showing the distribution of overlapping nucleotides between 
	                             mate-pairs for the Align step (top) and Reference step (bottom).
	                             Negative values for overlap indicate non-overlapping mate-pairs
	                             with the negative value being the number of gap characters between
	                             the ends of the two mate-pairs.")
	figures("assemble_error", "Histograms showing the distribution of paired-end assembly error 
	                           rates for the Align step (top) and identity to the reference germline 
	                           for the Reference step (bottom).")
	figures("assemble_pvalue", "Histograms showing the distribution of significance scores for 
	                            paired-end assemblies. P-values for the Align mode are shown in the top
	                            panel. E-values from the Reference step's alignment against the 
	                            germline sequences are shown in the bottom panel for both input files
	                            separately.")
	```
	
	```{r, echo=FALSE, warning=FALSE}
	assemble_log <- loadLogTable(file.path(".", "!{assemble}"))
	
	# Subset to align and reference logs
	align_fields <- c("ERROR", "PVALUE")
	ref_fields <- c("REFID", "GAP", "EVALUE1", "EVALUE2", "IDENTITY")
	align_log <- assemble_log[!is.na(assemble_log$ERROR), !(names(assemble_log) %in% ref_fields)]
	ref_log <- assemble_log[!is.na(assemble_log$REFID), !(names(assemble_log) %in% align_fields)]
	
	# Build log set
	assemble_list <- list()
	if (nrow(align_log) > 0) { assemble_list[["Align"]] <- align_log }
	if (nrow(ref_log) > 0) { assemble_list[["Reference"]] <- ref_log }
	plot_titles <- names(assemble_list)
	```
	
	# Paired-End Assembly
	
	Assembly of paired-end reads is performed using the AssemblePairs tool which 
	determines the read overlap in two steps. First, de novo assembly is attempted 
	using an exhaustive approach to identify all possible overlaps between the 
	two reads with alignment error rates and p-values below user-defined thresholds. 
	This method is denoted as the `Align` method in the following figures. 
	Second, those reads failing the first stage of de novo assembly are then 
	mapped to the V-region reference sequences to create a full length sequence, 
	padding with Ns, for any amplicons that have insufficient overlap for 
	de novo assembly. This second stage is referred to as the `Reference` step in the
	figures below.
	
	## Assembled sequence lengths
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="length", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_length")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="overlap", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_overlap")`
	
	## Alignment error rates and significance
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="error", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```
	
	`r figures("assemble_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plot_params <- list(titles=plot_titles, style="pvalue", sizing="figure")
	do.call(plotAssemblePairs, c(assemble_list, plot_params))
	```

	`r figures("assemble_pvalue")`

	EOF
	
	open OUT, ">AP_!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	"""
	echo -e 'AssemblePairs works only on pair-end reads.'
	"""
}
}


process Assemble_pairs_presto_render_rmarkdown {

input:
 file rmk from g3_19_rMarkdown0_g3_25
 file log_file from g3_15_logFile0_g3_25

output:
 file "*.html"  into g3_25_outputFileHTML00
 file "*csv" optional true  into g3_25_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}

boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            mc = (mc=="none") ? "" : "--n ${mc}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


process Build_Consensus_build_consensus {

input:
 set val(name),file(reads) from g1_11_reads0_g4_10
 val mate from g_6_mate_g4_10

output:
 set val(name),file("*_consensus-pass.fastq")  into g4_10_reads0_g_16
 set val(name),file("BC*")  into g4_10_logFile1_g4_12
 set val(name),file("out*")  into g4_10_logFile22

script:
failed = params.Build_Consensus_build_consensus.failed
nproc = params.Build_Consensus_build_consensus.nproc
barcode_field = params.Build_Consensus_build_consensus.barcode_field
primer_field = params.Build_Consensus_build_consensus.primer_field
act = params.Build_Consensus_build_consensus.act
copy_field = params.Build_Consensus_build_consensus.copy_field
mincount = params.Build_Consensus_build_consensus.mincount
minqual = params.Build_Consensus_build_consensus.minqual
minfreq = params.Build_Consensus_build_consensus.minfreq
maxerror = params.Build_Consensus_build_consensus.maxerror
prcons = params.Build_Consensus_build_consensus.prcons
maxgap = params.Build_Consensus_build_consensus.maxgap
maxdiv = params.Build_Consensus_build_consensus.maxdiv
dep = params.Build_Consensus_build_consensus.dep
//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

// args 
if(isCollectionOrArray_bc(args_values_bc)){
	args_1 = args_values_bc[0]
	args_2 = args_values_bc[1]
}else{
	args_1 = args_values_bc
	args_2 = args_values_bc
}

failed = (failed=="true") ? "--failed" : "" 


if(mate=="pair"){
	// files
	readArray = reads.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	
	"""
	BuildConsensus.py --version
	BuildConsensus.py -s $R1 ${args_1} --log BC_${name}_R1.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	BuildConsensus.py -s $R2 ${args_2} --log BC_${name}_R2.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	"""
}else{
	"""
	BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc} 2>&1 | tee -a out_${R1}_BC.log
	"""
}


}


process collapse_seq {

input:
 set val(name), file(reads) from g4_10_reads0_g_16

output:
 set val(name),  file("*_collapse-unique.fast*")  into g_16_reads0_g_14
 set val(name),  file("*_collapse-duplicate.fast*") optional true  into g_16_reads_duplicate11
 set val(name),  file("*_collapse-undetermined.fast*") optional true  into g_16_reads_undetermined22
 file "CS_*"  into g_16_logFile33

script:
max_missing = params.collapse_seq.max_missing
inner = params.collapse_seq.inner
fasta = params.collapse_seq.fasta
act = params.collapse_seq.act
uf = params.collapse_seq.uf
cf = params.collapse_seq.cf
nproc = params.collapse_seq.nproc
failed = params.collapse_seq.failed

inner = (inner=="true") ? "--inner" : ""
fasta = (fasta=="true") ? "--fasta" : ""
act = (act=="none") ? "" : "--act ${act}"
cf = (cf=="") ? "" : "--cf ${cf}"
uf = (uf=="") ? "" : "--uf ${uf}"
failed = (failed=="false") ? "" : "--failed"

"""
CollapseSeq.py -s ${reads} -n ${max_missing} ${fasta} ${inner} ${uf} ${cf} ${act} --log CS_${name}.log ${failed}
"""

}


process vdjbase_input {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_16_reads0_g_14

output:
 file "${chain}"  into g_14_germlineDb00

script:
chain = params.vdjbase_input.chain

"""
mkdir ${chain}
mv ${reads} ${chain}/${name}.fasta
"""

}


process Build_Consensus_parse_log_BC {

input:
 set val(name),file(log_file) from g4_10_logFile1_g4_12
 val mate from g_6_mate_g4_12

output:
 file "*table.tab"  into g4_12_logFile0_g4_14, g4_12_logFile0_g4_20

script:
readArray = log_file.toString()

"""
ParseLog.py -l ${readArray} -f BARCODE SEQCOUNT CONSCOUNT PRCONS PRFREQ ERROR
"""

}


process Build_Consensus_report_Build_Consensus {

input:
 val matee from g_6_mate_g4_14
 file log_files from g4_12_logFile0_g4_14

output:
 file "*.rmd"  into g4_14_rMarkdown0_g4_20




shell:

if(matee=="pair"){
	readArray = log_files.toString().split(' ')	
	R1 = readArray[0]
	R2 = readArray[1]
	name = R1-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	plot_titles <- c("Read 1", "Read 2")
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", 
	        paste("Histogram of UMI read group sizes (reads per UMI) for",  
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "The x-axis indicates the number of reads in a UMI group and the y-axis is the 
	               number of UMI groups with that size. The Consensus and Total bars are overlayed 
	               (not stacked) histograms indicating whether the distribution has been calculated 
	               using the total number of reads (Total) or only those reads used for consensus 
	               generation (Consensus)."))
	figures("cons_prfreq", 
	        paste("Histograms showing the distribution of majority primer frequency for all UMI read groups for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prsize", 
	        paste("Violin plots showing the distribution of UMI read group sizes by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when running
	               BuildConsensus. Meaning, only retained UMI groups."))
	figures("cons_error", 
	        paste("Histogram showing the distribution of UMI read group error rates for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom)."))
	figures("cons_prerror", 
	        paste("Violin plots showing the distribution of UMI read group error rates by majority primer for",
	              plot_titles[1], "(top) and", plot_titles[2], "(bottom).",
	              "Only groups with majority primer frequency over the PRFREQ threshold set when 
	               running BuildConsensus. Meaning, only retained UMI groups."))
	```
	
	```{r, echo=FALSE}
	consensus_log_1 <- loadLogTable(file.path(".", "!{R1}"))
	consensus_log_2 <- loadLogTable(file.path(".", "!{R2}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group primer frequencies
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prfreq", sizing="figure")
	```
	
	`r figures("cons_prfreq")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prsize", sizing="figure")
	```
	
	`r figures("cons_prsize")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log_1, consensus_log_2, titles=plot_titles, 
	                   style="prerror", sizing="figure")
	```
	
	`r figures("cons_prerror")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''

}else{
	
	readArray = log_files.toString().split(' ')
	R1 = readArray[0]
	name = R1-"_table.tab"
	'''
	#!/usr/bin/env perl
	
	
	my $script = <<'EOF';
	
	
	
		
	```{R, message=FALSE, echo=FALSE, results="hide"}
	# Setup
	library(prestor)
	library(knitr)
	library(captioner)
	
	if (!exists("tables")) { tables <- captioner(prefix="Table") }
	if (!exists("figures")) { figures <- captioner(prefix="Figure") }
	figures("cons_size", "Histogram of UMI read group sizes (reads per UMI). 
	                      The x-axis indicates the number of reads 
	                      in a UMI group and the y-axis is the number of UMI groups 
	                      with that size. The Consensus and Total bars are overlayed
	                      (not stacked) histograms indicating whether the distribution
	                      has been calculated using the total number of reads (Total)
	                      or only those reads used for consensus generation (Consensus).")
	figures("cons_error", "Histogram showing the distribution of UMI read group error rates.")
	```
	
	```{r, echo=FALSE}
	consensus_log <- loadLogTable(file.path(".", "!{R1}"))
	```
	
	# Generation of UMI Consensus Sequences
	
	Reads sharing the same UMI are collapsed into a single consensus sequence by
	the BuildConsensus tool. BuildConsensus considers several factors in determining
	the final consensus sequence, including the number of reads in a UMI group, 
	Phred quality scores (`Q`), primer annotations, and the number of mismatches 
	within a UMI group. Quality scores are used to resolve conflicting base calls in
	a UMI read group and the final consensus sequence is assigned consensus quality 
	scores derived from the individual base quality scores. The numbers of reads in a UMI
	group, number of matching primer annotations, and error rate (average base mismatches from 
	consensus) are used as strict cut-offs for exclusion of erroneous UMI read groups.
	Additionally, individual reads are excluded whose primer annotation differs from 
	the majority in cases where there are sufficient number of reads exceeding 
	the primer consensus cut-off.
	
	## Reads per UMI
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="size", sizing="figure")
	```
	
	`r figures("cons_size")`
	
	## UMI read group error rates
	
	```{r, echo=FALSE, warning=FALSE}
	plotBuildConsensus(consensus_log, style="error", sizing="figure")
	```
	
	`r figures("cons_error")`
	
	EOF
	
	open OUT, ">!{name}.rmd";
	print OUT $script;
	close OUT;
	
	'''
}

}


process Build_Consensus_presto_render_rmarkdown {

input:
 file rmk from g4_14_rMarkdown0_g4_20
 file log_file from g4_12_logFile0_g4_20

output:
 file "*.html"  into g4_20_outputFileHTML00
 file "*csv" optional true  into g4_20_csvFile11

"""

#!/usr/bin/env Rscript 

rmarkdown::render("${rmk}", clean=TRUE, output_format="html_document", output_dir=".")

"""
}


process metadata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.json$/) "metadata/$filename"}

output:
 file "*.json"  into g_12_jsonFile00

script:
metadata = params.metadata.metadata
"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

data <- read_json("${metadata}") 

versions <- lapply(1:length(data), function(i){
	
	docker <- data[i]
	tool <- names(data)[i]
	
	if(grepl("Custom", docker)){
		ver <- "0.0"
	}else{
		ver <- system(paste0(tool," --version"), intern = TRUE)
		ver <- gsub(paste0(tool,": "), "", ver)
	}
	ver
	
})

names(versions) <- names(data)

json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        software_versions = versions
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
