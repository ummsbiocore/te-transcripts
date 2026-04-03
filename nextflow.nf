$HOSTNAME = ""
params.outdir = 'results'  

def pathChecker(input, path, type){
	def cmd = "mkdir -p check && mv ${input} check/. "
	if (!input || input.empty()){
		input = file(path).getName().toString()
		cmd = "mkdir -p check && cd check && ln -s ${path} ${input} && cd .."
		if (path.indexOf('s3:') > -1 || path.indexOf('S3:') >-1){
			def recursive = (type == "folder") ? "--recursive" : ""
			cmd = "mkdir -p check && cd check && aws s3 cp ${recursive} ${path} ${workDir}/${input} && ln -s ${workDir}/${input} . && cd .."
		} else if (path.indexOf('gs:') > -1 || path.indexOf('GS:') >-1){
			if (type == "folder"){
				cmd = "mkdir -p check ${workDir}/${input} && cd check && gsutil rsync -r ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			} else {
				cmd = "mkdir -p check && cd check && gsutil cp ${path} ${workDir}/${input} && cp -R ${workDir}/${input} . && cd .."
			}
		} else if (path.indexOf('/') == -1){
			cmd = ""
		}
}
	return [cmd,input]
}
if (!params.groups_file){params.groups_file = ""} 
if (!params.compare_file){params.compare_file = ""} 
if (!params.rmsk_TE_gtf){params.rmsk_TE_gtf = ""} 
if (!params.gtf){params.gtf = ""} 
if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.run_DESeq2){params.run_DESeq2 = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
ch_empty_file_4 = file("$baseDir/.emptyfiles/NO_FILE_4", hidden:true)

g_3_1_g20_0 = params.groups_file && file(params.groups_file, type: 'any').exists() ? file(params.groups_file, type: 'any') : ch_empty_file_1
g_3_1_g20_1 = params.groups_file && file(params.groups_file, type: 'any').exists() ? file(params.groups_file, type: 'any') : ch_empty_file_1
g_4_2_g20_0 = params.compare_file && file(params.compare_file, type: 'any').exists() ? file(params.compare_file, type: 'any') : ch_empty_file_2
g_4_0_g20_1 = params.compare_file && file(params.compare_file, type: 'any').exists() ? file(params.compare_file, type: 'any') : ch_empty_file_2
g_11_3_g_25 = file(params.rmsk_TE_gtf, type: 'any')
g_12_2_g_25 = file(params.gtf, type: 'any')
if (params.reads){
Channel
	.fromFilePairs( params.reads,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_15_1_g17_28}
 (g_15_0_g17_18) = [g_15_1_g17_28]
 } else {  
	g_15_1_g17_28 = Channel.empty()
	g_15_0_g17_18 = Channel.empty()
 }

Channel.value(params.mate).set{g_16_1_g14_30}
(g_16_0_g14_31,g_16_1_g17_11,g_16_1_g17_16,g_16_1_g17_21,g_16_1_g17_24,g_16_0_g17_28,g_16_0_g17_31,g_16_1_g17_18,g_16_1_g17_23,g_16_1_g17_19,g_16_1_g17_20) = [g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30,g_16_1_g14_30]
Channel.value(params.run_DESeq2).set{g_21_3_g20_0}

//* params.gtf =  ""  //* @input
//* params.genome =  ""  //* @input
//* params.commondb =  ""  //* @input
//* params.genome_source =  ""  //* @input
//* params.gtf_source =  ""  //* @input
//* params.commondb_source =  ""  //* @input @optional

def downFile(path, task){
	println workDir
    if (path.take(1).indexOf("/") == 0){
      target=path
      if (task.executor == "awsbatch" || task.executor == "google-batch") {
      	a=file(path)
    	fname = a.getName().toString()
    	target = "${workDir}/${fname}"
    	if (!file(target).exists()){
    		a.copyTo(workDir)
    	}
      }
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${workDir}/${fname}"
      if (!file(target).exists()){
    		a.copyTo(workDir)
      } 
    }
    return target
}

def getLastName (str){
	if (str.indexOf("/") > -1){
		return  str.substring(str.lastIndexOf('/')+1,str.length())
	} 
	return ""
}

process Check_and_Build_Module_Check_Genome_GTF {


output:
 path "${newNameFasta}"  ,emit:g13_21_genome00_g13_52 
 path "${newNameGtf}"  ,emit:g13_21_gtfFile10_g13_53 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeSource = !file("${params.genome}").exists() ? params.genome_source : params.genome
genomeName = getLastName(genomeSource)

gtfSource = !file("${params.gtf}").exists() ? params.gtf_source : params.gtf
gtfName = getLastName(gtfSource)


newNameGtf = gtfName
newNameFasta = genomeName
if (gtfName.contains('.gz')) { newNameGtf =  newNameGtf - '.gz'  } 
if (genomeName.contains('.gz')) { newNameFasta =  newNameFasta - '.gz'  } 

runGzip = ""
if (gtfName.contains('.gz') || genomeName.contains('.gz')) {
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} 

slashCountGenome = params.genome_source.count("/")
cutDirGenome = slashCountGenome - 3;

slashCountGtf = params.gtf_source.count("/")
cutDirGtf = slashCountGtf - 3;

"""
if [ ! -e "${params.genome_source}" ] ; then
    echo "${params.genome_source} not found"
	if [[ "${params.genome_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.genome_source}"
		aws s3 cp ${params.genome_source} ${workDir}/${genomeName} && ln -s ${workDir}/${genomeName} ${genomeName}
	elif [[ "${params.genome_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.genome_source}"
		gsutil cp  ${params.genome_source} ${genomeName}
	else
		echo "Downloading genome with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGenome -R 'index.html*' -r --no-parent  ${params.genome_source}
	fi

else 
	ln -s ${params.genome_source} ${genomeName}
fi

if [ ! -e "${params.gtf_source}" ] ; then
    echo "${params.gtf_source} not found"
	if [[ "${params.gtf_source}" =~ "s3" ]]; then
		echo "Downloading s3 path from ${params.gtf_source}"
		aws s3 cp  ${params.gtf_source} ${workDir}/${gtfName} && ln -s ${workDir}/${gtfName} ${gtfName}
	elif [[ "${params.gtf_source}" =~ "gs" ]]; then
		echo "Downloading gs path from ${params.gtf_source}"
		gsutil cp  ${params.gtf_source} ${gtfName}
	else
		echo "Downloading gtf with wget"
		wget --no-check-certificate --secure-protocol=TLSv1 -l inf -nc -nH --cut-dirs=$cutDirGtf -R 'index.html*' -r --no-parent  ${params.gtf_source}
	fi

else 
	ln -s ${params.gtf_source} ${gtfName}
fi

$runGzip

"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.genome_sizes =  ""  //* @input

process Check_and_Build_Module_Check_chrom_sizes_and_index {

input:
 path genome

output:
 path "${genomeName}.chrom.sizes"  ,emit:g13_52_genomeSizes02_g13_54 

when:
params.run_Download_Genomic_Sources == "yes"

script:
genomeName  = genome.baseName
genome_sizes_dir = ""
if (params.genome_sizes.indexOf('/') > -1){
	genome_sizes_dir  = params.genome_sizes.substring(0, params.genome_sizes.lastIndexOf('/')) 
}

"""
if [ ! -e "${params.genome_sizes}" ] ; then
    echo "${params.genome_sizes} not found"
    cat ${genome} | awk '\$0 ~ ">" {print c; c=0;printf substr(\$1,2,100) "\\t"; } \$0 !~ ">" {c+=length(\$0);} END { print c; }' > ${genomeName}.chrom.sizes
    ##clean first empty line
    sed -i '1{/^\$/d}' ${genomeName}.chrom.sizes
    if [ "${genome_sizes_dir}" != "" ] ; then
    	mkdir -p ${genome_sizes_dir}
		cp -n ${genomeName}.chrom.sizes ${params.genome_sizes} 
	fi
else 
	cp ${params.genome_sizes} ${genomeName}.chrom.sizes
fi

"""




}

//* params.gtf2bed_path =  ""  //* @input
//* params.bed =  ""  //* @input

process Check_and_Build_Module_Check_BED12 {

input:
 path gtf

output:
 path "${gtfName}.bed"  ,emit:g13_53_bed03_g13_54 

container "${ params.IMAGE_BASE ? "${params.IMAGE_BASE}/rnaseq:4.0" : "quay.io/ummsbiocore/rnaseq:4.0" }"

when:
params.run_Download_Genomic_Sources == "yes"

script:
gtfName  = gtf.baseName
beddir = ""
if (params.bed.indexOf('/') > -1){
	beddir  = params.bed.substring(0, params.bed.lastIndexOf('/')) 
}
"""

if [ ! -e "${params.bed}" ] ; then
    echo "${params.bed} not found"
    perl ${params.gtf2bed_path} $gtf > ${gtfName}.bed
else 
	cp -n ${params.bed} ${gtfName}.bed
fi
if [ "${beddir}" != "" ] ; then
	mkdir -p ${beddir}
	cp -n ${gtfName}.bed ${params.bed} 
fi
"""




}


process Check_and_Build_Module_check_files {

input:
 path gtf
 path genome
 path genomeSizes
 path bed

output:
 path "*/${gtf2}" ,optional:true  ,emit:g13_54_gtfFile01_g14_21 
 path "*/${genome2}" ,optional:true  ,emit:g13_54_genome10_g14_21 
 path "*/${genomeSizes2}" ,optional:true  ,emit:g13_54_genomeSizes22 
 path "*/${bed2}" ,optional:true  ,emit:g13_54_bed33 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
stageInMode 'copy'

script:
(cmd1, gtf2) = pathChecker(gtf, params.gtf, "file")
(cmd2, genome2) = pathChecker(genome, params.genome, "file")
(cmd3, genomeSizes2) = pathChecker(genomeSizes, params.genome_sizes, "file")
(cmd4, bed2) = pathChecker(bed, params.bed, "file")
"""
$cmd1
$cmd2
$cmd3
$cmd4
"""
}

build_STAR_index = params.STAR_Module_Check_Build_STAR_Index.build_STAR_index

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 150
}
//* platform
//* platform
//* autofill

process STAR_Module_Check_Build_STAR_Index {

input:
 path genome
 path gtf

output:
 path "STARIndex"  ,emit:g14_21_starIndex00_g14_26 

when:
build_STAR_index == true && ((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)

script:
star_build_parameters = params.STAR_Module_Check_Build_STAR_Index.star_build_parameters
newDirName = "STARIndex" 
"""
if [ ! -e "${params.star_index}/SA" ] ; then
    echo "STAR index not found"
    mkdir -p $newDirName 
    STAR --runMode genomeGenerate ${star_build_parameters} --genomeDir $newDirName --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf}
else 
	ln -s ${params.star_index} STARIndex
fi

"""





}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process STAR_Module_check_STAR_files {

input:
 path star

output:
 path "*/${star2}" ,optional:true  ,emit:g14_26_starIndex02_g14_31 

container 'quay.io/ummsbiocore/pipeline_base_image:1.0'
stageInMode 'copy'

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
(cmd, star2) = pathChecker(star, params.star_index, "folder")
"""
$cmd
"""
}

//* params.run_Adapter_Removal =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Adapter_Removal"
//* @style @multicolumn:{seed_mismatches, palindrome_clip_threshold, simple_clip_threshold} @condition:{Tool_for_Adapter_Removal="trimmomatic", seed_mismatches, palindrome_clip_threshold, simple_clip_threshold}, {Tool_for_Adapter_Removal="fastx_clipper", discard_non_clipped}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 5
    $MEMORY = 5
}
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_Adapter_Removal {

input:
 tuple val(name), file(reads)
 val mate

output:
 tuple val(name), file("reads/*.fastq.gz")  ,emit:g17_18_reads01_g17_31 
 path "*.{fastx,trimmomatic}.log"  ,emit:g17_18_log_file10_g17_11 

errorStrategy 'retry'

when:
(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
phred = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.phred
Tool_for_Adapter_Removal = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Tool_for_Adapter_Removal
Adapter_Sequence = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
min_length = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length
seed_mismatches = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches
palindrome_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold
simple_clip_threshold = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
discard_non_clipped = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] 
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1] }
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads !{task.cpus} -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz reads/!{name}.2.fastq.gz unpaired/!{name}.2.fastq.unpaired.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper is not suitable for paired reads.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads !{task.cpus}  -phred${quality} !{file1} reads/!{name}.fastq.gz ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq.gz > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: inputs reads will be removed if they are located in the $workdir $inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary {

input:
 path logfile
 val mate

output:
 path "adapter_removal_summary.tsv"  ,emit:g17_11_outputFileTSV00 
 path "adapter_removal_detailed_summary.tsv" ,optional:true  ,emit:g17_11_outputFile11 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* @style @condition:{single_or_paired_end_reads="single", barcode_pattern1,remove_duplicates_based_on_UMI}, {single_or_paired_end_reads="pair", barcode_pattern1,barcode_pattern2}


process Adapter_Trimmer_Quality_Module_UMIextract {

input:
 tuple val(name), file(reads)
 val mate

output:
 tuple val(name), file("result/*.fastq.gz")  ,emit:g17_23_reads00_g17_19 
 path "${name}.*.log"  ,emit:g17_23_log_file10_g17_24 

container 'quay.io/ummsbiocore/fastq_preprocessing:1.0'

when:
params.run_UMIextract == "yes" 

script:
readArray = reads.toString().split(' ')
file2 = ""
file1 =  readArray[0]
if (mate == "pair") {file2 =  readArray[1]}


single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_UMIextract.single_or_paired_end_reads
barcode_pattern1 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern1
barcode_pattern2 = params.Adapter_Trimmer_Quality_Module_UMIextract.barcode_pattern2
UMIqualityFilterThreshold = params.Adapter_Trimmer_Quality_Module_UMIextract.UMIqualityFilterThreshold
phred = params.Adapter_Trimmer_Quality_Module_UMIextract.phred
remove_duplicates_based_on_UMI = params.Adapter_Trimmer_Quality_Module_UMIextract.remove_duplicates_based_on_UMI

"""
set +e
source activate umi_tools_env 2> /dev/null || true
mkdir result
if [ "${mate}" == "pair" ]; then
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --bc-pattern2='${barcode_pattern2}' \
                  --extract-method=regex \
                  --stdin=${file1} \
                  --stdout=result/${name}_R1.fastq.gz \
                  --read2-in=${file2} \
                  --read2-out=result/${name}_R2.fastq.gz\
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred} \
				  --log=${name}.umitools.log 


else
umi_tools extract --bc-pattern='${barcode_pattern1}' \
                  --log=${name}.umitools.log \
                  --extract-method=regex \
                  --stdin ${file1} \
                  --stdout result/${name}.fastq.gz \
				  --quality-filter-threshold=${UMIqualityFilterThreshold} \
				  --quality-encoding=phred${phred}
	if [ "${remove_duplicates_based_on_UMI}" == "true" ]; then		  
        mv result/${name}.fastq.gz  result/${name}_umitools.fastq.gz && gunzip result/${name}_umitools.fastq.gz
        ## only checks last part of the underscore splitted header for UMI
        awk '(NR%4==1){name=\$1;header=\$0;len=split(name,umiAr,"_");umi=umiAr[len];} (NR%4==2){total++;if(a[umi]!=1){nondup++;a[umi]=1;  print header;print;getline; print; getline; print;}} END{print FILENAME"\\t"total"\\t"nondup > "${name}.dedup.log"}' result/${name}_umitools.fastq > result/${name}.fastq
        rm result/${name}_umitools.fastq
        gzip result/${name}.fastq
	fi			  
fi
"""

}

//* params.run_Trimmer =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Trimmer"
//* @style @multicolumn:{trim_length_5prime,trim_length_3prime}, {trim_length_5prime_R1,trim_length_3prime_R1}, {trim_length_5prime_R2,trim_length_3prime_R2} @condition:{single_or_paired_end_reads="single", trim_length_5prime,trim_length_3prime}, {single_or_paired_end_reads="pair", trim_length_5prime_R1,trim_length_3prime_R1,trim_length_5prime_R2,trim_length_3prime_R2}

//* autofill
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_Trimmer {

input:
 tuple val(name), file(reads)
 val mate

output:
 tuple val(name), file("reads/*q.gz")  ,emit:g17_19_reads00_g17_20 
 path "*.log" ,optional:true  ,emit:g17_19_log_file10_g17_21 

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Adapter_Trimmer_Quality_Module_Trimmer.phred
single_or_paired_end_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime
trim_length_3prime = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Trimmer.remove_previous_reads



file1 =  reads[0] 
file2 = ""
if (mate == "pair") {file2 =  reads[1] }
rawFile1 = "_length_check1.fastq"
rawFile2 = "_length_check2.fastq"
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
runCmd("mkdir reads");
runCmd("zcat !{file1} | head -n 100 > !{rawFile1}");
if ("!{mate}" eq "pair") {
	runCmd("zcat !{file2} | head -n 100 > !{rawFile2}");
}
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength("!{rawFile1}");
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength("!{rawFile2}");
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength("!{rawFile1}");
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="gunzip -c $file | fastx_trimmer $quality -v $param -z -o $outfile  > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer skipped for $file \\n";
          runCmd("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer executed for $file \\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}


process Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary {

input:
 path logfile
 val mate

output:
 path "trimmer_summary.tsv"  ,emit:g17_21_outputFileTSV00 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.run_Quality_Filtering =   "no"   //* @dropdown @options:"yes","no" @show_settings:"Quality_Filtering"
//* @style @multicolumn:{window_size,required_quality}, {leading,trailing,minlen}, {minQuality,minPercent} @condition:{tool="trimmomatic", minlen, trailing, leading, required_quality_for_window_trimming, window_size}, {tool="fastx", minQuality, minPercent}

//* autofill
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_Quality_Filtering {

input:
 tuple val(name), file(reads)
 val mate

output:
 tuple val(name), file("reads/*.gz")  ,emit:g17_20_reads01_g14_31 
 path "*.{fastx,trimmomatic}_quality.log" ,optional:true  ,emit:g17_20_log_file10_g17_16 

when:
(params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering    

shell:
tool = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool
phred = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.phred
window_size = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size
required_quality_for_window_trimming = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming
leading = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading
trailing = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing
minlen = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen


// fastx parameters
minQuality = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality
minPercent = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent

remove_previous_reads = params.Adapter_Trimmer_Quality_Module_Quality_Filtering.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 ="";
if (nameAll.contains('.gz')) {
    file1 =  nameArray[0] 
    if (mate == "pair") {file2 =  nameArray[1]}
} 
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
runCmd("mkdir reads unpaired");
my $param = "SLIDINGWINDOW:"."!{window_size}".":"."!{required_quality_for_window_trimming}";
$param.=" LEADING:"."!{leading}";
$param.=" TRAILING:"."!{trailing}";
$param.=" MINLEN:"."!{minlen}";

my $quality="!{phred}";

print "INFO: fastq quality: $quality\\n";
     
if ("!{tool}" eq "trimmomatic") {
    if ("!{mate}" eq "pair") {
        runCmd("trimmomatic PE -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq.gz unpaired/!{name}.1.fastq.unpaired.gz reads/!{name}.2.fastq.gz unpaired/!{name}.2.fastq.unpaired.gz $param 2> !{name}.trimmomatic_quality.log");
    } else {
        runCmd("trimmomatic SE -phred${quality} !{file1} reads/!{name}.fastq.gz $param 2> !{name}.trimmomatic_quality.log");
    }
} elsif ("!{tool}" eq "fastx") {
    if ("!{mate}" eq "pair") {
        print("WARNING: Fastx option is not suitable for paired reads. This step will be skipped.");
        runCmd("mv !{file1} !{file2} reads/.");
    } else {
        runCmd("fastq_quality_filter  -Q $quality -q !{minQuality} -p !{minPercent} -v -i !{file1} -o reads/!{name}.fastq.gz > !{name}.fastx_quality.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: inputs reads will be removed if they are located in the workdir inputsdir\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            runCmd("rm -f $targetFile");
            print "INFO: $targetFile deleted.\\n";
        }
    }
}

##Subroutines
sub runCmd {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}


'''

}


process Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary {

input:
 path logfile
 val mate

output:
 path "quality_filter_summary.tsv"  ,emit:g17_16_outputFileTSV00 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i = 0;
chomp( my $contents = `ls *.log` );
my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapper   = "";
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_quality\\.log/){
        $mapper   = "fastx";
        $file =~ /(.*)\\.fastx_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    } elsif ($file =~ /(.*)\\.trimmomatic_quality\\.log/){
        $mapper   = "trimmomatic";
        $file =~ /(.*)\\.trimmomatic_quality\\.log/;
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Quality Filtering" ];
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "quality_filter_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

//* params.star_index =  ""  //* @input

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 50
}
//* platform
//* platform
//* autofill

process STAR_Module_Map_STAR {

input:
 val mate
 tuple val(name), file(reads)
 path star_index

output:
 tuple val(name), file("${name}Log.final.out")  ,emit:g14_31_outputFileOut00_g14_18 
 tuple val(name), file("${name}.flagstat.txt")  ,emit:g14_31_outputFileTxt11 
 tuple val(name), file("${name}Log.out")  ,emit:g14_31_logOut22 
 tuple val(name), file("${name}.bam")  ,emit:g14_31_mapped_reads30_g14_30 
 tuple val(name), file("${name}SJ.out.tab")  ,emit:g14_31_outputFileTab44 
 tuple val(name), file("${name}Log.progress.out")  ,emit:g14_31_progressOut55 
 tuple val(name), file("${name}Aligned.toTranscriptome.out.bam") ,optional:true  ,emit:g14_31_transcriptome_bam60_g14_15 

when:
(params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR

script:
params_STAR = params.STAR_Module_Map_STAR.params_STAR
sense_antisense = params.STAR_Module_Map_STAR.sense_antisense
transcriptomeSAM = ""
if (params.run_Salmon_after_STAR && params.run_Salmon_after_STAR == "yes" && params_STAR.indexOf("--quantMode") < 0){
	transcriptomeSAM = " --quantMode TranscriptomeSAM "
}

"""
STAR ${params_STAR} ${transcriptomeSAM} --genomeDir ${star_index} --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix ${name}
echo "Alignment completed."
if [ ! -e "${name}Aligned.toTranscriptome.out.bam" -a -e "${name}Aligned.toTranscriptome.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.toTranscriptome.out.sam > ${name}Aligned.toTranscriptome.out.bam
elif [ ! -e "${name}Aligned.out.bam" -a -e "${name}Aligned.out.sam" ] ; then
    samtools view -S -b ${name}Aligned.out.sam > ${name}Aligned.out.bam
fi
rm -rf *.sam
if [ -e "${name}Aligned.sortedByCoord.out.bam" ] ; then
    mv ${name}Aligned.sortedByCoord.out.bam ${name}.bam
elif [ -e "${name}Aligned.out.bam" ] ; then
    mv ${name}Aligned.out.bam ${name}.bam
fi

samtools flagstat ${name}.bam > ${name}.flagstat.txt
"""


}


process STAR_Module_STAR_Summary {

input:
 tuple val(name), file(alignSum)

output:
 path "*.tsv"  ,emit:g14_18_outputFile00_g14_11 
 val "star_alignment_sum"  ,emit:g14_18_name11_g14_11 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = '!{name}';


alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_star_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned (STAR)");
	push(@headers, "Unique Reads Aligned (STAR)");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($inputCount = `cat $file | grep 'Number of input reads' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($aligned = `cat $file | grep 'Uniquely mapped reads number' | awk '{sum+=\\$6} END {print sum}'`);
		chomp($multimapped = `cat $file | grep 'Number of reads mapped to multiple loci' | awk '{sum+=\\$9} END {print sum}'`);
		$multimappedSum += int($multimapped);
        $alignedSum += int($aligned);
        $inputCountSum += int($inputCount);
	}
	$tsv{$name} = [$name, $inputCountSum];
	push(@{$tsv{$name}}, $multimappedSum);
	push(@{$tsv{$name}}, $alignedSum);
}

sub runCommand {
    my ($com) = @_;
    if ($com eq ""){
		return "";
    }
    my $error = system(@_);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}


process STAR_Module_merge_tsv_files_with_same_header {

input:
 path tsv
 val outputFileName

output:
 path "${name}.tsv"  ,emit:g14_11_outputFileTSV00 

errorStrategy 'retry'
maxRetries 3

script:
name = outputFileName[0]
"""    
awk 'FNR==1 && NR!=1 {  getline; } 1 {print} ' *.tsv > ${name}.tsv
"""
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process STAR_Module_Merge_Bam_and_create_sense_antisense {

input:
 tuple val(oldname), file(bamfiles)
 val mate

output:
 path "*_sorted.bam.bai"  ,emit:g14_30_bam_index00_g_25 
 tuple val(oldname), file("*_sorted.bam")  ,emit:g14_30_bamFile11_g_25 

errorStrategy 'retry'
maxRetries 2

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi


'''
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32 
}
//* platform
//* platform
//* autofill

process TE_Count {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_counts.tsv$/) "TE_Count/$filename"}
input:
 path index
 tuple val(name), file(bam)
 path gtf
 path rmsk_TE_gtf

output:
 path "${name}_counts.tsv"  ,emit:g_25_outFileTSV00_g_23 

container 'quay.io/ummsbiocore/te-transcripts:1.0'

script:

stranded = params.TE_Count.stranded
mode = params.TE_Count.mode
iteration = params.TE_Count.iteration

//* @style @multicolumn:{stranded,mode,iteration}


"""
TEcount --BAM ${bam} \
        --GTF ${gtf} \
        --TE ${rmsk_TE_gtf} \
        --project ${name} \
        --sortByPos \
        --stranded ${stranded} \
        --mode ${mode} \
        --iteration ${iteration}
        
clean_TE_count.py --input ${name}.cntTable --output ${name}_counts.tsv
"""

}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 8
}
//* platform
//* platform
//* autofill

process TE_Summary {

input:
 path counts

output:
 path "counts.tsv"  ,emit:g_23_outputFileTSV00_g20_0 

container 'quay.io/ummsbiocore/te-transcripts:1.0'

script:
"""
collect_TE_counts.py --input-files ${counts} --output counts.tsv
"""
}


//* autofill
//* platform
//* platform
//* autofill

process STAR_Module_merge_transcriptome_bam {

input:
 tuple val(oldname), file(bamfiles)

output:
 tuple val(oldname), file("${oldname}.bam")  ,emit:g14_15_merged_bams00 
 tuple val(oldname), file("*_sorted*bai")  ,emit:g14_15_bam_index11 
 tuple val(oldname), file("*_sorted*bam")  ,emit:g14_15_sorted_bam22 

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
    samtools merge !{oldname}.bam !{bamfiles.join(" ")} && samtools sort -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
else
    mv !{bamfiles.join(" ")} !{oldname}.bam 2>/dev/null || true
    samtools sort  -o !{oldname}_sorted.bam !{oldname}.bam && samtools index !{oldname}_sorted.bam
fi
'''
}


process Adapter_Trimmer_Quality_Module_Umitools_Summary {

input:
 path logfile
 val mate

output:
 path "umitools_summary.tsv"  ,emit:g17_24_outputFileTSV00 

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.umitools\\.log/){
        $file =~ /(.*)\\.umitools\\.log/;
        my $mapper   = "umitools";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        my $dedupout;
        chomp( $in =`cat $file | grep 'INFO Input Reads:' | awk '{sum=\\$6} END {print sum}'` );
        chomp( $out =`cat $file | grep 'INFO Reads output:' | awk '{sum=\\$6} END {print sum}'` );
        my $deduplog = $name.".dedup.log";
        $headerHash{$mapOrder} = $mapper;
        if (-e $deduplog) {
            print "dedup log found\\n";
            chomp( $dedupout =`cat $deduplog | grep '$name' | awk '{sum=\\$3} END {print sum}'` );
            $tsv{$name}{$mapper} = [ $in, $out, $dedupout];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract", "Reads After Deduplication" ]; 
        } else {
            $tsv{$name}{$mapper} = [ $in, $out ];
            $headerText{$mapOrder} = [ "Total Reads", "Reads After Umiextract" ]; 
        }
        
        
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "umitools_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}


//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_FastQC {

input:
 val mate
 tuple val(name), file(reads)

output:
 path '*.{html,zip}'  ,emit:g17_28_FastQCout00 

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
"""
fastqc ${reads} 
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 10
}
//* platform
//* platform
//* autofill

process Adapter_Trimmer_Quality_Module_FastQC_after_Adapter_Removal {

input:
 val mate
 tuple val(name), file(reads)

output:
 path '*.{html,zip}'  ,emit:g17_31_FastQCout00 

when:
(params.run_FastQC && params.run_FastQC == "yes" && params.run_Adapter_Removal && params.run_Adapter_Removal == "yes")

script:
"""
fastqc ${reads} 
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 1 
}
//* platform
//* platform
//* autofill

process DESeq2_module_Prepare_DESeq2 {

input:
 path counts
 path groups_file
 path compare_file
 val run_DESeq2

output:
 path "*Rmd"  ,emit:g20_0_rMarkdown03_g20_1 

container 'quay.io/ummsbiocore/de_module:3.0'

when:
run_DESeq2 == 'yes' && groups_file && !groups_file.name.startsWith('NO_FILE') && compare_file && !compare_file.name.startsWith('NO_FILE')

script:

include_distribution = params.DESeq2_module_Prepare_DESeq2.include_distribution
include_all2all = params.DESeq2_module_Prepare_DESeq2.include_all2all
include_pca = params.DESeq2_module_Prepare_DESeq2.include_pca

filter_type = params.DESeq2_module_Prepare_DESeq2.filter_type
min_count = params.DESeq2_module_Prepare_DESeq2.min_count
min_samples = params.DESeq2_module_Prepare_DESeq2.min_samples
min_counts_per_sample = params.DESeq2_module_Prepare_DESeq2.min_counts_per_sample
excluded_events = params.DESeq2_module_Prepare_DESeq2.excluded_events

include_batch_correction = params.DESeq2_module_Prepare_DESeq2.include_batch_correction
batch_correction_column = params.DESeq2_module_Prepare_DESeq2.batch_correction_column
batch_correction_group_column = params.DESeq2_module_Prepare_DESeq2.batch_correction_group_column
batch_normalization_algorithm = params.DESeq2_module_Prepare_DESeq2.batch_normalization_algorithm

transformation = params.DESeq2_module_Prepare_DESeq2.transformation
pca_color = params.DESeq2_module_Prepare_DESeq2.pca_color
pca_shape = params.DESeq2_module_Prepare_DESeq2.pca_shape
pca_fill = params.DESeq2_module_Prepare_DESeq2.pca_fill
pca_transparency = params.DESeq2_module_Prepare_DESeq2.pca_transparency
pca_label = params.DESeq2_module_Prepare_DESeq2.pca_label

include_deseq2 = params.DESeq2_module_Prepare_DESeq2.include_deseq2
design = params.DESeq2_module_Prepare_DESeq2.design
fitType = params.DESeq2_module_Prepare_DESeq2.fitType
use_batch_corrected_in_DE = params.DESeq2_module_Prepare_DESeq2.use_batch_corrected_in_DE
apply_shrinkage = params.DESeq2_module_Prepare_DESeq2.apply_shrinkage
shrinkage_type = params.DESeq2_module_Prepare_DESeq2.shrinkage_type
include_volcano = params.DESeq2_module_Prepare_DESeq2.include_volcano
include_ma = params.DESeq2_module_Prepare_DESeq2.include_ma
include_heatmap = params.DESeq2_module_Prepare_DESeq2.include_heatmap

padj_significance_cutoff = params.DESeq2_module_Prepare_DESeq2.padj_significance_cutoff
fc_significance_cutoff = params.DESeq2_module_Prepare_DESeq2.fc_significance_cutoff
padj_floor = params.DESeq2_module_Prepare_DESeq2.padj_floor
fc_ceiling = params.DESeq2_module_Prepare_DESeq2.fc_ceiling

convert_names = params.DESeq2_module_Prepare_DESeq2.convert_names
count_file_names = params.DESeq2_module_Prepare_DESeq2.count_file_names
converted_name = params.DESeq2_module_Prepare_DESeq2.converted_name
num_labeled = params.DESeq2_module_Prepare_DESeq2.num_labeled
highlighted_genes = params.DESeq2_module_Prepare_DESeq2.highlighted_genes
include_volcano_highlighted = params.DESeq2_module_Prepare_DESeq2.include_volcano_highlighted
include_ma_highlighted = params.DESeq2_module_Prepare_DESeq2.include_ma_highlighted

//* @style @condition:{convert_names="true", count_file_names, converted_name},{convert_names="false"},{include_batch_correction="true", batch_correction_column, batch_correction_group_column, batch_normalization_algorithm, use_batch_corrected_in_DE},{include_batch_correction="false"},{include_deseq2="true", design, fitType, apply_shrinkage, shrinkage_type, include_volcano, include_ma, include_heatmap, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, convert_names, count_file_names, converted_name, num_labeled, highlighted_genes},{include_deseq2="false"},{apply_shrinkage="true", shrinkage_type},{apply_shrinkage="false"},{include_pca="true", transformation, pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{include_pca="false"} @multicolumn:{include_distribution, include_all2all, include_pca},{filter_type, min_count, min_samples, min_counts_per_sample},{include_batch_correction, batch_correction_column, batch_correction_group_column, batch_normalization_algorithm},{design, fitType, use_batch_corrected_in_DE, apply_shrinkage, shrinkage_type},{pca_color, pca_shape, pca_fill, pca_transparency, pca_label},{padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling},{convert_names, count_file_names, converted_name},{include_volcano, include_ma, include_heatmap}{include_volcano_highlighted,include_ma_highlighted} 

include_distribution = include_distribution == 'true' ? 'TRUE' : 'FALSE'
include_all2all = include_all2all == 'true' ? 'TRUE' : 'FALSE'
include_pca = include_pca == 'true' ? 'TRUE' : 'FALSE'
include_batch_correction = include_batch_correction == 'true' ? 'TRUE' : 'FALSE'
include_deseq2 = include_deseq2 == 'true' ? 'TRUE' : 'FALSE'
use_batch_corrected_in_DE = use_batch_corrected_in_DE  == 'true' ? 'TRUE' : 'FALSE'
apply_shrinkage = apply_shrinkage == 'true' ? 'TRUE' : 'FALSE'
convert_names = convert_names == 'true' ? 'TRUE' : 'FALSE'
include_ma = include_ma == 'true' ? 'TRUE' : 'FALSE'
include_volcano = include_volcano == 'true' ? 'TRUE' : 'FALSE'
include_heatmap = include_heatmap == 'true' ? 'TRUE' : 'FALSE'

excluded_events = excluded_events.replace("\n", " ").replace(',', ' ')

pca_color_arg = pca_color.equals('') ? '' : '--pca-color ' + pca_color
pca_shape_arg = pca_shape.equals('') ? '' : '--pca-shape ' + pca_shape
pca_fill_arg = pca_fill.equals('') ? '' : '--pca-fill ' + pca_fill
pca_transparency_arg = pca_transparency.equals('') ? '' : '--pca-transparency ' + pca_transparency
pca_label_arg = pca_label.equals('') ? '' : '--pca-label ' + pca_label

highlighted_genes = highlighted_genes.replace("\n", " ").replace(',', ' ')

excluded_events_arg = excluded_events.equals('') ? '' : '--excluded-events ' + excluded_events
highlighted_genes_arg = highlighted_genes.equals('') ? '' : '--highlighted-genes ' + highlighted_genes
include_ma_highlighted = include_ma_highlighted == 'true' ? 'TRUE' : 'FALSE'
include_volcano_highlighted = include_volcano_highlighted == 'true' ? 'TRUE' : 'FALSE'

"""
prepare_DESeq2.py \
--counts ${counts} --groups ${groups_file} --comparisons ${compare_file} --prefix inputs/ \
--include-distribution ${include_distribution} --include-all2all ${include_all2all} --include-pca ${include_pca} \
--filter-type ${filter_type} --min-counts-per-event ${min_count} --min-samples-per-event ${min_samples} --min-counts-per-sample ${min_counts_per_sample} \
--transformation ${transformation} ${pca_color_arg} ${pca_shape_arg} ${pca_fill_arg} ${pca_transparency_arg} ${pca_label_arg} \
--include-batch-correction ${include_batch_correction} --batch-correction-column ${batch_correction_column} --batch-correction-group-column ${batch_correction_group_column} --batch-normalization-algorithm ${batch_normalization_algorithm} \
--include-DESeq2 ${include_deseq2} --design '${design}' --fitType ${fitType} --use-batch-correction-in-DE ${use_batch_corrected_in_DE} --apply-shrinkage ${apply_shrinkage} --shrinkage-type ${shrinkage_type} \
--include-volcano ${include_volcano} --include-ma ${include_ma} --include-heatmap ${include_heatmap} \
--padj-significance-cutoff ${padj_significance_cutoff} --fc-significance-cutoff ${fc_significance_cutoff} --padj-floor ${padj_floor} --fc-ceiling ${fc_ceiling} \
--convert-names ${convert_names} --count-file-names ${count_file_names} --converted-names ${converted_name} --num-labeled ${num_labeled} \
${highlighted_genes_arg} --include-volcano-highlighted ${include_volcano_highlighted} --include-ma-highlighted ${include_ma_highlighted} \
${excluded_events_arg}
"""

}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 30 
}
//* platform
//* platform
//* autofill

process DESeq2_module_run_DESeq2_markdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /(.*.html|.*.Rmd|inputs|outputs)$/) "DESeq2/$filename"}
input:
 path comparison_file
 path group_file
 path count_file
 path rmd

output:
 path "{*.html,*.Rmd,inputs,outputs}"  ,emit:g20_1_rMarkdown00 

container 'quay.io/ummsbiocore/de_module:3.0'

script:

basename=rmd.baseName
finalname = basename + "_DE.Rmd"

"""
mkdir inputs
mkdir outputs

mv ${rmd} ${finalname}
mv ${count_file} inputs/${count_file}
mv ${group_file} inputs/${group_file}
mv ${comparison_file} inputs/${comparison_file}
cp inputs/${group_file} inputs/debrowser_metadata.txt
Rscript -e 'rmarkdown::render("${finalname}", "html_document")'
"""
}


workflow {


Check_and_Build_Module_Check_Genome_GTF()
g13_21_genome00_g13_52 = Check_and_Build_Module_Check_Genome_GTF.out.g13_21_genome00_g13_52
(g13_21_genome01_g13_54) = [g13_21_genome00_g13_52]
g13_21_gtfFile10_g13_53 = Check_and_Build_Module_Check_Genome_GTF.out.g13_21_gtfFile10_g13_53
(g13_21_gtfFile10_g13_54) = [g13_21_gtfFile10_g13_53]


Check_and_Build_Module_Check_chrom_sizes_and_index(g13_21_genome00_g13_52)
g13_52_genomeSizes02_g13_54 = Check_and_Build_Module_Check_chrom_sizes_and_index.out.g13_52_genomeSizes02_g13_54


Check_and_Build_Module_Check_BED12(g13_21_gtfFile10_g13_53)
g13_53_bed03_g13_54 = Check_and_Build_Module_Check_BED12.out.g13_53_bed03_g13_54

g13_21_gtfFile10_g13_54= g13_21_gtfFile10_g13_54.ifEmpty(ch_empty_file_1) 
g13_21_genome01_g13_54= g13_21_genome01_g13_54.ifEmpty(ch_empty_file_2) 
g13_52_genomeSizes02_g13_54= g13_52_genomeSizes02_g13_54.ifEmpty(ch_empty_file_3) 
g13_53_bed03_g13_54= g13_53_bed03_g13_54.ifEmpty(ch_empty_file_4) 


Check_and_Build_Module_check_files(g13_21_gtfFile10_g13_54,g13_21_genome01_g13_54,g13_52_genomeSizes02_g13_54,g13_53_bed03_g13_54)
g13_54_gtfFile01_g14_21 = Check_and_Build_Module_check_files.out.g13_54_gtfFile01_g14_21
g13_54_genome10_g14_21 = Check_and_Build_Module_check_files.out.g13_54_genome10_g14_21
g13_54_genomeSizes22 = Check_and_Build_Module_check_files.out.g13_54_genomeSizes22
g13_54_bed33 = Check_and_Build_Module_check_files.out.g13_54_bed33


STAR_Module_Check_Build_STAR_Index(g13_54_genome10_g14_21,g13_54_gtfFile01_g14_21)
g14_21_starIndex00_g14_26 = STAR_Module_Check_Build_STAR_Index.out.g14_21_starIndex00_g14_26

g14_21_starIndex00_g14_26= g14_21_starIndex00_g14_26.ifEmpty(ch_empty_file_1) 


if (!((params.run_STAR && (params.run_STAR == "yes")) || !params.run_STAR)){
g14_21_starIndex00_g14_26.set{g14_26_starIndex02_g14_31}
} else {

STAR_Module_check_STAR_files(g14_21_starIndex00_g14_26)
g14_26_starIndex02_g14_31 = STAR_Module_check_STAR_files.out.g14_26_starIndex02_g14_31
}


if (!((params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal)){
g_15_0_g17_18.set{g17_18_reads01_g17_31}
(g17_18_reads00_g17_23) = [g17_18_reads01_g17_31]
g17_18_log_file10_g17_11 = Channel.empty()
} else {

Adapter_Trimmer_Quality_Module_Adapter_Removal(g_15_0_g17_18,g_16_1_g17_18)
g17_18_reads01_g17_31 = Adapter_Trimmer_Quality_Module_Adapter_Removal.out.g17_18_reads01_g17_31
(g17_18_reads00_g17_23) = [g17_18_reads01_g17_31]
g17_18_log_file10_g17_11 = Adapter_Trimmer_Quality_Module_Adapter_Removal.out.g17_18_log_file10_g17_11
}


Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary(g17_18_log_file10_g17_11.collect(),g_16_1_g17_11)
g17_11_outputFileTSV00 = Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary.out.g17_11_outputFileTSV00
g17_11_outputFile11 = Adapter_Trimmer_Quality_Module_Adapter_Removal_Summary.out.g17_11_outputFile11


if (!(params.run_UMIextract == "yes")){
g17_18_reads00_g17_23.set{g17_23_reads00_g17_19}
g17_23_log_file10_g17_24 = Channel.empty()
} else {

Adapter_Trimmer_Quality_Module_UMIextract(g17_18_reads00_g17_23,g_16_1_g17_23)
g17_23_reads00_g17_19 = Adapter_Trimmer_Quality_Module_UMIextract.out.g17_23_reads00_g17_19
g17_23_log_file10_g17_24 = Adapter_Trimmer_Quality_Module_UMIextract.out.g17_23_log_file10_g17_24
}


if (!((params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer)){
g17_23_reads00_g17_19.set{g17_19_reads00_g17_20}
g17_19_log_file10_g17_21 = Channel.empty()
} else {

Adapter_Trimmer_Quality_Module_Trimmer(g17_23_reads00_g17_19,g_16_1_g17_19)
g17_19_reads00_g17_20 = Adapter_Trimmer_Quality_Module_Trimmer.out.g17_19_reads00_g17_20
g17_19_log_file10_g17_21 = Adapter_Trimmer_Quality_Module_Trimmer.out.g17_19_log_file10_g17_21
}


Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary(g17_19_log_file10_g17_21.collect(),g_16_1_g17_21)
g17_21_outputFileTSV00 = Adapter_Trimmer_Quality_Module_Trimmer_Removal_Summary.out.g17_21_outputFileTSV00


if (!((params.run_Quality_Filtering && (params.run_Quality_Filtering == "yes")) || !params.run_Quality_Filtering)){
g17_19_reads00_g17_20.set{g17_20_reads01_g14_31}
g17_20_log_file10_g17_16 = Channel.empty()
} else {

Adapter_Trimmer_Quality_Module_Quality_Filtering(g17_19_reads00_g17_20,g_16_1_g17_20)
g17_20_reads01_g14_31 = Adapter_Trimmer_Quality_Module_Quality_Filtering.out.g17_20_reads01_g14_31
g17_20_log_file10_g17_16 = Adapter_Trimmer_Quality_Module_Quality_Filtering.out.g17_20_log_file10_g17_16
}


Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary(g17_20_log_file10_g17_16.collect(),g_16_1_g17_16)
g17_16_outputFileTSV00 = Adapter_Trimmer_Quality_Module_Quality_Filtering_Summary.out.g17_16_outputFileTSV00


STAR_Module_Map_STAR(g_16_0_g14_31,g17_20_reads01_g14_31,g14_26_starIndex02_g14_31)
g14_31_outputFileOut00_g14_18 = STAR_Module_Map_STAR.out.g14_31_outputFileOut00_g14_18
g14_31_outputFileTxt11 = STAR_Module_Map_STAR.out.g14_31_outputFileTxt11
g14_31_logOut22 = STAR_Module_Map_STAR.out.g14_31_logOut22
g14_31_mapped_reads30_g14_30 = STAR_Module_Map_STAR.out.g14_31_mapped_reads30_g14_30
g14_31_outputFileTab44 = STAR_Module_Map_STAR.out.g14_31_outputFileTab44
g14_31_progressOut55 = STAR_Module_Map_STAR.out.g14_31_progressOut55
g14_31_transcriptome_bam60_g14_15 = STAR_Module_Map_STAR.out.g14_31_transcriptome_bam60_g14_15


STAR_Module_STAR_Summary(g14_31_outputFileOut00_g14_18.groupTuple())
g14_18_outputFile00_g14_11 = STAR_Module_STAR_Summary.out.g14_18_outputFile00_g14_11
g14_18_name11_g14_11 = STAR_Module_STAR_Summary.out.g14_18_name11_g14_11


STAR_Module_merge_tsv_files_with_same_header(g14_18_outputFile00_g14_11.collect(),g14_18_name11_g14_11.collect())
g14_11_outputFileTSV00 = STAR_Module_merge_tsv_files_with_same_header.out.g14_11_outputFileTSV00


STAR_Module_Merge_Bam_and_create_sense_antisense(g14_31_mapped_reads30_g14_30.groupTuple(),g_16_1_g14_30)
g14_30_bam_index00_g_25 = STAR_Module_Merge_Bam_and_create_sense_antisense.out.g14_30_bam_index00_g_25
g14_30_bamFile11_g_25 = STAR_Module_Merge_Bam_and_create_sense_antisense.out.g14_30_bamFile11_g_25


TE_Count(g14_30_bam_index00_g_25,g14_30_bamFile11_g_25,g_12_2_g_25,g_11_3_g_25)
g_25_outFileTSV00_g_23 = TE_Count.out.g_25_outFileTSV00_g_23


TE_Summary(g_25_outFileTSV00_g_23.collect())
g_23_outputFileTSV00_g20_0 = TE_Summary.out.g_23_outputFileTSV00_g20_0
(g_23_outputFileTSV02_g20_1) = [g_23_outputFileTSV00_g20_0]


STAR_Module_merge_transcriptome_bam(g14_31_transcriptome_bam60_g14_15.groupTuple())
g14_15_merged_bams00 = STAR_Module_merge_transcriptome_bam.out.g14_15_merged_bams00
g14_15_bam_index11 = STAR_Module_merge_transcriptome_bam.out.g14_15_bam_index11
g14_15_sorted_bam22 = STAR_Module_merge_transcriptome_bam.out.g14_15_sorted_bam22


Adapter_Trimmer_Quality_Module_Umitools_Summary(g17_23_log_file10_g17_24.collect(),g_16_1_g17_24)
g17_24_outputFileTSV00 = Adapter_Trimmer_Quality_Module_Umitools_Summary.out.g17_24_outputFileTSV00


Adapter_Trimmer_Quality_Module_FastQC(g_16_0_g17_28,g_15_1_g17_28)
g17_28_FastQCout00 = Adapter_Trimmer_Quality_Module_FastQC.out.g17_28_FastQCout00


Adapter_Trimmer_Quality_Module_FastQC_after_Adapter_Removal(g_16_0_g17_31,g17_18_reads01_g17_31)
g17_31_FastQCout00 = Adapter_Trimmer_Quality_Module_FastQC_after_Adapter_Removal.out.g17_31_FastQCout00



DESeq2_module_Prepare_DESeq2(g_23_outputFileTSV00_g20_0.first(),g_3_1_g20_0,g_4_2_g20_0,g_21_3_g20_0)
g20_0_rMarkdown03_g20_1 = DESeq2_module_Prepare_DESeq2.out.g20_0_rMarkdown03_g20_1


DESeq2_module_run_DESeq2_markdown(g_4_0_g20_1,g_3_1_g20_1,g_23_outputFileTSV02_g20_1.first(),g20_0_rMarkdown03_g20_1.flatten())
g20_1_rMarkdown00 = DESeq2_module_run_DESeq2_markdown.out.g20_1_rMarkdown00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
