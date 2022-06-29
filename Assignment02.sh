obterSampleName(){
'''
    Obtem o nome de cada amostra do .xml do sra
    e escreve num ficheiro.
'''
rm sample_names.txt   
ids=$(wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi/\?db\=sra\&term\="PRJNA522647"\&usehistory\=y\&retmax=1000 -O - | grep "<Id>" | sed "s/^....//" | rev | sed "s/.....//" | rev)
echo $ids > ./idList.txt
idList=$(cat ./idList.txt)

for id in $idList
do
    wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/\?db\=sra\&id\=$id\&rettype=fasta -O ./entry.xml
    entry=$(cat entry.xml)
    sample_name=$(echo $entry | sed "s/.*sample_name=\"//" | sed "s/\" spots=.*//")
    echo $sample_name >> sample_names.txt
done
}


obterSampleRegion(){
'''
    Obtem a região de cada amostra do .xml do sra
    e escreve num ficheiro.
'''    
rm region.txt
for id in $idList
do
    wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/\?db\=sra\&id\=$id\&rettype=fasta -O ./entry.xml
    entry=$(cat entry.xml)
    sample_region=$(echo $entry | grep -o "<SAMPLE_ATTRIBUTE>.*" | grep -o "<TAG>.*<TAG>" | grep -o "<VALUE>.*</VALUE>"  | sed 's/.*_name//' | grep -o "<VALUE>.*</VALUE>" | sed 's/<\VALUE>\.*//' | sed 's/<.*//' | tr -d ":" | tr -d "," | cut -d ' ' -f 1)
    
    echo $sample_region >> region.txt
done
}


obterAccList(){
'''
    Obtem a referência de cada amostra do .xml do sra
    e escreve num ficheiro.
'''    
rm sraAccList.txt
for id in $idList
do
    wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/\?db\=sra\&id\=$id\&rettype=fasta -O ./entry.xml
    entry=$(cat entry.xml)
    sraAccList=$(echo $entry | sed "s/.*accession=\"SRR/SRR/" | sed "s/\" alias=.*//")
    echo $sraAccList >> sraAccList.txt
done    
}


obterFicheiroSampleNameRegion(){
'''
    Escreve no ficheiro sampleNameRegion.indfile o nome das amostras e a região 
    correspondente.
'''    
    paste sample_names.txt region.txt | column -s $'\t' -t > sampleNameRegion.indfile
}


obterSequenciaDeReferencia(){
'''
    A partir do script python primeiroHomeWork_GuilhermeSilva_Marine_Miguel.py
    obtem a sequência de referência do nosso paiper no NCBI.
'''    
    python3 primeiroHomeWork_GuilhermeSilva_Marine_Miguel.py "Nucleotide" "NC_009628.2" > refSeq.fasta
}


dumpSequencias(){
'''
   A partir do sra-tools e do fiheiro com as referências das amostras
   Obtem todas as sequencias da base de dados do NCBI sra.
''' 
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate sra
    numeroDeCores=$(cat /proc/cpuinfo | grep processor | wc -l)
    paste sraAccList.txt sample_names.txt| while read line sample_name; do
    echo -e "$line"
    echo -e "$sample_name"
    fasterq-dump -p --mem 1G --split-3 --threads $numeroDeCores --skip-technical "$line"
    done
}


headSequencias(){
'''
   Cria um novo ficheiro onde somente cada amostra têm as primeiras 1000000 linhas
   e depois comprime o ficheiro com o nome de cada camelo.
''' 
    
    paste sraAccList.txt sample_names.txt| while read line sample_name; do
    echo -e "$line"
    echo -e "$sample_name"
    head -n 1000000 $line"_1".fastq | gzip > $sample_name"_1".fastq.gz
    head -n 1000000 $line"_2".fastq | gzip > $sample_name"_2".fastq.gz
    done
}


substituirNomeSequencias(){
'''
    Substitui o _R1_ e _R2_ nos nomes dos ficheiros de cada amostra onde está _1 e _2, correpondentemente.
'''
for file in ./*fastq.gz
do
    filename=$(echo "$file" | sed "s/_1/_R1_/" | sed "s/_2/_R2_/")
    echo $filename
    mv $file $filename
done
}


ipyradSequencias(){
'''
    Executa todos os passos do ipyrad de todas as sequências com os nomes corretos, após as alterações, diretamente dos 
    ficheiros comprimidos.
'''  
    source ~/miniconda3/etc/profile.d/conda.sh
    conda deactivate	
    conda activate ipyrad
    numeroDeCores=$(cat /proc/cpuinfo | grep processor | wc -l)
    ipyrad -p config_ipyrad.txt -s 1234567 -c $numeroDeCores -f -r
    echo $numeroDeCores
}


plotingAd(){
'''
    A partir do output do ipyrad é criado um AdmixturePlot com o método "STRUCTURE" a partir do .vcf que se obtem após o ipyrad
    O AdmixturePlot fica guardado numa pasta com os resultados.
'''   
    source ~/miniconda3/etc/profile.d/conda.sh
    conda deactivate
    conda activate structure		
    vcftools --vcf ./camels_project/camels_outfiles/camels.vcf --maf 0.005 --max-missing 0.60 --recode --out camels
    mv camels.recode.vcf camels.vcf
    structure_threader run -i camels.vcf -o ./results_camels -als ~/miniconda3/envs/structure/bin/alstructure_wrapper.R -K 6 -t 12 --ind sampleNameRegion.indfile
}


plotingStructure(){ 
'''
    A partir do output do ipyrad é criado um AdmixturePlot com o método "fastStructure" a partir do .vcf que se obtem após o ipyrad
    O AdmixturePlot fica guardado numa pasta com os resultados.
'''    
    structure_threader run -i ./camels_project/camels_outfiles/camels.str -fs ~/miniconda3/envs/structure/bin/fastStructure -K 6 -o ./camels_results -t 12 --ind sampleNameRegion.indfile
}


pcaModificarFicheiroStr(){
'''
    Ordena o ficheiro sampleNameRegion para  sample_region_for_str 
    Acresenta esse ficheiro ao .str dos camelos e cria o camels_pca.str
    Por fim move para a pasta dos outfiles do ipyrad
'''     
    sort sampleNameRegion.indfile > sample_region_for_str.indfile
    join sample_region_for_str.indfile camels_project/camels_outfiles/camels.str > camels_pca.str
    mv camels_pca.str camels_project/camels_outfiles/
}


executarScriptR(){
'''
    Executa o script pca.r que cria um PCA a partir do camels_pca.str
''' 
    sudo Rscript pca.R
}


obterSampleName
obterSampleRegion
obterFicheiroSampleNameRegion
obterAccList
dumpSequencias
headSequencias
substituirNomeSequencias
obterSequenciaDeReferencia
ipyradSequencias
plotingAd
plotingStructure
pcaModificarFicheiroStr
executarScriptR


