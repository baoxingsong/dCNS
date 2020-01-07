//
// Created by Baoxing Song on 2019-04-06.
//

#include "gffToCategory.h"




// for those duplication gff records and neighbour records with position over, try to merge them and generate a list of unique and non-overlapping records
void removeDuplications( std::map<std::string, std::vector<element>> & elements0, std::map<std::string, std::vector<element> > & elements ){
    for( std::map<std::string, std::vector<element>>::iterator it = elements0.begin(); it != elements0.end(); ++it ){
        elements[it->first] = std::vector<element>();
        std::vector<element> elementv = elements0[it->first];
        if (0 == elementv.size()) {
            continue;
        }

        std::sort(elementv.begin(), elementv.end(), [](element a, element b) {
            return a.start < b.start;
        });

        element e=elementv[0];
        for ( int i=1; i< elementv.size(); ++i ){
            if( (e.end) >= (elementv[i].end) ){

            } else if ( elementv[i].start <= e.end ) {
                e.end = elementv[i].end >  e.end ? elementv[i].end :  e.end;
            }else{
                elements[it->first].push_back(e);
                e=elementv[i];
            }
        }
    }
}

// update the weigth vector (categories)
// weigths is a vector of element and the key is the chromosome
// weigth is the weight of this category
void updateWeight( std::map<std::string, std::vector<element> > & elements, std::map<std::string, int16_t *> & weigths, const int16_t & weigth){
    int i, j;
    for( std::map<std::string, std::vector<element>>::iterator it = elements.begin(); it != elements.end(); ++it ){
        for( i=0; i<it->second.size(); ++i ){
            for( j=it->second[i].start; j<=it->second[i].end; ++j ) {
                weigths[it->first][j]=weigth;
            }
        }
    }
}

//read gff file and write the weight vector into a file
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, const std::string & outputFile){
    std::map<std::string, int16_t *> categories;
    readGffFileWithEveryThing (filePath,  chrSize, categories);

    std::ofstream ofile;
    ofile.open(outputFile);
    int elementIndex;
    for( std::map<std::string, int32_t>::iterator  it = chrSize.begin(); it != chrSize.end(); ++it ){
        ofile << ">" << it->first << "\n";
        elementIndex = 0;
        for( int16_t i=0; i<it->second; ++i ){
            ++elementIndex;
            if( 0 == elementIndex%1000 ){
                ofile << i << "\n";
            }else{
                ofile << i << "\t";
            }
        }
    }
    ofile.close();
    for( std::map<std::string, int16_t *>::iterator it=categories.begin(); it != categories.end(); ++it ){
        delete[] it->second;
    }
}



// read gff file and generate a weight vector
void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, std::map<std::string, int16_t *> & categories){

    std::map<std::string, std::vector <element>> introns0; // gene, transcript protein, intron
    std::map<std::string, std::vector <element>> exons0;
    std::map<std::string, std::vector <element>> cds0;
    std::map<std::string, std::vector <element>> transposons0;

    int32_t i;

    for( std::map<std::string, int32_t>::iterator it = chrSize.begin(); it!=chrSize.end(); ++it ){
        categories[it->first] =  new int16_t[it->second]; //(, 10); //10 is the weight of intergenetic
        for( i=0; i < it->second; ++i ){
            categories[it->first][i] = 10; //intergenetic region
        }
        introns0[it->first]=std::vector <element>();
        exons0[it->first]=std::vector <element>();
        cds0[it->first]=std::vector <element>();
        transposons0[it->first]=std::vector <element>();
    }

    std::set<std::string> intronsNickNames;
    intronsNickNames.insert("transcript");
    intronsNickNames.insert("pseudogenic_transcript");
    intronsNickNames.insert("mRNA");
    intronsNickNames.insert("miRNA");
    intronsNickNames.insert("tRNA");
    intronsNickNames.insert("ncRNA");
    intronsNickNames.insert("mRNA_TE_gene");
    intronsNickNames.insert("rRNA");
    intronsNickNames.insert("snoRNA");
    intronsNickNames.insert("snRNA");
    intronsNickNames.insert("lincRNA");
    intronsNickNames.insert("gene");
    intronsNickNames.insert("pseudogene");
    intronsNickNames.insert("transposable_element_gene");
    intronsNickNames.insert("lincRNA_gene");
    intronsNickNames.insert("tRNA_gene");
    intronsNickNames.insert("protein");
    intronsNickNames.insert("miRNA_gene");
    intronsNickNames.insert("ncRNA_gene");
    intronsNickNames.insert("lnc_RNA");
    intronsNickNames.insert("pre_miRNA");
    intronsNickNames.insert("SRP_RNA");

    std::set<std::string> exonsNickNames;
    exonsNickNames.insert("exon");
    exonsNickNames.insert("five_prime_utr");
    exonsNickNames.insert("five_prime_UTR");
    exonsNickNames.insert("three_prime_UTR");
    exonsNickNames.insert("three_prime_utr");

    std::set<std::string> cdssNickNames;
    cdssNickNames.insert("CDS");
    cdssNickNames.insert("start_codon");
    cdssNickNames.insert("stop_codon");

    std::set<std::string> transposableNickNames;
    transposableNickNames.insert("transposable_element");
    transposableNickNames.insert("transposon_fragment");
    transposableNickNames.insert("repeat_region");
    transposableNickNames.insert("LINE_element");
    transposableNickNames.insert("LTR_retrotransposon");
    transposableNickNames.insert("SINE_element");
    transposableNickNames.insert("helitron");
    transposableNickNames.insert("solo_LTR");
    transposableNickNames.insert("terminal_inverted_repeat_element");

    std::set<std::string> ignoreTypes;
    ignoreTypes.insert("chromosome");
    ignoreTypes.insert("contig");


    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }

    char splim='\t';
    std::string line;
    std::string chromosomeName;
    while (std::getline(infile, line)){
        if(line[0]=='#' || line.size()<9){
        }else {
            std::vector<std::string> elemetns;
            split(line, splim, elemetns);
            if(elemetns.size()==9){
                int start = stoi(elemetns[3]);
                int end = stoi(elemetns[4]);
                if (start > end) {
                    int temp = start;
                    start = end;
                    end = temp;
                }
                --start; // change the coordinates to 0 based
                --end;
                chromosomeName = elemetns[0];
                if (intronsNickNames.find(elemetns[2]) != intronsNickNames.end()) {
                    introns0[chromosomeName].push_back(element{start, end});
                } else if (exonsNickNames.find(elemetns[2]) != exonsNickNames.end()) {
                    exons0[chromosomeName].push_back(element{start, end});
                } else if ( cdssNickNames.find(elemetns[2]) != cdssNickNames.end() ) {
                    cds0[chromosomeName].push_back(element{start, end});
                } else if ( transposableNickNames.find(elemetns[2]) != transposableNickNames.end() ) { // ignore those elements
                    transposons0[chromosomeName].push_back(element{start, end});
                } else if ( ignoreTypes.find(elemetns[2]) != ignoreTypes.end() ) { // ignore those elements

                } else {
                    std::cout << "we could not analysis the line in the gff/gtf file: " << line << std::endl;
                }
            }
        }
    }
    time_t my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "gff file reading done" << std::endl;

    std::map<std::string, std::vector <element>> introns; // gene, transcript protein, intron
    std::map<std::string, std::vector <element>> exons;
    std::map<std::string, std::vector <element>> cds;
    std::map<std::string, std::vector <element>> transposons;

    removeDuplications(introns0, introns);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "gff file intron remove duplication done" << std::endl;
    removeDuplications(exons0, exons);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "gff file exons remove duplication done" << std::endl;
    removeDuplications(cds0, cds);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "gff file cds remove duplication done" << std::endl;
    removeDuplications(transposons0, transposons);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "gff file transposons remove duplication done" << std::endl;

    updateWeight( introns0, categories, 70);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update weigth introns done" << std::endl;
    updateWeight( exons, categories, 85);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update weigth exon done" << std::endl;
    updateWeight( cds0, categories, 100);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update weigth cds done" << std::endl;
    updateWeight( transposons0, categories, 1);
    my_time = time(NULL); printf("%s", ctime(&my_time));
    std::cout << "update weigth transposons done" << std::endl;
}

void readGffFileWithEveryThing (const std::string& filePath, std::map<std::string, int32_t> & chrSize, std::map<std::string, int16_t *> & categories, std::map<std::string, int16_t *> & categories_rev){
    readGffFileWithEveryThing (filePath, chrSize, categories);
    int i;
    for( std::map<std::string, int32_t>::iterator it = chrSize.begin(); it!=chrSize.end(); ++it ){
        categories_rev[it->first] =  new int16_t[it->second];
        for( i=0; i < it->second; ++i ){
            categories_rev[it->first][it->second-1-i] = categories[it->first][i];
        }
    }
}

