from genbankpy.parser import GenBankFastaWriter


sp_list = ['Emiliania huxleyi']
gbk_data_dir = '/home/robaina/Documents/GenBankpy/gbk_data/'

gbkwriter = GenBankFastaWriter.fromSpecies(species_list=sp_list,
                                           data_dir=gbk_data_dir,
                                           only_representatives=True)



gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['any']},
    output_fasta='results/allPeptidesEmiliania.faa', 
    sequence='protein'
)