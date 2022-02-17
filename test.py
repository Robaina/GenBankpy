from genbankpy.parser import GenBankFastaWriter, GBK


sp_list = ['Emiliania huxleyi']

gbkwriter = GenBankFastaWriter.fromSpecies(species_list=sp_list,
                                           data_dir='/home/robaina/Documents/GenBankpy/')