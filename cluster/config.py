import configparser

config = configparser.ConfigParser()

config.add_section('postgresql')
config.set('postgresql', 'number_of_clusters', '3')
config.set('postgresql', 'number_of_vector_hash_tables', '3')
config.set('postgresql', 'number_of_vector_hash_functions', '4')
config.set('postgresql', 'max_number_M_hypercube', '10')
config.set('postgresql', 'number_of_hypercube_dimensions', '3')
config.set('postgresql', 'number_of_probes', '2')

with open("configuration.ini", 'w') as configfile:
    config.write(configfile)