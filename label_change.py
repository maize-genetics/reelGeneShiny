import os


def change_label(fasta_file, dictionary):
    # Read the contents of the FASTA file
    with open(fasta_file, 'r') as file:
        lines = file.readlines()

    # Change the label names based on the dictionary
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith('>'):
            label = line[1:]  # Extract the label name
            # Use the dictionary to get the new label name
            new_label = dictionary.get(label, label)
            # Update the label name in the line
            lines[i] = '>' + new_label + '\n'

    # Create a new file with the modified label names
    file_name = os.path.basename(fasta_file)
    new_file = os.path.join(new_direct, file_name)
    with open(new_file, 'w') as file:
        file.writelines(lines)

    print(f"Created a new file: {new_file}")


# Example dictionary mapping old label names to new label names
label_dictionary = {
    'B73_REF_mRNA': 'B73_REF_mRNA',
    'B73_REF_FULL': 'B73_REF_FULL',
    'Ab-Traiperm_572-DRAFT-PanAnd-1': 'Andropogon burmanicus (Long)',
    'Ac-Pasquet1232-DRAFT-PanAnd-1': 'Andropogon chinensis (Long)',
    'Ag-CAM1351-DRAFT-PanAnd-1': 'Andropogon gerardi (Long)',
    'AN20T007': 'Heteropogon contortus (Short)',
    'AN20T012': 'Themeda avenacea (Short)',
    'AN20T014': 'Bothriochloa decipiens (Short)',
    'AN20T018': 'Mnesithea formosa (Short)',
    'AN20T024': 'Saccharum hildebrandtii (Short)',
    'AN20T168': 'Andropogon amethystinus (Short)',
    'AN20T174': 'Arthraxon antsirabensis (Short)',
    'AN20T184': 'Diheteropogon filifolius (Short)',
    'AN20T190': 'Hemarthria altissima (Short)',
    'AN20T197': 'Microstegium nudum (Short)',
    'AN21TNTL0089': 'Pseudopogonatherum contortum (Short)',
    'AN21TNTL0090': 'Phacelurus franksiae (Short)',
    'AN21TNTL0091': 'Sorghastrum fuscescens (Short)',
    'AN21TNTL0093': 'Eulalia monostachya (Short)',
    'AN21TNTL0095': 'Rottboellia afraurita (Short)',
    'AN21TNTL0097': 'Eremochloa lanceolata (Short)',
    'AN21TNTL0098': 'Eremochloa attenuata (Short)',
    'AN21TNTL0102': 'Cymbopogon bhutanicus (Short)',
    'AN21TNTL0106': 'Arthraxon lancifolius (Short)',
    'AN21TNTL0108': 'Arthraxon meeboldii (Short)',
    'AN21TNTL0109': 'Arthraxon microphyllus (Short)',
    'AN21TNTL0115': 'Clausospicula extensa (Short)',
    'AN21TNTL0121': 'Lasiurus scindicus (Short)',
    'AN21TNTL0124': 'Elionurus elegans (Short)',
    'AN21TNTL0126': 'Germainia capitata (Short)',
    'AN21TNTL0131': 'Tripsacum latifolium (Short)',
    'AN21TNTL0135': 'Thaumastochloa striata (Short)',
    'AN21TNTL0139': 'Dimeria aristata (Short)',
    'AN21TNTL0148': 'Exotheca abyssinica (Short)',
    'AN21TNTL0153': 'Microstegium geniculatum (Short)',
    'AN21TNTL0154': 'Microstegium fasciculatum (Short)',
    'AN21TNTL0157': 'Miscanthus nepalensis (Short)',
    'AN21TNTL0163': 'Ischaemum byrone (Short)',
    'AN21TNTL0167': 'Iseilema vaginiflorum (Short)',
    'AN21TNTL0169': 'Diectomis fastigiata (Short)',
    'AN21TNTL0176': 'Eriochrysis cayennensis (Short)',
    'AN21TNTL0178': 'Sorghum arundinaceum (Short)',
    'Av-Kellogg1287_8-REFERENCE-PanAnd-1': 'Andropogon virginicus (Long)',
    'Bl-K1279B-DRAFT-PanAnd-1': 'Bothriochloa laguroides (Long)',
    'Cc-PI314907-DRAFT-PanAnd-1': 'Cymbopogon citratus (Long)',
    'Cr-AUB069-DRAFT-PanAnd-1': 'Cymbopogon refractus (Long)',
    'Et-Layton_Zhong168-DRAFT-PanAnd-1': 'Elionurus tripsacoides (Long)',
    'Hc-AUB53_1-DRAFT-PanAnd-1': 'Heteropogon contortus (Long)',
    'Hp-KelloggPI404118-DRAFT-PanAnd-1': 'Hemarthria compressa (Long)',
    'Ir-Pasquet1136-DRAFT-PanAnd-1': 'Ischaemum rugosum (Long)',
    'MCRTL003': 'Euclasta condylotaicha (Short)',
    'MCRTL004': 'Dichanthium foveolatum (Short)',
    'MCRTL005': 'Polytoca digitata (Short)',
    'MCRTL008': 'Kerriochloa siamensis (Short)',
    'MCRTL009': 'Tripsacum australe (Short)',
    'MCRTL011': 'Urelytrum agropyroides (Short)',
    'MCRTL013': 'Hyperthelia dissoluta (Short)',
    'MCRTL015': 'Eulaliopsis binata (Short)',
    'MCRTL018': 'Parahyparrhenia siamensis (Short)',
    'MCRTL019': 'Eremochloa ciliaris (Short)',
    'MCRTL020': 'Zea diploperennis (MCRTL020) (Short)',
    'MCRTL022': 'Eremochloa eriopoda (Short)',
    'MCRTL023': 'Trachypogon chevaleri (Short)',
    'MCRTL026': 'Hyparrhenia bracteata (Short)',
    'MCRTL027': 'Schizachyrium brevifolium (Short)',
    'MCRTL031': 'Anadelphia scyphofera (Short)',
    'MCRTL033': 'Pseudosorghum fasiculare (Short)',
    'MCRTL038': 'Schizachyrium reedii (Short)',
    'MCRTL053': 'Hackelochloa granularis (Short)',
    'MCRTL059': 'Elymandra archaelymandra (Short)',
    'MCRTL063': 'Schizachyrium delicatum (Short)',
    'Pi-Clark-DRAFT-PanAnd-1': 'Pogonatherum paniceum  (Long)',
    'Rr-Malcomber3106-DRAFT-PanAnd-1': 'Rhytachne rottboellioides (Long)',
    'Rt-Layton_Zhong169-DRAFT-PanAnd-1': 'Rottboellia tuberculosa (Long)',
    'Sb-JGI-v3': 'Sorghum bicolor (Long)',
    'Sm-PI203595-DRAFT-PanAnd-1': 'Schizachyrium microstachyum (Long)',
    'Sn-CAM1369-DRAFT-PanAnd-1': 'Sorghastrum nutans (Long)',
    'Ss-CAM1384-DRAFT-PanAnd-1': 'Schizachyrium scoparium (Long)',
    'Td-FL_9056069_6-DRAFT-PanAnd-1': 'Tripsacum dactyloides (FL_9056069_6) (Long)',
    'Td-KS_B6_1-DRAFT-PanAnd-1': 'Tripsacum dactyloides (KS_B6_1) (Long)',
    'Te-Pasquet1246-DRAFT-PanAnd-1': 'Thelepogon elegans (Long)',
    'Tt-AUB21_1-DRAFT-PanAnd-1': 'Themeda triandra (Long)',
    'Ud-Pasquet1171-DRAFT-PanAnd-1': 'Urelytrum digitatum (Long)',
    'Vc-Pasquet1098-DRAFT-PanAnd-1': 'Vossia cuspidata (Long)',
    'Zd-Gigi-REFERENCE-PanAnd-1': 'Zea diploperennis (Gigi) (Long)',
    'Zd-Momo-REFERENCE-PanAnd-1': 'Zea diploperennis (Momo) (Long)',
    'Zh-RIMHU001-REFERENCE-PanAnd-1': 'Zea mays subsp. huehuetenangensis (Long)',
    'Zl-RIL003-REFERENCE-PanAnd-1': 'Zea luxurians (Long)',
    'Zn-PI615697-REFERENCE-PanAnd-1': 'Zea nicaraguensis (Long)',
    'Zv-TIL01-REFERENCE-PanAnd-1': 'Zea mays subsp. parviglumis (TIL01) (Long)',
    'Zv-TIL11-REFERENCE-PanAnd-1': 'Zea mays subsp. parviglumis (TIL11) (Long)',
    'Zx-TIL18-REFERENCE-PanAnd-1': 'Zea mays subsp. mexicana (TIL18) (Long)',
    'Zx-TIL25-REFERENCE-PanAnd-1': 'Zea mays subsp. mexicana (TIL25) (Long)'
}

# Directory containing the FASTA files
fasta_directory = 'Buckler/transcripts/'

# Directory for new files
new_direct = 'Buckler/transcripts/new/'

# Create the new directory if it doesn't exist
os.makedirs(new_direct, exist_ok=True)

# Process each FASTA file in the directory
for file_name in os.listdir(fasta_directory):
    if file_name.endswith('.fa'):
        file_path = os.path.join(fasta_directory, file_name)
        change_label(file_path, label_dictionary)
