import json
def metadata_to_JSON(metadata_path, in_min, out_max, graphpath, count_kmers, outpath):
    # Read the metadata file
    with open(metadata_path, 'r') as file:
        lines = file.readlines()
    # Extract column names from the first line
    columns = lines[0].strip().split('\t')[1:]
    # Create a dictionary to store the row names for each category
    features = {column: {} for column in columns}
    # Process each line in the file
    for line in lines[1:]:
        row_data = line.strip().split('\t')
        row_name = row_data[0]
        # Process each column/category
        for i, feature in enumerate(columns):
            category = row_data[i + 1]
            # Add the row name to the respective category
            if category!='NA' and category!='NaN':
                if category not in features[feature]:
                    features[feature][category] = []
                features[feature][category].append(row_name)
    # Save the category data as separate JSON files
    for feature in features:
        for category in features[feature]:
            file_name = f'{feature}_{category}.json'
            with open(outpath+"/"+file_name, 'w') as json_file:
                data = {
                    "groups": [
                        {
                            "experiments": [
                                {
                                    "name": feature+"_"+category,
                                    "in_min_fraction": in_min,
                                    "out_max_fraction": out_max,
                                    "count_kmers": count_kmers,
                                    "in": [graphpath + "/" + f + ".fasta.gz" for f in features[feature][category]],
                                    "out": [graphpath + "/" + j +".fasta.gz" for j in [f for key, f in features[feature].items() if key!=category][0]]
                                }
                            ]
                        }
                    ]
                }
                json.dump(data, json_file, indent=4)