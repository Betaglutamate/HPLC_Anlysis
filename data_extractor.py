import os
import pandas as pd

# assign directory
directory = 'data/20230322_Trehalose_Acetate'
 
# iterate over files in
all_dfs = []

for root, dirs, files in os.walk(directory):
    for filename in files:
        current_file =(os.path.join(root, filename))
        current_df = pd.read_csv(current_file)
        beginning_data = current_df[current_df['[Header]'].str.contains('Peak Table')].index[0]
        end_data = current_df[current_df['[Header]'].str.contains('Compound Results')].index[0]
        new_df = current_df.iloc[beginning_data +2: end_data]
        if not new_df.empty:
            split = new_df['[Header]'].str.split("\t", expand = True)
            headers = split.iloc[0].values
            split.columns = headers
            split = split[1:]
            split['filename'] = filename
            all_dfs.append(split)

final_df = pd.concat(all_dfs)

final_df["compound"] = "unknown"




#             for line_number, line in enumerate(file_extracting):
#                 if "Compound Results" in line:
#                     line_start = line_number
#                 if "Group Results" in line:    
#                     line_end = line_number-1
#                 for line_number, line in enumerate(file_extracting):
#                     if line_number >+ line_start and line_number < line_end:
#                         lines.append(line)
                
                


#                     test_dict['']

                
#                     print(line_number)
#                     print(filename) 