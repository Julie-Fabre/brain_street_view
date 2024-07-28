import requests
import json

def find_connectivity_experiments(regions, mouse_line='', primary_injection=True):
    experiment_ids = []

    for region in regions:
        # Build the query URL
        base_url = 'http://api.brain-map.org/api/v2/data/query.json?criteria=service::mouse_connectivity_injection_structure'
        
        # Add region of interest
        base_url += f'[injection_structures$eq{region}]'
        
        # Mouse line
        if mouse_line:
            base_url += f'[transgenic_lines$eq{mouse_line}]'
        
        # Primary structure or not
        primary = 'true' if primary_injection else 'false'
        full_url = f'{base_url}[primary_structure_only$eq{primary}]'

        # Get data from Allen
        try:
            response = requests.get(full_url)
            response.raise_for_status()  # Raises an HTTPError for bad requests
            result = response.json()
        except requests.RequestException as e:
            print(f"Request failed: {e}")
            continue
        except json.JSONDecodeError:
            print("Failed to decode JSON response")
            continue

        if not result.get('success'):
            print(f'Query failed!\n{result.get("msg", "No message")}\nAt URL: {full_url}\n')
            continue

        # Get the experiment IDs
        for item in result.get('msg', []):
            experiment_ids.append(item['id'])

        print(f'Found {len(result.get("msg", []))} experiments in {region}')

    return experiment_ids


import os
import numpy as np
import pandas as pd
import bsv  # Assuming you have a Python module for bsv functions

def fetch_connectivity_data(experiment_ids, save_location, file_name, 
                            normalization_method, subtract_other_hemisphere, 
                            grouping_method=None, allen_atlas_path=None):
    
    projection_grid_size = [132, 80, 114]

    # Filepaths
    file_path_imgs = os.path.join(save_location, f"{file_name}_{normalization_method}_sub{int(subtract_other_hemisphere)}.npy")
    file_path_injection_summary = os.path.join(save_location, f"{file_name}_injectionSummary.csv")

    if not os.path.exists(file_path_imgs) or not file_name:
        combined_injection_info = pd.DataFrame()

        for exp_id in experiment_ids:
            save_dir = os.path.join(save_location, str(exp_id))
            os.makedirs(save_dir, exist_ok=True)

            # Summary (structure.ionizes)
            summary_file_path = os.path.join(save_dir, 'injectionSummary_all.npy')
            if not os.path.exists(summary_file_path):
                status = bsv.fetch_connectivity_summary(exp_id, save_dir)
                if not status:
                    continue

            injection_info = np.load(summary_file_path, allow_pickle=True).item()
            
            # Combine injection info
            if combined_injection_info.empty:
                combined_injection_info = pd.DataFrame(injection_info)
            else:
                combined_injection_info = pd.concat([combined_injection_info, pd.DataFrame(injection_info)], ignore_index=True)

        # Grouping method
        if grouping_method == 'brainRegion' and allen_atlas_path:
            st = pd.read_csv(os.path.join(allen_atlas_path, 'structure_tree_safe_2017.csv'))
            output_acronyms = [st.loc[st['id'] == id, 'acronym'].values[0] for id in combined_injection_info['structure_id']]
            
            sorted_values = st.loc[st['id'].isin(combined_injection_info['structure_id']), 'acronym']
            sort_indices = sorted_values.index
            groups, group_id = np.unique(sorted_values, return_inverse=True)
            group_id_original_order = np.zeros_like(group_id)
            group_id_original_order[sort_indices] = group_id
        else:
            print("Grouping method not recognized - skipping grouping.")
            groups = np.ones(len(experiment_ids))
            group_id_original_order = np.ones(len(experiment_ids))

        number_of_groups = len(np.unique(groups))
        combined_projection = np.zeros(projection_grid_size + [number_of_groups])

        # Get raw images
        for i, exp_id in enumerate(experiment_ids):
            this_group = group_id_original_order[i]

            # Raw data
            raw_file_path = os.path.join(save_location, str(exp_id), 'density.raw')
            if not os.path.exists(raw_file_path):
                status = bsv.fetch_connectivity_images(exp_id, os.path.join(save_location, str(exp_id)))
                if not status:
                    continue

            with open(raw_file_path, 'rb') as f:
                experiment_projection = np.fromfile(f, dtype=np.float32).reshape(projection_grid_size)

            if subtract_other_hemisphere:
                experiment_projection_tmp = np.zeros((132, 80, 114))
                if injection_info['max_voxel_z'] <= 114 / 2:  # left
                    for i_ml in range(57):
                        experiment_projection_tmp[:, :, i_ml] = experiment_projection[:, :, i_ml] - experiment_projection[:, :, 113-i_ml]
                elif injection_info['max_voxel_z'] >= 114 / 2:  # right
                    for i_ml in range(57):
                        experiment_projection_tmp[:, :, 113-i_ml] = experiment_projection[:, :, 113-i_ml] - experiment_projection[:, :, i_ml]

                experiment_projection = experiment_projection_tmp

            combined_projection[:,:,:,this_group] += experiment_projection

        # Normalize images to number of projections
        for i_group in range(number_of_groups):
            these_groups = np.sum(group_id_original_order == i_group)
            combined_projection[:,:,:,i_group] /= these_groups

        # Save files
        if file_name:
            np.save(file_path_imgs, combined_projection)
            combined_injection_info.to_csv(file_path_injection_summary, index=False)
    else:
        combined_projection = np.load(file_path_imgs)
        combined_injection_info = pd.read_csv(file_path_injection_summary)

    return combined_projection, combined_injection_info