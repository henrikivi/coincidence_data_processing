import os
import glob
import py_io.co_cli as cocli
import py_io.co_io as coio
import peak_funcs.peak_funcs as pfunc

"""
This script is for analysing an entire dataset. It needs the -data_folder and -param_file arguments to run.

Example:
python3 analyse_dataset.py -data_folder "dataset_1" -param_file "DCE_params_ds1"

The data folder needs to be set up such that the .mas-, .prs.-, and .trp-files for different energies are
organised to separate subfolders.
"""

if __name__ == "__main__":
    print("Loading...")
    # Get arguments from CLI
    parsey = cocli.parsing()
    parsed_args, _ = parsey.parse_known_args()

    # Extract data folder and parameter file
    data_folder = parsed_args.data_folder  
    param_file = parsed_args.param_file

    # Load parameters
    single_params, pair_params, triple_params, shift_params = coio.load_params(param_file)

    # Create a results folder
    results_folder = f"{data_folder}_results"
    os.makedirs(results_folder, exist_ok=True)

    # Iterate through each measurement folder in the data folder
    for measurement_folder in os.listdir(data_folder):
        measurement_path = os.path.join(data_folder, measurement_folder)

        # Check if it's a directory
        if not os.path.isdir(measurement_path):
            print(f"Warning: {measurement_path} is not a directory. Skipping this path.")
            continue
        
        print(f"Processing folder: {measurement_folder}")

        for row in shift_params:
            if row[0] == measurement_folder:
                shif = pfunc.calc_shift_amount(
                    row[1], #t1
                    row[2], #t2
                    shift_params[0][1], #reft1
                    shift_params[0][2] #reft2
                    )
                break
        print(f"Shift for dataset {measurement_folder}: {shif}\n")

        #Initialise shifted parameter arrays (copy param_file data)
        single_params_shif = [row[:] for row in single_params]
        pair_params_shif = [row[:] for row in pair_params]
        triple_params_shif = [row[:] for row in triple_params]

        #Apply shift
        idx_singles_update = [2,3,4,5]
        for row in single_params_shif:
            for idx in idx_singles_update:
                row[idx] += shif

        idx_pairs_update = [4,5,6,7]
        for row in pair_params_shif:
            for idx in idx_pairs_update:
                row[idx] += shif

        idx_trips_update = [6,7,8,9,10,11]
        for row in triple_params_shif:
            for idx in idx_trips_update:
                row[idx] += shif

        # print(single_params[0])
        # print(single_params_shif[0])
        # print()
        # print(pair_params[0])
        # print(pair_params_shif[0])
        # print()
        # print(triple_params[0])
        # print(triple_params_shif[0])

        measurement_results_folder = os.path.join(results_folder, measurement_folder)
        os.makedirs(measurement_results_folder, exist_ok=True)

        # Search for the .mas file in the measurement folder
        mas_files = glob.glob(os.path.join(measurement_path, "*.mas"))
        if not mas_files:
            print(f"Warning: No .mas files found in {measurement_path}. Skipping this folder.")
            continue

        # Extract the name of the .mas file for processing and load experiment data
        full_file = os.path.splitext(mas_files[0])[0]
        co_results = coio.load_eiei_data(full_file)

        single_results = []
        for row in single_params_shif:
            peak_count = pfunc.count_singles(co_results, row[2], row[3])
            bg_count = pfunc.count_singles(co_results, row[4], row[5])
            single_results.append([row[1], peak_count, bg_count])  # return m/z, peak, background

        pairs_results = []
        for row in pair_params_shif:
            pair_coords = [row[4], row[5], row[6], row[7], row[8]]
            selpairs, pair_counts, pair_means = pfunc.select_pairs(co_results, pair_coords, row[9])
            if pair_counts == 0:  # If no pairs, don't calculate gradient
                pairs_results.append([row[2], row[3], 0, None, None])
                continue

            t1, t2 = zip(*selpairs)
            grad_pairs, grad_pairs_err = pfunc.peak_fit_error_x_y(t1, t2, pair_means[0], pair_means[1])

            if row[10]:  # Correct dead-time, save raw time-diff csv:s
                hist_vals, bin_edges = pfunc.calc_time_diff_hist(t1, t2)
                pfunc.save_hist_csv(hist_vals, bin_edges,f"{row[2]}+{row[3]}_DTC.csv",measurement_results_folder)

                padded = pfunc.dead_time_padding(hist_vals, 34, 5)

                pairs_results.append([row[2], row[3], sum(padded), grad_pairs, grad_pairs_err])
                continue

            pairs_results.append([row[2], row[3], pair_counts, grad_pairs, grad_pairs_err])

        triples_results = []
        for row in triple_params_shif:
            triple_coords = [row[8], row[9], row[10], row[11], row[12]]
            trip_pairs, trip_counts, trip_means = pfunc.select_triples(co_results, row[6], row[7], triple_coords, row[13])
            #print("\nCounts,[mean(t1),mean(t2)]:")
            #print(trip_counts, trip_means)
            if trip_counts == 0:
                triples_results.append([row[3], row[4],row[5], 0, None, None])
                continue

            t1_trip, t2_trip = zip(*trip_pairs)
            grad_trips, grad_trips_err = pfunc.peak_fit_error_x_y(t1_trip, t2_trip, trip_means[0], trip_means[1])

            if row[14]:  # Correct dead-time, save uncorrected time-diff csv:s
                hist_vals, bin_edges = pfunc.calc_time_diff_hist(t1_trip, t2_trip)
                pfunc.save_hist_csv(hist_vals, bin_edges,f"{row[2]}+{row[3]}_DTC.csv",measurement_results_folder)

                padded = pfunc.dead_time_padding(hist_vals, 34, 5)

                triples_results.append([row[3], row[4], row[5], sum(padded), grad_trips, grad_trips_err])
                continue

            triples_results.append([row[3], row[4], row[5], trip_counts, grad_trips, grad_trips_err])

        # Save results for this measurement
        coio.save_results(single_results,pairs_results,triples_results,f"{measurement_results_folder}.csv")
