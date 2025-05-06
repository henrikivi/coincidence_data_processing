import os
import py_io.co_io as coio
import py_io.co_cli as cocli
import peak_funcs.peak_funcs as pfunc

"""
This script is for analysing a single measurement. It needs the -filepath and -param_file arguments
when executing the script.

Example:
python3 analyse_files.py -filepath "200eV/280125aa" -param_file "DCE_params_ds1"
"""

if __name__ == "__main__":
    print("Loading...")
    # get arguments from CLI
    parsey = cocli.parsing()
    parsed_args, _ = parsey.parse_known_args()
    #print(parsed_args)
    #print(os.getcwd())

    file_path = parsed_args.filepath
    full_file = os.path.join(os.getcwd(), file_path)
    co_results = coio.load_eiei_data(full_file)

    results_folder = f"{file_path.split('/')[0]}_results"
    os.makedirs(results_folder, exist_ok=True)

    param_file = parsed_args.param_file
    single_params, pair_params, triple_params, _ = coio.load_params(param_file)
    #print(single_params, pair_params, triple_params)
    #print(single_params[2][3])

    single_results = []
    for row in single_params:
        peak_count = pfunc.count_singles(co_results, row[2], row[3])
        bg_count = pfunc.count_singles(co_results, row[4], row[5])
        single_results.append([row[1], peak_count, bg_count]) # return m/z, peak, background
    
    #print(single_results)
    
    #print(pair_params)
    #print(co_results["pairs"])
    
    #This block of code will return a list of lists, containing the counts and gradient of the
    #pair peaks. The entries will be of form [m/z(1),m/z(2),counts,gradient,gradient_error]
    #If dead time is to be corrected, counts not calculated.
    pairs_results = []
    for row in pair_params:
        pair_coords = [row[4],row[5],row[6],row[7],row[8]]
        #print("\nCoords:",pair_coords)
        selpairs, pair_counts, pair_means = pfunc.select_pairs(co_results,pair_coords,row[9])
        #print("Counts,[mean(t1),mean(t2)]:")
        #print(pair_counts, pair_means)
        if pair_counts == 0: #If no pairs, don't calculate gradient
            pairs_results.append([row[2],row[3],0,None,None])
            continue

        t1, t2 = zip(*selpairs)
        grad_pairs, grad_pairs_err = pfunc.peak_fit_error_x_y(t1, t2, pair_means[0],pair_means[1])

        if row[10] == True: #If dead time needs to be corrected, save time-diff csv plots and counts as None
            hist_vals, bin_edges = pfunc.calc_time_diff_hist(t1, t2)
            pfunc.save_hist_csv(hist_vals, bin_edges,f"{row[2]}+{row[3]}_DTC.csv",results_folder)

            padded = pfunc.dead_time_padding(hist_vals,34,5)
            
            pairs_results.append([row[2],row[3],sum(padded),grad_pairs,grad_pairs_err])
            continue

        pairs_results.append([row[2],row[3],pair_counts,grad_pairs,grad_pairs_err])
    
    
    triples_results = []
    for row in triple_params:
        triple_coords = [row[8],row[9],row[10],row[11],row[12]]
        #print("\nCoords:",triple_coords)
        trip_pairs, trip_counts, trip_means = pfunc.select_triples(co_results,row[6],row[7],triple_coords,row[13])
        #print("Counts,[mean(t1),mean(t2)]:")
        #print(trip_counts, trip_means)

        if trip_counts == 0:
            triples_results.append([row[2],row[3],0,None,None])
            continue

        t1_trip, t2_trip = zip(*trip_pairs)
        grad_trips, grad_trips_err = pfunc.peak_fit_error_x_y(t1_trip, t2_trip, trip_means[0],trip_means[1])
        if row[14] == True: #DTC Check, same as in pairs
            hist_vals, bin_edges = pfunc.calc_time_diff_hist(t1_trip, t2_trip)
            pfunc.save_hist_csv(hist_vals, bin_edges,f"{row[2]}+{row[3]}_DTC.csv",results_folder)

            padded = pfunc.dead_time_padding(hist_vals,34,5)

            triples_results.append(row[3],row[4],row[5],sum(padded),grad_trips,grad_trips_err)
            continue
        triples_results.append([row[3],row[4],row[5],trip_counts,grad_trips,grad_trips_err])
    
    #print(pairs_results)
    #print("\n",triples_results)

    for row in single_results:
        print(row[0], row[1])

    for row in pairs_results:
        print(row[0],row[1],row[2])

    print()

    for row in triples_results:
        print(row[0],row[1],row[2],row[3])


    coio.save_results(single_results,pairs_results,triples_results,f"{results_folder}.csv",results_folder)

