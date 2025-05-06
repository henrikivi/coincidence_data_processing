"""Module of useful functions to look at pipico peaks


"""
import os
import math
import numpy as np
from scipy.ndimage import center_of_mass
import csv


def count_singles(results:dict, t1:int, t2:int):
    """
    Function to extract counts between t1 & t2 in singles array

    Parameters
    ----------
    results: loaded from a eiei data set, it must contain a dictionary of singles
    t1, t2: lower and upper time values for peak boundaries

    """
    counts = sum(results["singles"][t1:t2+1])
    return counts


def select_pairs(results, box_coord, box_shape=False):
    """Function to select pairs that fall within a box region

    Parameters
    ----------
    results : dict
        Results loaded from a eiei data set, it must contain a dictionary of pairs
    box_coord : list
        list of the box coordinates to make selection in correct units, is in the form
        box_coord[0]: x1
        box_coord[1]: y1
        box_coord[2]: x2
        box_coord[3]: y2
        box_coord[4]: gradient

    box_shape : bool, optional
        If true is a flat topped box, by default False

    Returns
    -------
    lists
        Selected pairs, total counts and the averages inside the box
    """
    spairs = []
    # Get other values for box
    Trap = [0, 0]
    # Work out line constants
    Trap[0] = box_coord[1] - (box_coord[4] * box_coord[0])  # ca
    Trap[1] = box_coord[3] - (box_coord[4] * box_coord[2])  # cb
    
    # print(Trap)
    counts = 0
    means = [0, 0]
    for pair in results["pairs"]:
        t1, t2 = pair  # Extract x, y

        if box_shape:
            # Flat-topped trapezoid
            if box_coord[1] < t2 < box_coord[3]:  # y1 < y < y2
                x_left = (t2 - Trap[0]) / box_coord[4]  # Left boundary x given y
                x_right = (t2 - Trap[1]) / box_coord[4]  # Right boundary x given y
                if x_left < t1 < x_right:  # Check if within the trapezoidal region
                    spairs.append(pair)
                    counts += 1
                    means[0] += t1
                    means[1] += t2
        else:
            # Flat-sided trapezoid
            if box_coord[0] < t1 < box_coord[2]:  # x1 < x < x2
                ct = t2 - (box_coord[4] * t1)  # Calculate ct
                if Trap[0] < ct < Trap[1]:  # Check if within trapezoidal region
                    spairs.append(pair)
                    counts += 1
                    means[0] += t1
                    means[1] += t2
        
    # Compute means if at least one pair was selected
    if counts > 0:
        means[0] /= counts
        means[1] /= counts

    return spairs, counts, means


def printdecstats(descrip, axisnumber, axislabel, neatflag):
    """
    Function to print out some descriptive stats for the least squares fit
    Inputs

    """
    if neatflag:
        print(f"{axislabel} Mean: {descrip.mean[axisnumber]}")
        print(f"{axislabel} Var: {descrip.variance[axisnumber]}")
        print(f"{axislabel} Std: {math.sqrt(descrip.variance[axisnumber])}")
        print(f"{axislabel} Kurtosis: {descrip.kurtosis[axisnumber]}")
        print(f"{axislabel} Skewness: {descrip.skewness[axisnumber]}")
    else:
        print("Mean \t Var \t Std \t Kurtosis \t Skewness \t")
        print(
            f"{descrip.mean[axisnumber]} \t {descrip.variance[axisnumber]} \t {math.sqrt(descrip.variance[axisnumber])} \t {descrip.kurtosis[axisnumber]} \t {descrip.skewness[axisnumber]} \t"
        )


def peak_fit_error_x_y(x, y, mean_x, mean_y):
    """Fit the gradient of the data assuming errors in both x and y

    Parameters
    ----------
    x : float
        the x values
    y : float
        the y values
    mean_x : float
        mean of the x values
    mean_y : float
        mean of the y values

    Returns
    -------
    float
        the gradient and error of the fit
    """
    sumu = 0
    sumv = 0
    sumuv = 0
    if len(x) > 1:
        for xv, yv in zip(x, y):
            ui = xv - mean_x
            vi = yv - mean_y

            sumu = sumu + (ui * ui)
            sumv = sumv + (vi * vi)

            sumuv = sumuv + ui * vi
        
        if sumuv == 0 or sumu * sumv == 0:
            return 0, 0

        grad = (sumv - sumu + math.sqrt((sumv - sumu) ** 2 + 4 * (sumuv) ** 2)) / (
            2 * sumuv
        )
        rb = sumuv / (math.sqrt(sumu * sumv))
        grad_err = (grad / rb) * (math.sqrt((1 - rb**2) / len(x)))

        return grad, grad_err
    else:
        return 0, 0


def line_of_fit(x, y, mean_x, mean_y, variance_x, variance_y, covariance):
    """
    Calculate the line of best fit through a peak
    """

    # gradient
    grad = covariance / variance_x
    # intercept
    inte = mean_y - (grad * mean_x)

    # errors
    # b = (sumv - sumu + Math.Sqrt((sumv - sumu) ^ 2 + 4 * (sumuv) ^ 2)) / (2 * sumuv)
    b = (
        (variance_y - variance_x)
        + math.sqrt((variance_y - variance_x) ** 2 + 4 * (covariance) ** 2)
    ) / (2 * covariance)
    rb = covariance / (math.sqrt(variance_x * variance_y))
    db = (b / rb) * (math.sqrt((1 - rb**2) / len(x)))
    # dint = meanY - (db * meanX)
    print(f"\nGradient of line: {grad} \u00b1 {db}")
    # print("\n Gradient of line: {0}, {1}".format(b, db))
    print(f"Intercept of line: {inte}")

    return grad, inte, db


def peak_diff(t1, t2):
    """Take difference in time and histrogram for a peak"""
    diff_list = []
    for time1, time2 in zip(t1, t2):
        diff_list.append(time1 - time2)

    # plt.hist(diff_list)
    return diff_list


def get_1D_ranges(dlist, padding):
    """get the max and min of a 1D list to allow setting of histogram range

    Args:
        dlist (float): list of floats
        padding (int): Any padding to add on either end

    Returns:
        _type_: _description_
    """

    binmin = min(dlist) - padding
    binmax = max(dlist) + padding

    return binmin, binmax


def get_1D_histogram(dlist, padding, binsize):
    """Take 1D list of time differences
    returns the histogram of this

    Args:
        dlist (float): list of floats
        padding (int): Any padding to add on either end
        binsize (int): Size of histogram bin

    Returns:
        _type_: _description_
    """
    binmin, binmax = get_1D_ranges(dlist, padding)
    hist_x = np.arange(binmin, binmax + binsize, binsize)
    # get histogram of time differences
    peakd = np.histogram(dlist, bins=hist_x)

    return peakd

def calc_time_diff_hist(t1, t2, bin_size=1):
    # Compute time differences
    diff_list = [time1 - time2 for time1, time2 in zip(t1, t2)]
    
    if not diff_list:
        raise ValueError("Error: No time differences to process.")

    # Get min and max values
    bin_min, bin_max = min(diff_list), max(diff_list)

    # Ensure bins include 0
    bin_min = min(bin_min, 0)

    # Define bin edges with step of 1
    bins = np.arange(bin_min, bin_max + bin_size, bin_size)

    # Compute histogram
    hist_values, bin_edges = np.histogram(diff_list, bins=bins)

    return hist_values, bin_edges

def save_hist_csv(hist_values, bin_edges, filename="time_diff_histogram.csv", folder="."):
    os.makedirs(folder,exist_ok=True)

    file_path = os.path.join(folder, filename)

    # Save to CSV
    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Bin Edge", "Count"])  # Column headers
        writer.writerows(zip(bin_edges[:-1], hist_values))

    print(f"Histogram saved to {file_path}")

def save_time_diff_csv(t1, t2, filename="time_diff_histogram.csv", folder=".", bin_size=1):
    """Compute time differences, bin them in bins of 1, and save as CSV.

    Args:
        t1 (list of float): First time series.
        t2 (list of float): Second time series.
        filename (str): Name of the output CSV file.
        bin_size (int): Bin width for histogram (set to 2).
    """
    # Compute time differences
    diff_list = [time1 - time2 for time1, time2 in zip(t1, t2)]
    
    if not diff_list:
        raise ValueError("Error: No time differences to process.")

    # Get min and max values
    bin_min, bin_max = min(diff_list), max(diff_list)

    # Ensure bins include 0
    bin_min = min(bin_min, 0)

    # Define bin edges with step of 1
    bins = np.arange(bin_min, bin_max + bin_size, bin_size)

    # Compute histogram
    hist_values, bin_edges = np.histogram(diff_list, bins=bins)

    os.makedirs(folder,exist_ok=True)

    file_path = os.path.join(folder, filename)

    # Save to CSV
    with open(file_path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Bin Edge", "Count"])  # Column headers
        writer.writerows(zip(bin_edges[:-1], hist_values))

    print(f"Histogram saved to {file_path}")

def dead_time_padding(hist_values, start_idx: int, steps: int):
    #Check that length of list is not larger than the number of steps to be averaged
    end_idx = min(start_idx + steps, len(hist_values))
    
    val_sum = sum(hist_values[start_idx:end_idx])
    val_avg = val_sum / (end_idx - start_idx) if end_idx > start_idx else 0


    for j in range(start_idx - 1, -1, -1):
        hist_values[j] = val_avg

    return hist_values
    

def get_fwhm(ylist, xlist, method):
    """calculate fwhm of the 1D histogram


    Args:
        ylist (float): y axis of histogram
        xlist (float): x axis of histogram
        method (string): name of method to use to get FWHM
    """

    if method == "average":
        hei = np.average(ylist)
        halfmax = hei / 2
        crosses_half = np.where(ylist >= halfmax)[0]
        fwhm = xlist[crosses_half[-1]] - xlist[crosses_half[0]]
    else:
        fwhm = -1

    return fwhm


def generate_diff(blist, tmin, tmax):
    """Take the 2D time histogram and produce 1D diff histogram
     and centre of mass of peak

    Args:
        blist (float): output from np.2dhist
        tmin (_type_): min of difference spectrum
        tmax (_type_): max of difference
    """
    lrange = tmax - tmin
    inten = np.zeros(lrange)  # define as just big for now
    # iterate over the 2D array to create difference plot
    for index1, xv in enumerate(blist[1][:-1]):
        xv = int(xv)
        for index2, yv in enumerate(blist[2][:-1]):
            yv = int(yv)
            index = int(yv - xv)
            value = blist[0][index1][index2]
            inten[index] = inten[index] + value

    return inten


def generate_diff2(blist):
    """Take the 2D time histogram and produce 1D diff histogram
     and centre of mass of peak

    Args:
        blist (float): output from np.2dhist

    """
    diff_list = []  # define as just big for now
    # iterate over the 2D array to create difference plot
    for index1, xv in enumerate(blist[1][:-1]):
        xv = int(xv)
        for index2, yv in enumerate(blist[2][:-1]):
            yv = int(yv)
            index = yv - xv
            value = int(blist[0][index1][index2])
            for _ in range(0, value):
                diff_list.append(index)

    return diff_list


def generate_proj(blist):
    """get x projection and counts

    Args:
        blist (_type_): _description_
    """
    pro_list = []
    # iterate over the 2D array to create difference plot
    for index1, xv in enumerate(blist[1][:-1]):
        xv = int(xv)
        for index2, yv in enumerate(blist[2][:-1]):
            yv = int(yv)
            index = xv
            value = int(blist[0][index1][index2])
            for _ in range(0, value):
                pro_list.append(index)

    return pro_list


def fit_unigauss(x, amp, a, xc, sigma):
    """Returns a convolution of a gaussian and a square topped peak

    Args:
        x (float): postion on x axis
        amp (float): the height of the peak
        a (float): width of the flat topped peak
        xc (float): centre of peak
        sigma (float): width of gaussian
    """
    uni = (
        amp
        * (1 / (4 * a))
        * (
            math.erf(((x - xc) + a) / (sigma * math.sqrt(2)))
            - math.erf(((x - xc) - a) / (sigma * math.sqrt(2)))
        )
    )

    return uni


def create_rot_mat(radians):
    """Create a rotation matrix

    Args:
        radians (float): Angle in radians to perform rotation
    Returns:
        cos_rad (float): cos(Theta)
        sin_rad (float): sin(Theta)
    """
    # x, y = xy
    # offset_x, offset_y = origin
    # adjusted_x = (x - offset_x)
    # adjusted_y = (y - offset_y)
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    # qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    # qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y

    return cos_rad, sin_rad


def apply_rot_mat(blist, trig_val, origin=(0, 0)):
    """Function to apply rotation matrix to 2D histogram

    Args:
        blist  (float): output from np.2dhist
        trig_val (float): cos and sin theta for rotation
        origin (tuple, optional): Origin of rotation. Defaults to (0, 0).

    Returns:
       rlist (float): the rotated output
    """
    # create zeroed array to fill
    rlist = np.zeros_like(blist[0])
    # get origin
    offset_x, offset_y = origin
    # iterate over the 2D array to create difference plot
    for index1, xv in enumerate(blist[1][:-1]):
        xv = int(xv)
        for index2, yv in enumerate(blist[2][:-1]):
            yv = int(yv)
            # shift coordinates based on origin
            adjusted_x = index1 - offset_x
            adjusted_y = index2 - offset_y
            # apply rotation - convert to int as will be index
            qx = int(offset_x + trig_val[0] * adjusted_x + trig_val[1] * adjusted_y)
            qy = int(offset_y + -trig_val[1] * adjusted_x + trig_val[0] * adjusted_y)
            value = blist[0][index1][index2]
            if (
                (qx >= 0)
                and (qx < len(blist[1][:-1]))
                and (qy >= 0)
                and (qy < len(blist[2][:-1]))
            ):
                rlist[qx][qy] = rlist[qx][qy] + value
    return (rlist, blist[1], blist[2])


def get_sigma_p(xlist, ylist, fcentre, xpad=200, ypad=200, binsize=1):
    """function to calculate sigma p for transformed pairs data

    Args:
        xlist (_float_): t1 values for the peak
        ylist (_float_): t2 values for the peak
        fcentre (_Bool_): Whether to find centre or not
        xpad (int): Padding for x in binning (200 is default)
        ypad (int): Padding for y in binning (200 is default)
        binsize (int):  Size of histogram bins (1 is default)

    Returns:
        _type_: _description_
        sig_p_d (_float_): Value of sigma_d for this data
        xAcen (_float_): X Value of the centre of the transformed peak
        yAcen (_float_): Y Value of the centre of the transformed peak
    """

    # get the ranges for the histogram
    bint1, bint2 = get_bin_ranges(xlist, ylist, xpad, ypad, binsize)
    peakd = np.histogram2d(xlist, ylist, bins=(bint1, bint2))
    # check sigma P
    # have to work out if h(x,y) =h(-x,-y)
    # If we assume that the CoM of the peak gives x=0 and y=0 then scale indices accordingly
    # h(x,y) - h(-x,-y)

    if fcentre:
        # get CoM
        cm = center_of_mass(peakd[0])
        # Turn into integers
        xcen = int(cm[0])
        ycen = int(cm[1])
        xAcen = peakd[1][int(cm[0])]
        yAcen = peakd[2][int(cm[1])]
    else:
        xcen = 0
        ycen = 0
        xAcen = 0
        yAcen = 0

    # Zero sum and sigma
    sump = 0
    sig_p_d = 0
    sig_bad = 0
    # iterate through the histrogram
    badp = 0
    badhist = []
    for indexX, __ in enumerate(peakd[1][:-1]):
        for indexY, __ in enumerate(peakd[2][:-1]):
            # get number of counts, N
            sump = sump + peakd[0][indexX][indexY]
            # get opposite indices
            # process is index - cent to shift into centred coordinates
            # negate this to get opposite index in centred coordinates
            # then + cent to return to index of list
            # therefore cent-(index-cent), or 2*cent-index
            negIndX = 2 * xcen - indexX
            negIndY = 2 * ycen - indexY
            if peakd[0][indexX][indexY] > 0:
                # check it is a non-zero value
                if (negIndX >= len(peakd[0])) or (negIndY >= len(peakd[0][0])):
                    badp = badp + 1
                    sig_bad = sig_bad + (
                        peakd[0][indexX][indexY] * peakd[0][indexX][indexY]
                    )
                    badhist.append(
                        [peakd[0][indexX][indexY], peakd[1][indexX], peakd[2][indexY]]
                    )
                else:
                    try:
                        sig_p_d = sig_p_d + (
                            (peakd[0][indexX][indexY]) - (peakd[0][negIndX][negIndY])
                        ) * ((peakd[0][indexX][indexY]) - (peakd[0][negIndX][negIndY]))
                    except IndexError:
                        print(
                            f"{len(peakd[0])},{len(peakd[1])},{len(peakd[2])},{indexX},{indexY},{negIndX},{negIndY}"
                        )

    sig_p_d = sig_p_d / sump
    # think we are double counting?

    return sig_p_d, xAcen, yAcen, badp, badhist


def get_sigma_d(xlist, ylist, fcentre, kin, kfin, knum, xpad=300, ypad=300, binsize=1):
    """Calculates sigma d which checks the symmetry with the following transform
    h(x,y) =  h(k*y,x/k)

    Args:
        xlist (float): List of the x values to use
        ylist (float): List of the y values to use
        fcentre (bool): Whether to find the centre of the peak (should be true normally)
        kin (float): Initial momentum ratio
        kfin (float): Final momentum ratio
        knum (int): Number of momentum ratios to calculate
        xpad (int, optional): Padding in the x dimension to allow transformation to not go out of bounds. Defaults to 300.
        ypad (int, optional): Padding in the y dimension to allow transformation to not go out of bounts. Defaults to 300.
        binsize (int, optional): Bin size of histogram. Defaults to 1.

    Returns:
        : A mixture ,includes calculate sigma d, the x and y centres of the peak, the list of momenta,how many bad points and sum of bad points
    """
    bint1, bint2 = get_bin_ranges(xlist, ylist, xpad, ypad, binsize)
    peakd = np.histogram2d(xlist, ylist, bins=(bint1, bint2))
    # check sigma d
    # have to work out if h(x,y) =h(ky,x/k)
    # Where k = p/q
    # If we assume that the CoM of the peak gives x=0 and y=0 then scale indices accordingly
    if fcentre:
        # get CoM
        cm = center_of_mass(peakd[0])
        # Turn into integers
        xcen = int(cm[0])
        ycen = int(cm[1])
        xAcen = peakd[1][int(cm[0])]
        yAcen = peakd[2][int(cm[1])]
    else:
        xcen = 0
        ycen = 0
        xAcen = 0
        yAcen = 0

    # Zero sum and sigma
    sump = 0
    sig_p_d = 0
    sig_bad = 0
    # iterate through the histrogram for a given k
    # kval = kin
    sig = []
    ks = []
    sbad = []
    sigy_bad = []
    for kval in np.linspace(kin, kfin, num=knum):
        sig_p_d = 0
        sum_bad = 0
        sump = 0
        for indexX, __ in enumerate(peakd[1][:-1]):
            for indexY, __ in enumerate(peakd[2][:-1]):
                if peakd[0][indexX][indexY] > 0:
                    # get number of counts, N
                    sump = sump + peakd[0][indexX][indexY]
                    # get scaled indices
                    newIndX = int(((indexY - ycen) / kval) + xcen)
                    newIndY = int(((indexX - xcen) * kval) + ycen)
                    if (
                        (newIndX < len(peakd[1][:-1]))
                        and (newIndX > 0)
                        and (newIndY < len(peakd[2][:-1]))
                        and (newIndY > 0)
                    ):
                        sig_p_d = sig_p_d + (
                            (peakd[0][indexX][indexY]) - (peakd[0][newIndX][newIndY])
                        ) * ((peakd[0][indexX][indexY]) - (peakd[0][newIndX][newIndY]))
                        # print(f"Old: {peakd[1][indexX]} , {peakd[2][indexY]}")
                        # print(f"New: {peakd[1][newIndX]} , {peakd[2][newIndY]}")
                    else:
                        sig_bad = sig_bad + (
                            peakd[0][indexX][indexY] * peakd[0][indexX][indexY]
                        )
                        sum_bad = sum_bad + 1

        sig_p_d = sig_p_d / sump
        sig.append(sig_p_d)
        ks.append(kval)
        sbad.append(sum_bad)
        sigy_bad.append(sig_bad)
    # think we are double counting?
    return sig, xAcen, yAcen, ks, sbad, sigy_bad


def get_bin_ranges(xlist, ylist, xbuff, ybuff, bintsize):
    """Sets up ranges for histograms
    buffers are to pad the histogram

    Args:
        xlist (list): times to be binned
        ylist (list): times to be binned
        xbuff (int): padding in x
        ybuff (int): padding in y
        bintsize (int): siz of the bins
    """
    # set the ranges based on data
    bt1max = math.ceil(max(xlist)) + xbuff
    bt1min = math.floor(min(xlist)) - xbuff
    bt1range = bt1max - bt1min
    bt2max = math.ceil(max(ylist)) + ybuff
    bt2min = math.floor(min(ylist)) - ybuff
    bt2range = bt2max - bt2min
    # set the binsize
    # bintsize = 1
    # generate bins to avoid numpy not being integer
    bint1 = np.arange(bt1min, bt1max, step=bintsize)
    bint2 = np.arange(bt2min, bt2max, step=bintsize)

    return bint1, bint2

def filter_triples(results, t1, t2):
    tuples = []
    
    for triplet in results.get("triples"):  # Iterate through triples
        a, b, c = triplet  # Unpack the values
        
        # Check if at least one integer is in range [t1, t2]
        if t1 <= a <= t2 or t1 <= b <= t2 or t1 <= c <= t2:
            # Add the remaining two integers to the list
            tuples.append((b, c)) if t1 <= a <= t2 else \
            tuples.append((a, c)) if t1 <= b <= t2 else \
            tuples.append((a, b))
    
    return {"pairs": tuples}  # Return a dictionary so that select_pairs can be used

def select_triples(results, t1, t2, box_coord, box_shape=False):
    """Function to select pairs that fall within a box region

    Parameters
    ----------
    results : dict
        Results loaded from a eiei data set, must contain a dictionary of triples
    t1, t2 : int
        lower and upper time values for singles peak boundaries 
    box_coord : list
        list of the box coordinates to make selection in correct units, is in the form
        box_coord[0]: x1
        box_coord[1]: y1
        box_coord[2]: x2
        box_coord[3]: y2
        box_coord[4]: gradient

    box_shape : bool, optional
        If true is a flat topped box, by default False

    Returns
    -------
    lists
        Selected pairs, total counts and the averages inside the box
    """

    pairs = filter_triples(results, t1, t2)
    spairs, counts, means = select_pairs(pairs, box_coord, box_shape)
    return spairs, counts, means

def calc_shift_amount(t1,t2,t1_ref,t2_ref):
    val1 = t1 - t1_ref
    val2 = t2 - t2_ref
    if val1 * val2 < 0:
        return 0
    avg = (val1 + val2) / 2
    if avg >= 0:
        return int(math.floor(avg))
    else:
        return int(math.ceil(avg))
    

