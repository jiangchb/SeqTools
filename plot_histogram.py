#
# Support functions for plotting histograms
#
import math, os, re, sys

# depricated
#import cogent.maths.stats.test as stats

plot_params = {} # a global hashtable with parameters about the plot

# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )

def roundup(x, step):
    """Rounds the value x up the nearest step, where step could, e.g., be 100 or 1000"""
    return int(math.ceil(x / float(step))) * step

def rounddown(x, step):
    return int(math.floor(x / float(step))) * step

def set_params(data, fixed_max=None, fixed_nbins=None):
    """Determines min, max, and bar size based on the input data."""

    if fixed_max != None:
        maxval = float(fixed_max)
    else:
        maxval = max(data)
    
    minval = min(data)

    if fixed_nbins != None:
        plot_params["NBINS"] = int( fixed_nbins )
    else:
        plot_params["NBINS"] = 20
    plot_params["MIN_BIN_FLOOR"] = minval
    if plot_params["MIN_BIN_FLOOR"] > 0.0:
        plot_params["MIN_BIN_FLOOR"] = 0.0
    plot_params["MAX_BIN_CEILING"] = (1.01 * maxval)

    val_range = plot_params["MAX_BIN_CEILING"] - plot_params["MIN_BIN_FLOOR"]
    
    plot_params["BIN_SIZE"] = val_range / plot_params["NBINS"]
    
    # Now adjust the bin size to be a nice round number
    if plot_params["BIN_SIZE"] <= 10:
        plot_params["BIN_SIZE"] = roundup(plot_params["BIN_SIZE"], 1)
        plot_params["MIN_BIN_FLOOR"] = rounddown( plot_params["MIN_BIN_FLOOR"], 1)
        plot_params["MAX_BIN_CEILING"] = roundup( plot_params["MAX_BIN_CEILING"], 1)
    elif plot_params["BIN_SIZE"] < 100:
        plot_params["BIN_SIZE"] = roundup(plot_params["BIN_SIZE"], 10)
        plot_params["MIN_BIN_FLOOR"] = rounddown( plot_params["MIN_BIN_FLOOR"], 10)
        plot_params["MAX_BIN_CEILING"] = roundup( plot_params["MAX_BIN_CEILING"], 10)
    elif plot_params["BIN_SIZE"] < 1000:
        plot_params["BIN_SIZE"] = roundup(plot_params["BIN_SIZE"], 100)
        plot_params["MIN_BIN_FLOOR"] = rounddown( plot_params["MIN_BIN_FLOOR"], 100)
        plot_params["MAX_BIN_CEILING"] = roundup( plot_params["MAX_BIN_CEILING"], 100)
    elif plot_params["BIN_SIZE"] < 10000:
        plot_params["BIN_SIZE"] = roundup(plot_params["BIN_SIZE"], 500)
        plot_params["MIN_BIN_FLOOR"] = rounddown( plot_params["MIN_BIN_FLOOR"], 500)
        plot_params["MAX_BIN_CEILING"] = roundup( plot_params["MAX_BIN_CEILING"], 500)
    
    n_neg_bins = int( roundup( abs(plot_params["MIN_BIN_FLOOR"]) / plot_params["BIN_SIZE"], 1.0) )
    n_pos_bins = int( roundup( abs(plot_params["MAX_BIN_CEILING"]) / plot_params["BIN_SIZE"], 1.0) )
    plot_params["NBINS"] = n_neg_bins + n_pos_bins
    
    if plot_params["MIN_BIN_FLOOR"] < 0.0:
        plot_params["MIN_BIN_FLOOR"] = n_neg_bins * plot_params["BIN_SIZE"] * -1
    plot_params["MAX_BIN_CEILING"] = n_pos_bins * plot_params["BIN_SIZE"]
    
    print "\n. Dynamic plot parameters: N bins=",plot_params["NBINS"], "minimum value=", plot_params["MIN_BIN_FLOOR"], "max. value=", plot_params["MAX_BIN_CEILING"], "bin width=", plot_params["BIN_SIZE"]

    return plot_params

def barplot1(values, filekeyword, plot_params, xlab="", ylab=""):    
    bins = [] # a list of bin floors
    for i in range(0, plot_params["NBINS"]):
        bins.append(plot_params["MIN_BIN_FLOOR"] + i*plot_params["BIN_SIZE"] )
     
    bin_values = {} # key = bin index, value = count of values in that bin
    for i in range(0, plot_params["NBINS"]):
        bin_values[i] = 0
 
    for val in values:
        for i in range(0, plot_params["NBINS"]):    
            if val >= bins[i]:
                if i == plot_params["NBINS"]-1: # last bin
                    bin_values[i] += 1
                elif val < bins[i+1]:
                    bin_values[i] += 1

    normalized_bin_values = {}
    for i in bin_values:
        normalized_bin_values[i] = float(bin_values[i]) / values.__len__()
        
    tablepath = filekeyword + ".table.txt"
    fout = open(tablepath, "w")
    for i in range(0, plot_params["NBINS"]):
        if i < plot_params["NBINS"]-1:
            fout.write(bins[i].__str__() + "\t")
        else:
            fout.write(bins[i].__str__() )
    fout.write("\n")
    for i in range(0, plot_params["NBINS"]):
        if i < plot_params["NBINS"]-1:
            fout.write( normalized_bin_values[i].__str__() + "\t")
        else:
            fout.write( normalized_bin_values[i].__str__() )            
    fout.write("\n")
    fout.close()

    pdfpath = filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"

    # X axis labels:
    cranstr += "bins <- c("
    for i in range(0, plot_params["NBINS"]):
        cranstr += "'" + bins[i].__str__() 
        if i == plot_params["NBINS"]-1:
            cranstr += "+"
        cranstr += "',"
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars),axes = FALSE, axisnames = FALSE, beside=TRUE, ylim=range(0,1.0), xlab=\"" + xlab + "\", ylab=\"" + ylab + "\");\n"
    cranstr += "text(barx, par(\"usr\")[3]-0.05, labels=bins, srt=45, adj=1, xpd=TRUE, cex=0.9)\n"
    #cranstr += "axis(1, at=barx, labels = FALSE)\n"
    cranstr += "axis(2)\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)

def plot_histogram(data, keyword, xlab="", ylab="", fixed_max=None, fixed_nbins=None):
    """data is an array of floating point values"""
    plot_params = set_params(data, fixed_max=fixed_max, fixed_nbins=fixed_nbins)
    barplot1(data, keyword + ".histo", plot_params, xlab=xlab, ylab=ylab)
    