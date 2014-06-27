#
# Support functions for plotting histograms
#
import math, os, re, sys
import cogent.maths.stats.test as stats

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

def set_params(data, fixed_max=None):
    """Determines min, max, and bar size based on the input data."""

    if fixed_max != None:
        maxval = fixed_max
    else:
        maxval = None  
        if maxval == None:
            maxval = max(data)
        if max(data) > maxval:
            maxval = max(data)
    
    minval = None  
    if minval == None:
        minval = min(data)
    if min(data) < minval:
        minval = min(data)
    
    plot_params["NBINS"] = 20
    plot_params["MIN_BIN_FLOOR"] = minval - (minval)
    plot_params["MIN_BIN_CEILING"] = maxval + (1.01 * maxval)
    plot_params["BIN_SIZE"] = (plot_params["MIN_BIN_CEILING"] - plot_params["MIN_BIN_FLOOR"]) / plot_params["NBINS"]
    
    print "\n. Dynamic plot parameters: minimum value=", plot_params["MIN_BIN_FLOOR"], "max. value=", plot_params["MIN_BIN_CEILING"], "bin width=", plot_params["BIN_SIZE"]

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
        
    tablepath = "barplot.table." + filekeyword + ".txt"
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

    pdfpath = "barplot." + filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"

    # X axis labels:
    cranstr += "bins <- c("
    for i in range(0, plot_params["NBINS"]):
        cranstr += bins[i].__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, ylim=range(0,1.0), names.arg=bins, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\");\n"
    cranpath = "barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)

def plot_histogram(data, keyword, xlab="", ylab="", fixed_max=None):
    """data is an array of floating point values"""
    plot_params = set_params(data, fixed_max=fixed_max)
    barplot1(data, "histogram." + keyword, plot_params, xlab=xlab, ylab=ylab)
    
    
    
    