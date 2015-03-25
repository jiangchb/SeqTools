import math, os, re, sys

from scipy import stats as scipystats


def matplot_scatter1(filekeyword, xvalues, yvalues, names=None, xmin=None, ymin=None, xmax=None, ymax=None, logy=False, logx=False, xlab="", ylab="", format="pdf", colors=["blue", "red", "purple", "black"], sizes=[10,18,18,18], markers=["o", "D", "v", "<"], alphas=[0.4,0.4,0.4,0.4] ):
    """xvalues and yvalues are lists of lists. Each list is series of points.
        colors is a list of colors, used to color series 0, series 1, series 2, etc."""
    
    """Sanity Check"""
    if xvalues.__len__() != yvalues.__len__():
        print "\n. Error 11: there are different numbers of X and Y data series."
        print xvalues.__len__(), "vs.", yvalues.__len__()
        exit()
        
    for ii in range(0, xvalues.__len__()):
        if xvalues[ii].__len__() != yvalues[ii].__len__():
            print "\n. Error 18, The x and y series have different lengths."
            print xvalues[ii].__len__(), yvalues[ii].__len__()
            exit()
    
    if colors.__len__() < xvalues.__len__():
        pass
    if sizes.__len__() < xvalues.__len__():
        pass
    if alphas.__len__() < xvalues.__len__():
        pass
    
    if names == None:
        for ii in range(0, xvalues.__len__()):
            names.append("Series " + ii.__str__())
    
    """If the user didn't specify xmin and ymin, then let's
        auto-determine these values from the first data series."""
    if xmin == None:
        xmin = min( map(min, xvalues) )
    if ymin == None:
        ymin = min (map (min, yvalues) )
    if xmax == None:
        xmax = max (map (max, xvalues) )
    if ymax == None:
        ymax = max (map( max, yvalues) )
    
    import numpy as np
    import scipy as sp
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from matplotlib.gridspec import GridSpec
    import pylab as P
        
    """Setup the figure"""
    fig1 = plt.figure(figsize=[8,6])
    gs = GridSpec(100, 100, bottom=0.10,left=0.15,right=0.85)
    ax1 = fig1.add_subplot(gs[15:100,0:85])
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    if logx:
        ax1.set_xscale("log")
    if logy:
        ax1.set_yscale("log")
    ax1.set_xlabel(xlab, fontsize=11)
    ax1.set_ylabel(ylab, fontsize=11)
    ax1.set_xticks([xmin,xmax])
    ax1.set_yticks([ymin,ymax])
    ax1.set_xticklabels([xmin.__str__(), xmax.__str__()])
    ax1.set_yticklabels([ymin.__str__(), ymax.__str__()] )
    ax1.tick_params(axis='x',length=3,width=1.5)
    
    """Plot the data"""
    returned_plots = []
    for ii in range(0, xvalues.__len__()):
        h = ax1.scatter(xvalues[ii], yvalues[ii], color=colors[ii], alpha=alphas[ii], edgecolors='none', s=sizes[ii], marker=markers[ii])
        returned_plots.append(h)

    """Figure Legend"""
    plt.legend(returned_plots,
               names,
               loc='upper left',
               ncol=1,
               fontsize=9, fancybox=True, shadow=True)

    """Top histogram"""
    ax2 = fig1.add_subplot(gs[0:14,0:85])
    ax2.set_xlim(xmin,xmax)
    if logx:
        n, bins, patches = ax2.hist(xvalues[0], bins = 10 ** np.linspace(np.log10(xmin), np.log10(xmax), 50), color="#DDDDDD", edgecolor = "none")
        ax2.set_xscale("log")
    else:
        n, bins, patches = ax2.hist(xvalues[0], bins = np.linspace(xmin, xmax, 50), color="#DDDDDD", edgecolor = "none")
    
    """Right-side histogram"""    
    ax3 = fig1.add_subplot(gs[15:100,86:100])
    ax3.set_ylim(ymin,ymax)
    if logy:
        n, bins, patches = ax3.hist(yvalues[0], orientation="horizontal", bins = 10 ** np.linspace(np.log10(ymin), np.log10(ymax), 50), color="#DDDDDD", edgecolor = "none")
        ax3.set_yscale("log")
    else:
        n, bins, patches = ax3.hist(yvalues[0], orientation="horizontal", bins = np.linspace(ymin, ymax, 50), color="#DDDDDD", edgecolor = "none")

    """Labels for the histograms"""
    ax2.set_frame_on(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)

    ax3.set_frame_on(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    
    """Finally, save the image."""
    fig1.savefig(filekeyword + "." + format, format=format, dpi=600)

def scatter1(values_ii, values_jj, filekeyword, xlab="", ylab="", force_square=False):    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    cranstr += "pdf(\"" + pdfpath + "\", width=6, height=6);\n"    

    print "\n. Writing a scatterplot to", pdfpath

    # X values
    cranstr += "x<-c("
    for v in values_ii:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"

    # Y values
    cranstr += "y<-c("
    for v in values_jj:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "plot(x, y, cex=0.5, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""

    maxx = max(values_ii)
    maxy = max(values_jj)
    lim = max( [maxx, maxy] )    
    if force_square:
        cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
    
    cranstr += ");\n"
    
    if force_square:
        cranstr += "abline(0,1)\n"
    
    """Pearson's linear value correlation."""
    corr_valsa = []
    corr_valsb = []
    for ww in range(0, values_ii.__len__()):
        if values_ii[ww] != 0 and values_jj[ww] != 0:
            corr_valsa.append( values_ii[ww] )
            corr_valsb.append( values_jj[ww] )
    (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + lim.__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
    
    """Spearman's non-linear non-parametric rank correlation."""
    (rho, pvalue) = scipystats.spearmanr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + (0.92*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
    
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save --slave < " + cranpath)
    
    return cranpath

def scatter1xn(values, filekeyword, title="", xlab="", ylab="", force_square=False):
    """Makes a multi-panel scatterplot, with 1 row and N columns."""
    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    """values[set name][0] = X value list, values[set name][1] = y value list."""
    pdfpath = filekeyword + ".pdf"
    cranstr += "pdf(\"" + pdfpath + "\", width=" + (3*values.keys()).__len__().__str__() + ", height=3);\n"    

    print "\n. Writing a scatterplot to", pdfpath

    cranstr += "par(mar=c(1.8,1.8,0.4,0.4), oma=c(1.5,2,1,1)  );\n"

    colwidth = 1.0 / values.__len__()

    sets = values.keys()
    for ii in range(0, sets.__len__()):
        setname = sets[ii]
        values_ii = values[setname][0]
        values_jj = values[setname][1]

        cranstr += "par( fig=c(" + (ii*colwidth).__str__() + "," + ((ii+1)*colwidth).__str__() + ", 0, 1)"
        if ii != 0:
            cranstr += ", new=TRUE"
        cranstr += ");\n"

        # X values
        cranstr += "x<-c("
        for v in values_ii:
            cranstr += v.__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"
    
        # Y values
        cranstr += "y<-c("
        for v in values_jj:
            cranstr += v.__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"
        
        cranstr += "plot(x, y, cex=0.5, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""

        maxx = max(values_ii)
        maxy = max(values_jj)
        lim = max( [maxx, maxy] )        
        if force_square:
            cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
        
        cranstr += ");\n"
        
        if force_square:
            cranstr += "abline(0,1)\n"
        
        """Pearson's linear value correlation."""
        corr_valsa = []
        corr_valsb = []
        for ww in range(0, values_ii.__len__()):
            if values_ii[ww] != 0 and values_jj[ww] != 0:
                corr_valsa.append( values_ii[ww] )
                corr_valsb.append( values_jj[ww] )
        (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
        cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + lim.__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
        
        """Spearman's non-linear non-parametric rank correlation."""
        (rho, pvalue) = scipystats.spearmanr( values_ii, values_jj )
        cranstr += "text(" + ((lim-min(values_ii))/2).__str__() + ", " + (0.92*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
    
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=1.8);\n"
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save --slave < " + cranpath)
    
    return cranpath

def scatter_idr_nxm(width, height, values, names, filekeyword, title="", xlab="", ylab="", force_square=True):    
    """
    This method generates a PDF grid of scatterplots. Each grid is height*height, and there 3 grids on each page.
    The leftmost grid is a scatterplot of the ranked values in values[0:height-1] compared to each other.
    The middle grid is the Psi(t) plots from the IDR analysis.
    The rightmost grid is the Psi'(t) plots.
    
    If width > height, then it must be a factor of height. For example, width can be 2xheight, 3xheight,
    but not 1.1xheight.
    For each |height| sized set of arrays in values, a PDF page with three grids (described above) will be
    plotted.

    values[ii] = list of data. There should be 'width' number of entries in values.
    
    NOTE: in making the scatterplots, paired values of (0,0) will be ignored from the plot.
    This is a decision that is unique to the IDR analysis.
    """  
        
    if width <= 1:
        print "Error plot_scatter.py 159: The method scatter_idr_nxm requires a width of at least 2."
        exit()
    
    if names.__len__() != width:
        print "\n. ERROR plot_scatter.py 427, you called scatter_nxm without enough names."
        print names
        print names.__len__()
        exit()
    
    for ii in range(0, width-1):
        if values[ii].__len__() != values[ii+1].__len__():
            print "\n. ERROR scatter_nxm was called with value arrays of unequal length."
            print names
            exit()
        if values[ii].__len__() == 0:
            print "\n. ERROR scatter_idr_nxm can't deal with empty value arrays, i.e., you have no data to plot."
            print "\n. Check your call of scatter_idr_nxm for " + filekeyword
            exit()
        
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    tablepaths = {} # key = tablepath for IDR output data, value = tuple(the name of the comparison, ii, jj)
    
    """Setup the document."""
    pdfpath = filekeyword + ".pdf"
    print "\n. Computing IDR and plotting to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=" + (12*height).__str__() + ", height=" + (3*height).__str__() + ");\n"    
    cranstr += "par(mar=c(1.8,2.8,2.8,1), oma=c(1.5,2,1,1)  );\n"  
    colwidth = 1.0 / (4.0 * float(height-1))
    if height == 1:
        colwidth = 0.25
    rowheight = 1.0 / float(height-1)
    if height == 1:
        rowheight = 1.0
    gridsize = ( (height-1)*(height-1) - (height-2)*(0.5*(height-1))) # the number of scatterplots in each grid.
    if height == 1:
        gridsize = 1.0
    
    
    """value_pairs is a list of pairs (ii,jj) where ii and jj are indices into values.
    The idea is that the values for each pair should be compared."""
    value_pairs = [] 
    for ii in range(0, width):
        mod = 0
        for qq in range(1, width):
            if qq%height == 0:
                if ii >= qq:
                    mod = qq
        for jj in range(mod+(ii%height), mod+height):
            if ii != jj:
                value_pairs.append( (ii,jj) )
    
    ii_jj_compname = {}
    ii_jj_idmap = {}
    for (ii,jj) in value_pairs:
        if ii not in ii_jj_compname:
            ii_jj_compname[ii] = {}
        ii_jj_compname[ii][jj] = names[ii] + "-" + names[jj]
        if ii not in ii_jj_idmap:
            ii_jj_idmap[ii] = {}
        ii_jj_idmap[ii][jj] = {}
    
    """This first loop is to determine total_count.
    total_count is used solely for displaying a progress bar,
    and is not essential to the algorithm's accuracy."""
    total_count = value_pairs.__len__()

    """Now we iterate through pairs of arrays in values."""
    for count in range(0, value_pairs.__len__()):
        #if count >= gridsize and count%gridsize == 0:
        #    cranstr += "par(mar=c(1.8,2.8,2.8,1), oma=c(1.5,2,1,1)  );\n"
        #    cranstr += "plot.new();\n"
        if ii == jj: # We can't do IDR on data that is exactly identical.
            continue

        """OK, Progress bar display"""
        sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
        sys.stdout.flush()

        """Extract the values we're going to work with"""        
        ii = value_pairs[count][0]
        jj = value_pairs[count][1]
        values_ii = values[ii]
        values_jj = values[jj]
        
        """Sanity Check:"""
        if values_ii.__len__() != values_jj.__len__():
            print "Error: different sizes for values_ii and values_jj. Mark point 225."
            exit()

        """Next we build the x and y array strings for the R script.
        Note that we reject all x,y pairs in the paired data where both x=0 and y=0."""
        r_ids2erg = {}
        #genen_rn = {} # key = gene index in the R script, gene index in the values

        """x and y"""
        genecount = 0
        xcranstr = "x<-c("
        ycranstr = "y<-c("
        """xcranstr and ycranstr are local string variables for just this inner loop."""
        for xx in range(0, values_ii.__len__()):
            if values_ii[xx] != 0 and values_jj[xx] != 0:
                xcranstr += values_ii[xx].__str__() + ","
                ycranstr += values_jj[xx].__str__() + ","
                """The gene with x and y value == genecount is actually
                the xx'th gene listed in values_ii and values_jj.
                This translation will be useful later,
                when we read IDR values from the R output."""
                ii_jj_idmap[ii][jj][genecount] = xx
                genecount += 1
        xcranstr = re.sub(",$", "", xcranstr)
        xcranstr += ");\n"
        ycranstr = re.sub(",$", "", ycranstr)
        ycranstr += ");\n"
        cranstr += xcranstr
        cranstr += ycranstr
        
        """rank"""
        cranstr += "rankx <- rank(-x);\n"
        cranstr += "ranky <- rank(-y);\n"
        cranstr += "lim <- max( max(rankx), max(ranky) );\n"
        
        """lim is the plotting range upper limit."""
        maxa = values_ii.__len__()
        maxb = values_jj.__len__()
        lim = max( [maxa, maxb] )    
        
        """Plot rank in the leftmost grid"""
        parstr = "par( fig=c(" + ((ii%(height))*colwidth).__str__() + ","
        parstr += ((ii%(height))*colwidth + colwidth).__str__() + ", "
        parstr += ( ( (jj%height - 1)%(height-1) )*rowheight).__str__() + "," 
        parstr += ( ( (jj%height - 1)%(height-1) )*rowheight + rowheight).__str__() + ")"
        if count > 0 and count%gridsize != 0:
            parstr += ", new=TRUE" # don't call new=TRUE for the first plot.
        parstr += ");\n"
        
        cranstr += parstr
        
        cranstr += "plot(rankx, ranky, cex=0.5, cex.lab=0.8, xlab=\"rank\", ylab=\"rank\""
        if force_square:
            cranstr += ", xlim=range(0,lim), ylim=range(0,lim)"
        col = "black"
        pch = "1"
        if ii < height:
            col = "mediumblue"
            pch = "5"
        elif ii < 2*height:
            col = "red3"
            pch = "5"
        elif ii < 3*height:
            col = "green4"
            pch = "5"
        elif ii < 4*height:
            col = "darkorchid3"
            pch = "5"
        cranstr += ", col=\"" + col + "\""
        cranstr += ", pch=" + pch
        cranstr += ", las=1"
        cranstr += ");\n"             
                    
        if force_square:
            cranstr += "abline(0,1)\n"
        
        """Leftmost grid: write custom write axis labels"""
        if (jj+1)%height==0:
            cranstr += "mtext(\"rank of " + names[ii] + "\", side=1, line=2, col=\"black\", cex=1);\n"
        if ii%height == 0:
            cranstr += "mtext(\"rank of " + names[jj] + "\", side=2, line=2.8, col=\"black\", cex=1);\n"
        
        """Use the idr library in R to compute IDR estimates."""
        cranstr += "library(idr);\n"
        cranstr += "rankx <- rank(x);\n"
        cranstr += "ranky <- rank(y);\n"
        cranstr += "mu <- 2.6;\n"
        cranstr += "sigma <- 1.3;\n"
        cranstr += "rho <- 0.8;\n"
        cranstr += "p <- 0.7;\n"
        cranstr += "idr.out <- est.IDR( cbind(x, y), mu, sigma, rho, p, eps=0.001);\n"
        cranstr += "uv <- get.correspondence(rankx, ranky, seq(0.01, 0.99, by=1/length(rankx) ) );\n"
        
        """Plot the psi values from IDR."""
        cranstr += "par( fig=c(" + (height*colwidth - colwidth + (ii%height)*colwidth).__str__() + ","
        cranstr += (height*colwidth + (ii%height)*colwidth).__str__() + ", "
        cranstr += ( ( (jj%(height) - 1)%(height-1) )*rowheight).__str__() + "," 
        cranstr += ( ( (jj%(height) - 1)%(height-1) )*rowheight + rowheight).__str__() + ")"
        if ii > 0 or jj > 0:
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
        cranstr += ");\n"
        
        cranstr += "plot(uv$psi.n$t, uv$psi.n$value, cex=0.5, xlab=\"t\", ylab=\"psi\", xlim=c(0, max(uv$psi.n$t)),ylim=c(0, max(uv$psi.n$value)), cex.lab=0.8"
        col = "black"
        pch = "1"
        if ii < height:
            col = "mediumblue"
            pch = "5"
        elif ii < 2*height:
            col = "red3"
            pch = "5"
        elif ii < 3*height:
            col = "green4"
            pch = "5"
        elif ii < 4*height:
            col = "darkorchid3"
            pch = "5"
        cranstr += ", col=\"" + col + "\""
        cranstr += ", pch=" + pch
        cranstr += ", las=1"
        #cranstr += ", main=\"" + ii.__str__() + ":" + jj.__str__() + "\""
        cranstr += ");\n"
        cranstr += "lines(uv$psi.n$smoothed.line, lwd=4);\n"
        cranstr += "abline(coef=c(0,1), lty=3);\n"
        
        if (jj+1)%height == 0 and ii==0:
            cranstr += "mtext(\"t\", side=1, line=2, col=\"black\", cex=1.2);\n"
            cranstr += "mtext(\"psi\", side=2, line=2.5, col=\"black\", cex=1.2);\n"
            """Write labels for top row"""
            cranstr += "mtext(\"IDR Correspondence Curve\", side=3, line=1, col=\"black\", cex=1.2);\n"
        
        """Plot the psi-prime values from IDR"""
        cranstr += "par( fig=c(" + ( (2*height*colwidth) - 2*colwidth + (ii%height)*colwidth).__str__() + ","
        cranstr += ( (2*height*colwidth) - colwidth + (ii%height)*colwidth).__str__() + ", "
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight).__str__() + "," 
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight + rowheight).__str__() + ")"
        
        if ii > 0 or jj > 0:
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
        cranstr += ");\n"
         
        cranstr += "plot(uv$dpsi.n$t, uv$dpsi.n$value, cex=0.5, xlab=\"t\", ylab=\"psi'\", xlim=c(0, max(uv$dpsi.n$t)),ylim=c(0, max(uv$dpsi.n$value)), cex.lab=1.0"
        col = "black"
        pch = "1"
        if ii < height:
            col = "mediumblue"
            pch = "5"
        elif ii < 2*height:
            col = "red3"
            pch = "5"
        elif ii < 3*height:
            col = "green4"
            pch = "5"
        elif ii < 4*height:
            col = "darkorchid3"
            pch = "5"
        cranstr += ", col=\"" + col + "\""
        cranstr += ", pch=" + pch
        cranstr += ", las=1"
        cranstr += ");\n"
        cranstr += "lines(uv$dpsi.n$smoothed.line, lwd=4);\n"
        cranstr += "abline(h=1, lty=3);\n"
        
        """Write custome axis labels."""
        if (jj+1)%height == 0 and ii==0:
            cranstr += "mtext(\"t\", side=1, line=2, col=\"black\", cex=1.2);\n"
            cranstr += "mtext(\"psi'\", side=2, line=2.5, col=\"black\", cex=1.2);\n"
            cranstr += "mtext(\"IDR Change of Correspondence Curve\", side=3, line=1, col=\"black\", cex=1.2);\n"
        
        if count == 0:
            cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=2.2);\n"
    
        """Plot the IDR values for the x replicate."""
        #cranstr += "par( fig=c(" + ( (3*height*colwidth) - 3*colwidth + (ii%height)*colwidth).__str__() + ","
        cranstr += "par( fig=c(" + (3*(height-1)*colwidth + (ii%height)*colwidth).__str__() + ","
        cranstr += (3*(height-1)*colwidth + (ii%height)*colwidth + colwidth).__str__() + ","
        #cranstr += ( (2*height*colwidth) - colwidth + (ii%height)*colwidth).__str__() + ", "
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight).__str__() + "," 
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight + 0.5*rowheight).__str__() + ")"
        if ii > 0 or jj > 0:
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
        cranstr += ");\n" 

        cranstr += "plot(x, idr.out$idr, cex=0.5, xlab=\"" + names[ii] + "\", ylab=\"idr\", xlim=c(0, max(x)),"
        #cranstr += "ylim=c(0, max(idr.out$idr)),"
        cranstr += "cex.lab=1.0"
        col = "black"
        pch = "1"
        if ii < height:
            col = "mediumblue"
            pch = "5"
        elif ii < 2*height:
            col = "red3"
            pch = "5"
        elif ii < 3*height:
            col = "green4"
            pch = "5"
        elif ii < 4*height:
            col = "darkorchid3"
            pch = "5"
        cranstr += ", col=\"" + col + "\""
        cranstr += ", pch=" + pch
        cranstr += ", las=1"
        #cranstr += ", log='y'"
        cranstr += ");\n"

        
        """Write custom axis labels."""
        if (jj+1)%height == 0 and ii==0:
            cranstr += "mtext(\"fold-enrichment\", side=1, line=2, col=\"black\", cex=1.2);\n"
            cranstr += "mtext(\"IDR\", side=2, line=2.5, col=\"black\", cex=1.2);\n"
            cranstr += "mtext(\"IDR vs. FE\", side=3, line=1, col=\"black\", cex=1.2);\n"
        
        if count == 0:
            cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=2.2);\n"

        #"""Plot the IDR values for the y replicate."""
        cranstr += "par( fig=c(" + (3*(height-1)*colwidth + (ii%height)*colwidth).__str__() + ","
        cranstr += (3*(height-1)*colwidth + (ii%height)*colwidth + colwidth).__str__() + ","
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight + 0.5*rowheight).__str__() + "," 
        cranstr += ( ( (jj%height - 1)%(height-1) )*rowheight + rowheight).__str__() + ")"
        if ii > 0 or jj > 0:
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
        cranstr += ");\n" 

        cranstr += "plot(y, idr.out$idr, cex=0.5, xlab=\"" + names[jj] + "\", ylab=\"idr\", xlim=c(0, max(x)),"
        #cranstr += "ylim=c(0, max(idr.out$idr)),"
        cranstr += "cex.lab=1.0"
        col = "black"
        pch = "1"
        if ii < height:
            col = "mediumblue"
            pch = "5"
        elif ii < 2*height:
            col = "red3"
            pch = "5"
        elif ii < 3*height:
            col = "green4"
            pch = "5"
        elif ii < 4*height:
            col = "darkorchid3"
            pch = "5"
        cranstr += ", col=\"" + col + "\""
        cranstr += ", pch=" + pch
        cranstr += ", las=1"
        #cranstr += ", log='y'"
        cranstr += ");\n"
        
        """Write the IDR data to a text file, so that Python can import the data."""
        tablepath = filekeyword + ".ii=" + ii.__str__() + ".jj=" + jj.__str__() + ".tmp"
        tablepaths[tablepath] = (ii_jj_compname[ii][jj], ii, jj)
        cranstr += "write.table( data.frame(idr.out$idr, idr.out$IDR), \"" + tablepath + "\", sep=\"\t\");\n"
    
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave --silent --quiet < " + cranpath)
    
    #print tablepaths
    #exit()
    
    idr_stats = read_idr_results(tablepaths, ii_jj_idmap)
    if idr_stats == False:
        print "\n. Error: the scatterplot with IDR data could not be created."
        print ". See previous error messages for more informaiton."
        exit()
    return (cranpath, sinkpath, idr_stats, value_pairs) 

def read_idr_results(tablepaths, ii_jj_idmap):
    """Read the IDR tables, save the data in idr_stats, and then destroy the tables.
    
    Returns idr_stats:
        idr_stats[gene number][ii][jj] = local IDR
        
    Return False if any errors occur.
    """
    idr_stats = {} # hash of hashes, key = gene N in values, value = hash, key = name of comparison, value = IDR for the gene
    compnames = []
    for tablepath in tablepaths:
        #print "447:", tablepath
        if False == os.path.exists(tablepath):
            print "\n. Error in read_idr_results:"
            print ". The R output table cannot be found at", tablepath
            return False
        compname = tablepaths[tablepath][0]
        ii = tablepaths[tablepath][1]
        jj = tablepaths[tablepath][2]
        #print "454:", compname, ii, jj
        idmap = ii_jj_idmap[ii][jj] # idmap[gene ID in R script] = gene ID in the genes array
        if compname not in compnames:
            compnames.append( compname )
        fin = open(tablepath, "r")
        for line in fin.xreadlines():
            if line.__len__() > 2:
                tokens = line.split()
                if tokens.__len__() == 3:
                    genen = int( re.sub("\"", "", tokens[0]) ) - 1
                    genen = idmap[ genen ]
                    lidr = float(tokens[1])
                    idr = float(tokens[2])
                    if genen not in idr_stats:
                        idr_stats[ genen ] = {}
                    if ii not in idr_stats[genen]:
                        idr_stats[genen][ii] = {}
                    idr_stats[genen][ii][jj] = lidr
        fin.close()
        os.system("rm " + tablepath)
    return idr_stats 

def scatter_nxm(width, height, values, names, filekeyword, title="", force_square=False, plot_as_rank = [], skip_identity = False, skip_zeros = False, unit_labels=[], xlab=None, ylab=None):    
    """
    Creates a multi-panel collection of scatterplots. The dimensions are N scatterplots by M scatterplots.
    n = width, number of scatterplots
    m = height, number of scatterplots.
    
    values[ii] = list of data for data series ii. values.__len__() should equal width
    
    skip_zeros: If True, then any values that equal zero will be removed from the analysis.
        This essentially allows for data from species-specific or replicate-specific genes
        to be removed.
    
    NOTE: in making the scatterplots, data pairs where both the x-axis and y-axis series are (0,0) 
    will be ignored from the plot.
    
    plot_as_rank is a list of indices into the values array which should be plotted as rank, rather than their raw value.
    
    If skip_identify = True, then the multi-panel plot will not include same-vs.-same scatterplots, and the dimensions
    will be updated accordingly.
    """    
    if names.__len__() != width:
        print "\n. ERROR scatter_nxm was called with 'width' != to the number of names."
        print names
        print names.__len__()
        exit()
    
    for ii in range(0, width-1):
        if values[ii].__len__() != values[ii+1].__len__():
            print "\n. ERROR scatter_nxm was called with value arrays of unequal length."
            print names
            exit()
        if values[ii].__len__() == 0:
            print "\n. ERROR scatter_nxm can't deal with empty value arrays, i.e., you have no data to plot."
            print "\n. Check your call of scatter_nxm for " + filekeyword
            return None
    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    
    if skip_zeros:
        no_zero_values = []
        for ii in range(0, width-1):
            no_zero_values.append( [] )
        
        """For each gene"""
        for jj in range(0, values[0].__len__()):
            found_zero = False
            """For each data set"""
            for ii in range(0, width-1):
                if values[ii][jj] == 0.0:
                    found_zero = True
            """If we found no zeros at this gene, then let's keep it."""
            if found_zero == False:
                for ii in range(0, width-1):
                    no_zero_values[ii].append( values[ii][jj] )
        values = no_zero_values
                
    #cranstr += "require(mvtnorm);\n"
    
    pdfpath = filekeyword + ".pdf"
    print "\n. Writing a scatterplot to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=" + (3*width).__str__() + ", height=" + (3*height).__str__() + ");\n"    
    cranstr += "par(mar=c(1.8,2,2.8,1), oma=c(1.5,2,1,1)  );\n"
    colwidth = 1.0 / float(width)
    rowheight = 1.0 / float(height)
    

    """This first loop is to determine total_count"""
    total_count = 0
    for ii in range(0, width):
        mod = 0
        for qq in range(1,width):
            if qq%height == 0:
                if ii >= qq:
                    mod = qq
        for jj in range(mod+(ii%height), mod+height):
            total_count += 1 

    """Now we do the real work"""
    count = 0
    for ii in range(0, width):
        mod = 0
        for qq in range(1,width):
            if qq%height == 0:
                if ii >= qq:
                    mod = qq

        for jj in range(mod+(ii%height), mod+height):       
            count += 1     
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += (ii*colwidth + (0.85*colwidth)).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__()+ "," 
            cranstr += ( (jj%height)*rowheight + 0.86*rowheight).__str__()+ ")"
            if ii > 0 or jj > 0:
                cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
            
            #print (ii*colwidth), ((ii+1)*colwidth), ( (jj%4)*colwidth), (( (jj%4)+1)*colwidth)
            
            if ii in plot_as_rank and jj in plot_as_rank:
                xlim = values[ii].__len__()
                ylim = values[jj].__len__()
                xlim = max(xlim,ylim)
                ylim = xlim
            else:                
                xlim = max( values[ii] )
                ylim = max( values[jj] )
            
            # X values
            found_one = False # did we find at least one valid (x,y) point in the graph that's worth plotting?
            cranstr += "x<-c("
            for xx in range(0, values[ii].__len__()):
                v = values[ii][xx]
                if v != 0 or values[jj][xx] != 0:
                #for v in values_ii:
                    found_one = True
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
            
            """Sanity check: 'x' has more than one value."""
            if found_one == False:
                """Skip to the next ii,jj, iteration."""
                continue
        
            # Y values
            cranstr += "y<-c("
            for yy in range(0, values[jj].__len__()):
                v = values[jj][yy]
                if v != 0 or values[ii][yy] != 0:
                #for v in values_ii:
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
        
            if ii in plot_as_rank and jj in plot_as_rank:
                cranstr += "x <- rank(-x);\n"
                cranstr += "y <- rank(-y);\n"
            
            cranstr += "xlimit <- max(x);\n"
            cranstr += "ylimit <- max(y);\n"
            
            if xlab == None and unit_labels.__len__() > ii:
                xlab = unit_labels[ii]
            else:
                xlab = ""
            
            if ylab == None and unit_labels.__len__() > jj:
                ylab = unit_labels[jj]
            else:
                ylab = ""
            
            cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\", cex=0.5,"
            if force_square:
                cranstr += ", xlim=range(0,max(xlimit,ylimit)), ylim=range(0,max(ylimit,xlimit))"
            else:
                cranstr += ", xlim=range(0,xlimit), ylim=range(0,ylimit)"
            
            col = "black"
            pch = "1"
            if ii < height:
                col = "mediumblue"
                pch = "5"
            elif ii < 2*height:
                col = "red3"
                pch = "5"
            elif ii < 3*height:
                col = "green4"
                pch = "5"
            elif ii < 4*height:
                col = "darkorchid3"
                pch = "5"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ", las=1"
            cranstr += ");\n"
            
            """Write labels across left-side margin"""
            if ii == 0:
                cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=0.8);\n"
            
            
            """Write some summary statistics into the plotting area."""
            
            """Pearson's linear value correlation."""
            corr_valsa = []
            corr_valsb = []
            for ww in range(0, values[ii].__len__()):
                if values[ii][ww] != 0 and values[jj][ww] != 0:
                    corr_valsa.append( values[ii][ww] )
                    corr_valsb.append( values[jj][ww] )
            if corr_valsa.__len__() > 0 and corr_valsb.__len__() > 0:
                (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
                cranstr += "text(" + ((xlim-min(corr_valsa))/2).__str__() + ", " + (0.97*ylim).__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
                
                """Spearman's non-linear non-parametric rank correlation."""
                (rho, pvalue) = scipystats.spearmanr( corr_valsa, corr_valsb )
                cranstr += "text(" + ((xlim-min(corr_valsa))/2).__str__() + ", " + (0.87*ylim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.7);\n"
                
            if force_square:
                cranstr += "abline(0,1)\n"
            
            """Add histograms to this panel."""
            xseq_delta = xlim/float(10.0)
            yseq_delta = ylim/float(10.0)
            cranstr += "hist_a <- hist(x, plot=FALSE, breaks=c( seq(0," + xlim.__str__() + "," + xseq_delta.__str__() + ") ) );\n"#, length.out=20);\n"#, breaks=seq(from=0,to=" + lim.__str__() + "));\n"
            cranstr += "hist_b <- hist(y, plot=FALSE, breaks=c( seq(0," + ylim.__str__() + "," + yseq_delta.__str__() + ") ) );\n"#, length.out=20);\n"#, breaks=seq(from=0,to=" + lim.__str__() + "));\n"
            
            """The top histogram"""
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += (ii*colwidth + (0.85*colwidth)).__str__() + ", "
            cranstr += ( (jj%height)*rowheight + 0.55*rowheight).__str__()+ "," 
            cranstr += ( (jj%height)*rowheight + rowheight).__str__()+ ")"
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
             
            cranstr += "barplot(hist_a$density, axes=FALSE, space=0, col='" + col + "');\n" # barplot

            """Write labels for top row"""
            if (jj+1)%height==0:
                cranstr += "mtext(\"" + names[ii] + "\", side=3, line=1, col=\"black\", cex=0.8);\n"
 
            """The right-side histogram"""
            cranstr += "par( fig=c(" + (ii*colwidth + 0.66*colwidth).__str__() + ","
            cranstr += (ii*colwidth + colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__()+ "," 
            cranstr += ( (jj%height)*rowheight + 0.85*rowheight).__str__()+ ")"
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"

            cranstr += "barplot(hist_b$density, axes=FALSE, space=0, horiz=TRUE, col='" + col + "');\n" # barplot

    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=1.5);\n"
    
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave < " + cranpath)
    return cranpath 


