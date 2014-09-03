import math, os, re, sys

from scipy import stats as scipystats

def scatter1(values_a, values_b, filekeyword, xlab="", ylab="", force_square=False):    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    cranstr += "pdf(\"" + pdfpath + "\", width=6, height=6);\n"    

    print "\n. Writing a scatterplot to", pdfpath

    # X values
    cranstr += "x<-c("
    for v in values_a:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"

    # Y values
    cranstr += "y<-c("
    for v in values_b:
        cranstr += v.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""

    maxx = max(values_a)
    maxy = max(values_b)
    lim = max( [maxx, maxy] )    
    if force_square:
        cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
    
    cranstr += ");\n"
    
    if force_square:
        cranstr += "abline(0,1)\n"
    
    """Pearson's linear value correlation."""
    corr_valsa = []
    corr_valsb = []
    for ww in range(0, values_a.__len__()):
        if values_a[ww] != 0 and values_b[ww] != 0:
            corr_valsa.append( values_a[ww] )
            corr_valsb.append( values_b[ww] )
    (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + lim.__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
    
    """Spearman's non-linear non-parametric rank correlation."""
    (rho, pvalue) = scipystats.spearmanr( corr_valsa, corr_valsb )
    cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + (0.92*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
    
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
        values_a = values[setname][0]
        values_b = values[setname][1]

        cranstr += "par( fig=c(" + (ii*colwidth).__str__() + "," + ((ii+1)*colwidth).__str__() + ", 0, 1)"
        if ii != 0:
            cranstr += ", new=TRUE"
        cranstr += ");\n"

        # X values
        cranstr += "x<-c("
        for v in values_a:
            cranstr += v.__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"
    
        # Y values
        cranstr += "y<-c("
        for v in values_b:
            cranstr += v.__str__() + ","
        cranstr = re.sub(",$", "", cranstr)
        cranstr += ");\n"
        
        cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""

        maxx = max(values_a)
        maxy = max(values_b)
        lim = max( [maxx, maxy] )        
        if force_square:
            cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
        
        cranstr += ");\n"
        
        if force_square:
            cranstr += "abline(0,1)\n"
        
        """Pearson's linear value correlation."""
        corr_valsa = []
        corr_valsb = []
        for ww in range(0, values_a.__len__()):
            if values_a[ww] != 0 and values_b[ww] != 0:
                corr_valsa.append( values_a[ww] )
                corr_valsb.append( values_b[ww] )
        (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
        cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + lim.__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
        
        """Spearman's non-linear non-parametric rank correlation."""
        (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
        cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + (0.92*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\");\n"
    
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
    height, width = number of scatterplots
    values[ii] = list of data. There should be 'width' number of entries in values.
    
    NOTE: in making the scatterplots, cases where both the x-axis and y-axis series are (0,0) will
    be ignored from the plot.
    
    """  
    
    #print "156:", width, height
      
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
            print "\n. ERROR scatter_nxm can't deal with empty value arrays, i.e., you have no data to plot."
            exit()
        
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    print "\n. Computing IDR and plotting to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=" + (9*height).__str__() + ", height=" + (3*height).__str__() + ");\n"    
    cranstr += "par(mar=c(1.8,2.8,2.8,1), oma=c(1.5,2,1,1)  );\n"
    #cranstr += "par(mar=c(4,6,2.8,1), oma=c(1.5,2,1,1)  );\n"
    colwidth = 1.0 / (3.0 * float(height))
    rowheight = 1.0 / float(height)
    gridsize = (height*height - (height-1)*(0.5*height))
    #print "179:", colwidth, rowheight
    
    """This first loop is to determine total_count"""
    total_count = 0
    for ii in range(0, width):
        mod = 0
        for qq in range(1, width):
            if qq%height == 0:
                if ii >= qq:
                    mod = qq
        for jj in range(mod+(ii%height), mod+height):
            total_count += 1 

    """Now we do the real work"""
    count = 0
    for ii in range(0, width):
        mod = 0
        for qq in range(1, width):
            if qq%height == 0:
                if ii >= qq:
                    mod = qq

        for jj in range(mod+(ii%height), mod+height):       
            count += 1
            
            #print "202:", ii, jj, count

            if count > gridsize and (count-1)%gridsize == 0:
                cranstr += "par(mar=c(1.8,2.8,2.8,1), oma=c(1.5,2,1,1)  );\n"
                #cranstr += "par(mar=c(4,6,2.8,1), oma=c(1.5,2,1,1)  );\n"
                cranstr += "plot.new();\n"
            
            if ii == jj:
                continue

                 
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
                                    
            values_a = values[ii]
            values_b = values[jj]
            
            if values_a.__len__() != values_b.__len__():
                print "Error: different sizes for values_a and values_b. Mark point 225."
                exit()
    
            """ "X values"""
            cranstr += "x<-c("
            for xx in range(0, values_a.__len__()):
                v = values_a[xx]
                if v != 0 or values_b[xx] != 0:
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
            
            """Y values"""
            cranstr += "y<-c("
            for yy in range(0, values_b.__len__()):
                v = values_b[yy]
                if v != 0 or values_a[yy] != 0:
                #for v in values_a:
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
            
            cranstr += "rankx <- rank(-x);\n"
            cranstr += "ranky <- rank(-y);\n"
            cranstr += "lim <- max( max(rankx), max(ranky) );\n"
            
            maxa = values_a.__len__()
            maxb = values_b.__len__()
            lim = max( [maxa, maxb] )    
            
            """Plot the ranked data"""
            cranstr += "par( fig=c(" + ((ii%height)*colwidth).__str__() + ","
            cranstr += ((ii%height)*colwidth + colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__() + "," 
            cranstr += ( (jj%height)*rowheight + rowheight).__str__() + ")"
            #print "252:", count, gridsize
            if count > 1:
                cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
            
            #print ii, colwidth
            #print ii*colwidth, (ii*colwidth + colwidth), ( (jj%height)*rowheight), ( (jj%height)*rowheight + rowheight)
            
            cranstr += "plot(rankx, ranky, cex.lab=0.8, xlab=\"rank\", ylab=\"rank\""
            if force_square:
                cranstr += ", xlim=range(0,lim), ylim=range(0,lim)"
            col = "black"
            pch = "1"
            if ii < height:
                col = "mediumblue"
                pch = "3"
            elif ii < 2*height:
                col = "red3"
                pch = "5"
            elif ii < 3*height:
                col = "green4"
                pch = "1"
            elif ii < 4*height:
                col = "darkorchid3"
                pch = "23"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ", las=1"
            #cranstr += ", main=\"" +ii.__str__() + ":" + jj.__str__() + "\""
            cranstr += ");\n"             
                        
            #"""Write labels across left-side margin"""
            #if ii == 0:
            #    cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=1.7);\n"
            #
            if force_square:
                cranstr += "abline(0,1)\n"
            
            """Custom write axis labels"""
            if (jj+1)%height==0:
                cranstr += "mtext(\"" + names[ii] + "\", side=1, line=2, col=\"black\", cex=1);\n"
            if ii == 0:
                cranstr += "mtext(\"" + names[jj] + "\", side=3, line=1, col=\"black\", cex=1);\n"
            
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
            
            #"""Write labels for top row"""
            #if (jj+1)%height==0:
            #    cranstr += "mtext(\"" + names[ii] + " vs. " + names[jj] + "\", side=3, line=1, col=\"black\", cex=1.2);\n"
            
            """Plot the psi values from IDR."""
            cranstr += "par( fig=c(" + (height*colwidth + (ii%height)*colwidth).__str__() + ","
            cranstr += (height*colwidth + (ii%height)*colwidth + colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__() + "," 
            cranstr += ( (jj%height)*rowheight + rowheight).__str__() + ")"
            if ii > 0 or jj > 0:
                cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
            
            cranstr += "plot(uv$psi.n$t, uv$psi.n$value, xlab=\"t\", ylab=\"psi\", xlim=c(0, max(uv$psi.n$t)),ylim=c(0, max(uv$psi.n$value)), cex.lab=0.8"
            col = "black"
            pch = "1"
            if ii < height:
                col = "mediumblue"
                pch = "3"
            elif ii < 2*height:
                col = "red3"
                pch = "5"
            elif ii < 3*height:
                col = "green4"
                pch = "1"
            elif ii < 4*height:
                col = "darkorchid3"
                pch = "23"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ", las=1"
            cranstr += ");\n"
            cranstr += "lines(uv$psi.n$smoothed.line, lwd=4);\n"
            cranstr += "abline(coef=c(0,1), lty=3);\n"
            
            if (jj+1)%height == 0 and ii==0:
                cranstr += "mtext(\"t\", side=1, line=2, col=\"black\", cex=1.2);\n"
                cranstr += "mtext(\"psi\", side=2, line=2.5, col=\"black\", cex=1.2);\n"
                """Write labels for top row"""
                cranstr += "mtext(\"IDR Correspondence Curve\", side=3, line=1, col=\"black\", cex=1.2);\n"
            
            """Plot the psi-prime values from IDR"""
            cranstr += "par( fig=c(" + ( (2*height*colwidth) + (ii%height)*colwidth).__str__() + ","
            cranstr += ( (2*height*colwidth) + (ii%height)*colwidth + colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__() + "," 
            cranstr += ( (jj%height)*rowheight + rowheight).__str__() + ")"
            if ii > 0 or jj > 0:
                cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
             
            cranstr += "plot(uv$dpsi.n$t, uv$dpsi.n$value, xlab=\"t\", ylab=\"psi'\", xlim=c(0, max(uv$dpsi.n$t)),ylim=c(0, max(uv$dpsi.n$value)), cex.lab=1.0"
            col = "black"
            pch = "1"
            if ii < height:
                col = "mediumblue"
                pch = "3"
            elif ii < 2*height:
                col = "red3"
                pch = "5"
            elif ii < 3*height:
                col = "green4"
                pch = "1"
            elif ii < 4*height:
                col = "darkorchid3"
                pch = "23"
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
            
            """Write the IDR data to a text file, so that Python can import the data."""
            tablepath = filekeyword + ".xls"
            #cranstr += "write.table(idr.out, \"" + tablepath + "\", sep=\"\t\");\n"
            
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=2.2);\n"
    
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave < " + cranpath)
    return cranpath 

def scatter_nxm(width, height, values, names, filekeyword, title="", xlab="", ylab="", force_square=True):    
    """
    height, width = number of scatterplots
    values[ii] = list of data. There should be 'width' number of entries in values.
    
    NOTE: in making the scatterplots, cases where both the x-axis and y-axis series are (0,0) will
    be ignored from the plot.
    
    """    
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
            print "\n. ERROR scatter_nxm can't deal with empty value arrays, i.e., you have no data to plot."
            exit()
    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
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
            
            values_a = values[ii]
            values_b = values[jj]
    
            # X values
            cranstr += "x<-c("
            for xx in range(0, values_a.__len__()):
                v = values_a[xx]
                if v != 0 or values_b[xx] != 0:
                #for v in values_a:
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
        
            # Y values
            cranstr += "y<-c("
            for yy in range(0, values_b.__len__()):
                v = values_b[yy]
                if v != 0 or values_a[yy] != 0:
                #for v in values_a:
                    cranstr += v.__str__() + ","
            cranstr = re.sub(",$", "", cranstr)
            cranstr += ");\n"
            
            cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\""
            
            maxa = max(values_a)
            maxb = max(values_b)
            
            lim = max( [maxa, maxb] )
            if force_square:
                cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
            
            col = "black"
            pch = "1"
            if ii < height:
                col = "mediumblue"
                pch = "3"
            elif ii < 2*height:
                col = "red3"
                pch = "5"
            elif ii < 3*height:
                col = "green4"
                pch = "1"
            elif ii < 4*height:
                col = "darkorchid3"
                pch = "23"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ", las=1"
            cranstr += ");\n"
            
            """Write labels across left-side margin"""
            if ii == 0:
                cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=1.7);\n"
            
            
            """Write some summary statistics into the plotting area."""
            
            """Pearson's linear value correlation."""
            corr_valsa = []
            corr_valsb = []
            for ww in range(0, values_a.__len__()):
                if values_a[ww] != 0 and values_b[ww] != 0:
                    corr_valsa.append( values_a[ww] )
                    corr_valsb.append( values_b[ww] )
            (rho, pvalue) = scipystats.pearsonr( corr_valsa, corr_valsb )
            cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + (0.97*lim).__str__() + ", \"Prs R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.9);\n"
            
            """Spearman's non-linear non-parametric rank correlation."""
            (rho, pvalue) = scipystats.spearmanr( corr_valsa, corr_valsb )
            cranstr += "text(" + ((lim-min(values_a))/2).__str__() + ", " + (0.87*lim).__str__() + ", \"Spr R=%.2f"%rho + ", P=%.2f"%pvalue + "\", cex=0.9);\n"
            #cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.95*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\", cex=1.2);\n"
            #cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.85*max( [maxa, maxb] )).__str__() + ", \"" + ii.__str__() + "," + jj.__str__() + "\");\n"
            
            if force_square:
                cranstr += "abline(0,1)\n"
            
            """Add histograms to this panel."""
            seq_delta = lim/float(10.0)
            cranstr += "hist_a <- hist(x, plot=FALSE, breaks=c( seq(0," + lim.__str__() + "," + seq_delta.__str__() + ") ) );\n"#, length.out=20);\n"#, breaks=seq(from=0,to=" + lim.__str__() + "));\n"
            cranstr += "hist_b <- hist(y, plot=FALSE, breaks=c( seq(0," + lim.__str__() + "," + seq_delta.__str__() + ") ) );\n"#, length.out=20);\n"#, breaks=seq(from=0,to=" + lim.__str__() + "));\n"
            
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
                cranstr += "mtext(\"" + names[ii] + "\", side=3, line=1, col=\"black\", cex=1.7);\n"
 
            """The right-side histogram"""
            cranstr += "par( fig=c(" + (ii*colwidth + 0.66*colwidth).__str__() + ","
            cranstr += (ii*colwidth + colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__()+ "," 
            cranstr += ( (jj%height)*rowheight + 0.85*rowheight).__str__()+ ")"
            cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"

            cranstr += "barplot(hist_b$density, axes=FALSE, space=0, horiz=TRUE, col='" + col + "');\n" # barplot

    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=2.2);\n"
    
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave < " + cranpath)
    return cranpath 


