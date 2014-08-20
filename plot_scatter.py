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


def scatter_nxm(width, height, values, names, filekeyword, title="", xlab="", ylab="", force_square=True):    
    """
    height, width = number of scatterplots
    values[ii] = list of data. There should be 'width' number of entries in values.
    
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
            cranstr += ((ii+1)*colwidth).__str__() + ", "
            cranstr += ( (jj%height)*rowheight).__str__()+ "," 
            cranstr += (( (jj%height)+1)*rowheight).__str__()+ ")"
            if ii > 0 or jj > 0:
                cranstr += ", new=TRUE" # don't call new=TRUE for the first plot.
            cranstr += ");\n"
            
            #print (ii*colwidth), ((ii+1)*colwidth), ( (jj%4)*colwidth), (( (jj%4)+1)*colwidth)
            
            values_a = values[ii]
            values_b = values[jj]
    
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
            cranstr += ");\n"
            
            if ii == 0: #left most column
                cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=1.7);\n"
            if (jj+1)%height==0:
                cranstr += "mtext(\"" + names[ii] + "\", side=3, line=1, col=\"black\", cex=1.7);\n"
            
            """Write some summary statistics into the plotting area."""
            
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
            #cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.95*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\", cex=1.2);\n"
            #cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.85*max( [maxa, maxb] )).__str__() + ", \"" + ii.__str__() + "," + jj.__str__() + "\");\n"
            
            if force_square:
                cranstr += "abline(0,1)\n"
                
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=2.2);\n"
    
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave < " + cranpath)
    return cranpath 


