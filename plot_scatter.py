import math, os, re, sys

from scipy import stats as scipystats

def scatter1(values_a, values_b, filekeyword, xlab="", ylab="", force_square=False):    
    pdfpath = filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=6, height=6);\n"    

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
    
    if force_square:
        maxx = max(values_a)
        maxy = max(values_b)
        lim = max( [maxx, maxy] )
        cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
    
    cranstr += ");\n"
    
    if force_square:
        cranstr += "abline(0,1)\n"
    
    (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
    cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + max(values_b).__str__() + ", \"R=%.3f"%rho + ", P=" + pvalue.__str__() + "\");\n"
    
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)
    
    return cranpath

def scatter1xn(values, filekeyword, xlab="", ylab="", force_square=False):
    """Makes a multi-panel scatterplot, with 1 row and N columns."""
    
    """values[set name][0] = X value list, values[set name][1] = y value list."""
    pdfpath = filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=9, height=3);\n"    

    print "\n. Writing a scatterplot to", pdfpath

    cranstr += "par(mar=c(1.8,1.8,0.4,0.4), oma=c(1.5,2,1,1)  );\n"

    colwidth = 1.0 / values.__len__()

    sets = values.keys()
    for ii in range(0, sets.__len__()):
        setname = sets[ii]
        values_a = values[setname][0]
        values_b = values[setname][1]

        cranstr += "par( fig=c(" + (ii*colwidth).__str__() + "," + ((ii+1)*colwidth).__str__() + ", 0, 1), new=TRUE);\n"

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
        
        if force_square:
            maxx = max(values_a)
            maxy = max(values_b)
            lim = max( [maxx, maxy] )
            cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
        
        cranstr += ");\n"
        
        if force_square:
            cranstr += "abline(0,1)\n"
        
        (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
        cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + max(values_b).__str__() + ", \"R=%.3f"%rho + ", P=" + pvalue.__str__() + "\");\n"
    
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)
    
    return cranpath

def scatter4x4(values, names, filekeyword, title="", xlab="", ylab="", force_square=True):
    """Values[ii] = list of data. There should be 4 sets.
    """    
    pdfpath = filekeyword + ".pdf"
    print "\n. Writing a scatterplot to", pdfpath
    cranstr = "pdf(\"" + pdfpath + "\", width=12, height=12);\n"    
    cranstr += "par(mar=c(1.8,1.8,1.0,0.4), oma=c(1.5,2,1,1)  );\n"
    colwidth = 0.25
    for ii in range(0, 4):
        for jj in range(ii, 4):
            print ii, jj
            
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += ((ii+1)*colwidth).__str__() + ", "
            cranstr += (jj*colwidth).__str__()+ "," 
            cranstr += ((jj+1)*colwidth).__str__()+ "), new=TRUE);\n"
            
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
            
            if force_square:
                lim = max( [maxb, maxb] )
                cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
            
            cranstr += ");\n"
            
            (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
            cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.9*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\");\n"
            
            
            if force_square:
                cranstr += "abline(0,1)\n"
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=0);\n"
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save < " + cranpath)
    
    