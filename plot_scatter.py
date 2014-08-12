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
        
        if force_square:
            maxx = max(values_a)
            maxy = max(values_b)
            lim = max( [maxx, maxy] )
            cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
        
        cranstr += ");\n"
        
        if force_square:
            cranstr += "abline(0,1)\n"
        
        (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
        cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + max(values_b).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\");\n"
    
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=-0.8, cex=1.8);\n"
    cranstr += "dev.off();\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save --slave < " + cranpath)
    
    return cranpath

def scatter4x4(values, names, filekeyword, title="", xlab="", ylab="", force_square=True):
    """Values[ii] = list of data. There should be 4 sets.
    """ 
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
       
    pdfpath = filekeyword + ".pdf"
    print "\n. Writing a scatterplot to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=12, height=12);\n"    
    cranstr += "par(mar=c(1.8,1.8,1.0,0.4), oma=c(1.5,2,1,1)  );\n"
    colwidth = 0.25
    for ii in range(0, 4):
        for jj in range(ii, 4):
            print ii, jj
            
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += ((ii+1)*colwidth).__str__() + ", "
            cranstr += (jj*colwidth).__str__()+ "," 
            cranstr += ((jj+1)*colwidth).__str__()+ ")"
            if ii != 0 and jj != 0:
                cranstr += ", new=TRUE"
            cranstr += ");\n"
            
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
                lim = max( [maxa, maxb] )
                cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
            
            cranstr += ");\n"
            
            (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
            cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.9*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\");\n"
            cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.7*max( [maxa, maxb] )).__str__() + ", \"" + ii.__str__() + "," + jj.__str__() + "\");\n"
            
            if force_square:
                cranstr += "abline(0,1)\n"
    cranstr += "mtext(\"" + title + "\", side=3, outer=TRUE, line=0);\n"
    cranstr += "dev.off();\n"
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save --slave < " + cranpath)

def scatter8x4(values, names, filekeyword, title="", xlab="", ylab="", force_square=True):
    """Values[ii] = list of data. There should be 4 sets.
    """    
    if names.__len__() != 8:
        print "\n. ERROR plot_scatter.py 193, you called scatter12x4 without enough names."
        print names
        print names.__len__()
        exit()
    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    print "\n. Writing a scatterplot to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=24, height=12);\n"    
    cranstr += "par(mar=c(1.8,2,2.8,1), oma=c(1.5,2,1,1)  );\n"
    colwidth = 0.125
    
    for ii in range(0, 8):
        mod = 0
        if ii > 3:
            mod = 4
        for jj in range(mod+(ii%4), mod+4):
            sys.stdout.write(".")
            sys.stdout.flush()
            
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += ((ii+1)*colwidth).__str__() + ", "
            cranstr += ( (jj%4)*0.25).__str__()+ "," 
            cranstr += (( (jj%4)+1)*0.25).__str__()+ ")"
            if ii != 0 and jj != 0:
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
            
            if force_square:
                lim = max( [maxa, maxb] )
                cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
            
            col = "black"
            pch = "1"
            if ii < 4:
                col = "mediumblue"
                pch = "3"
            else:
                col = "red3"
                pch = "1"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ");\n"
            
            if ii == 0: #left most column
                cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=1.7);\n"
            if jj == 3 or jj == 7:
                cranstr += "mtext(\"" + names[ii] + "\", side=3, line=1, col=\"black\", cex=1.7);\n"
            
            """Write some summary statistics into the plotting area."""
            (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
            cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.95*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\", cex=1.2);\n"
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

def scatter12x4(values, names, filekeyword, title="", xlab="", ylab="", force_square=True):
    """Values[ii] = list of data. There should be 4 sets.
    """    
    if names.__len__() != 12:
        print "\n. ERROR plot_scatter.py 177, you called scatter12x4 without enough names."
        print names
        print names.__len__()
        exit()
    
    sinkpath = filekeyword + ".out"
    cranstr = "sink(\"" + sinkpath + "\", append=FALSE, split=FALSE);\n"
    
    pdfpath = filekeyword + ".pdf"
    print "\n. Writing a scatterplot to", pdfpath
    cranstr += "pdf(\"" + pdfpath + "\", width=36, height=12);\n"    
    cranstr += "par(mar=c(1.8,2,2.8,1), oma=c(1.5,2,1,1)  );\n"
    colwidth = 0.08333
    
    for ii in range(0, 12):
        mod = 0
        if ii > 3:
            mod = 4
        if ii > 7:
            mod = 8
        for jj in range(mod+(ii%4), mod+4):
            sys.stdout.write(".")
            sys.stdout.flush()
            
            cranstr += "par( fig=c(" + (ii*colwidth).__str__() + ","
            cranstr += ((ii+1)*colwidth).__str__() + ", "
            cranstr += ( (jj%4)*0.25).__str__()+ "," 
            cranstr += (( (jj%4)+1)*0.25).__str__()+ ")"
            if ii != 0 and jj != 0:
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
            
            if force_square:
                lim = max( [maxa, maxb] )
                cranstr += ", xlim=range(0," + lim.__str__() + "), ylim=range(0," + lim.__str__() + ")"
            
            col = "black"
            pch = "1"
            if ii < 4:
                col = "mediumblue"
                pch = "3"
            elif ii < 8:
                col = "green4"
                pch = "5"
            else:
                col = "red3"
                pch = "1"
            cranstr += ", col=\"" + col + "\""
            cranstr += ", pch=" + pch
            cranstr += ");\n"
            
            if ii == 0: #left most column
                cranstr += "mtext(\"" + names[jj] + "\", side=2, line=2, col=\"black\", cex=1.7);\n"
            if jj == 3 or jj == 7 or jj == 11:
                cranstr += "mtext(\"" + names[ii] + "\", side=3, line=1, col=\"black\", cex=1.7);\n"
            
            """Write some summary statistics into the plotting area."""
            (rho, pvalue) = scipystats.spearmanr( values_a, values_b )
            cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + (0.95*max( [maxa, maxb] )).__str__() + ", \"R=%.3f"%rho + ", P=%.3f"%pvalue + "\", cex=1.2);\n"
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
    