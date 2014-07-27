import math, os, re, sys

from scipy import stats

def scatter1(values_a, values_b, filekeyword, xlab="", ylab=""):    
    pdfpath = filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=6, height=6);\n"    

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
        
    cranstr += "plot(x, y, xlab=\"" + xlab + "\", ylab=\"" + ylab + "\");\n"
    
    (rho, pvalue) = stats.spearmanr( values_a, values_b )
    cranstr += "text(" + ((max(values_a)-min(values_a))/2).__str__() + ", " + max(values_b).__str__() + ", \"R=%.3f"%rho + ", P=" + pvalue.__str__() + "\");\n"
    
    cranpath = filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)