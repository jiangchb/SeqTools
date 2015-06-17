import re, os, sys, time
from annotation_db import *

def write_html_results(con, vizdbpath):
    """Creates an HTML results page that provides URLs to all
        the files created by the pipeline."""
    vizcon = lite.connect(vizdbpath, timeout=1)
    
    html = [] # a list of HTML lines that will be printed to the output file.    
    html.append("<!DOCTYPE html>")
    html.append("<html>")
    html.append("<head>")
    html.append("<title>" + get_setting("project_name", con) + "</title>")
    html.append("</head>")
    html.append("<body>")
    
    html.append("<table>")
    html.append("<tr>")
    html.append("<th>Replicate Name</th>")
    html.append("<th>Species</th>")
    html.append("<th>Reads</th>")
    html.append("<th>Fold Enrichment</th>")
    html.append("<th>Peaks</th>")    
    html.append("</tr>")
    vcur = vizcon.cursor()
    sql = "select id, name, species from Replicates"
    vcur.execute(sql)
    x = vcur.fetchall()
    for ii in x:
        html.append("<tr>")
        html.append("<td>" + ii[1].__str__() + "</td>")
        html.append("<td>" + ii[2].__str__() + "</td>")
        html.append("<td>fastq</td>")
        html.append("<td></td>")
        html.append("<td>fastq</td>")
        html.append("</tr>")
    html.append("</table>")

    html.append("</body>")
    html.append("</html>")
    
    fout = open("results." + get_setting("project_name", con).__str__() + ".html", "w")
    for hl in html:
        fout.write(hl + "\n")
    fout.close()