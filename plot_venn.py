        
def plot_venn_diagram(data, filekeyword):
    """data[key] = list of values in key's set."""
    """Note that this method will only print a Venn diagram for up to the first 3 series in the data.
    In other words, this method cannot make 4-way Venn diagrams."""
    from matplotlib import pyplot as plt
    import numpy as np
    
    pdfpath = filekeyword + ".venn.pdf"
    
    # 3-way Venn:
    if data.keys().__len__() >= 3:
        k = data.keys()
        from matplotlib_venn import venn3, venn3_circles
        f = plt.figure(figsize=(4,4))
        set0 = set( data[k[0]] )
        set1 = set( data[k[1]] )
        set2 = set( data[k[2]] )
        venn3( [set0, set1, set2], (k[0], k[1], k[2]) )
        print "\n. Writing a 3-way Venn diagram to", pdfpath
        f.savefig(pdfpath, format="pdf")

    # 2-way Venn:
    elif data.keys().__len__() == 2:
        k = data.keys()
        from matplotlib_venn import venn2, venn2_circles
        f = plt.figure(figsize=(4,4))
        set0 = set( data[k[0]] )
        set1 = set( data[k[1]] )
        venn2( [set0, set1], (k[0], k[1]) )
        print "\n. Writing a 2-way Venn diagram to", pdfpath
        f.savefig(pdfpath, format="pdf")

    return pdfpath