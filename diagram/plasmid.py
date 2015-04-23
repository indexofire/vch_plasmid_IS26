from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SeqIO.read("../sequences/0001.gb", "genbank")

gd_diagram = GenomeDiagram.Diagram("Vibrio cholerae  plasmid 2012EL-2176")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
   # print type(feature.qualifiers.get('product'))
    if feature.type != "gene":
        #print type(feature.qualifiers.get('product'))
        if isinstance(feature.qualifiers.get('product'), list) and feature.qualifiers.get('product')[0] == 'transposase':
            #print feature.qualifiers.get('product')
            if len(gd_feature_set) % 2 == 0:
                color = colors.blue
            else:
                color = colors.lightblue
            gd_feature_set.add_feature(feature, color=color, label=True)

gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=1, start=0, end=len(record))
gd_diagram.write("linear.pdf", "PDF")
#gd_diagram.write("plasmid_linear.eps", "EPS")
#gd_diagram.write("plasmid_linear.svg", "SVG")
